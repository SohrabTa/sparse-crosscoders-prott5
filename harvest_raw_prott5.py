"""
Harvest raw ProtT5 per-residue embeddings for ProteinGym matched assays.

For each assay in --matches:
  1. Reconstruct WT sequence from the first variant.
  2. Run ProtT5 once on WT, store per-residue hidden states from --layer.
  3. Run ProtT5 on each variant (batched), store per-residue hidden states.
  4. Save a per-assay .npz with:
       mutant[N]                   : variant identifier string
       position[N]                 : centroid of mutated positions (1-based)
       n_mutations[N]              : number of mutated positions (≥1)
       DMS_score[N]                : experimental fitness score
       emb_wt_mean[D]              : mean-pooled WT embedding (one per assay)
       emb_wt_at_pos[N, D]         : WT embedding at the variant's mutated pos
       emb_mut_mean[N, D]          : mean-pooled mutant embedding
       emb_mut_at_pos[N, D]        : mutant embedding at the mutated pos
       layer[scalar]               : which ProtT5 encoder block (0-indexed)

Powers the H-T2 ceiling test (Adams 2025-style ridge probe): does raw ProtT5
encode fitness at this layer, independent of any SAE/crosscoder?

The script is intentionally crosscoder-free — it should run in any env with
torch + transformers + sentencepiece + pyyaml. No InterPLM or crosscode
imports required.

Usage:
    uv run --project repos/crosscode python repos/sparse-crosscoders-prott5/harvest_raw_prott5.py \\
        --matches data/proteingym_concept_matches.csv \\
        --dms_dir data/DMS_ProteinGym_substitutions \\
        --output_dir data/proteingym/raw_prott5_layer12 \\
        --layer 12 --batch_size 8
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
from transformers import T5EncoderModel, T5Tokenizer

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("harvest")


MUTANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z\*])$")
PROTT5_NAME = "Rostlab/prot_t5_xl_uniref50"
N_LAYERS = 24
D_MODEL = 1024


def parse_mutant_positions(mutant_field: str) -> List[Tuple[str, int, str]]:
    out: List[Tuple[str, int, str]] = []
    for token in str(mutant_field).split(":"):
        m = MUTANT_RE.match(token.strip())
        if not m:
            return []
        out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


class ProtT5SingleLayerHarvester:
    """Loads ProtT5 once and exposes a single hook on the chosen encoder block.

    Returns per-residue hidden states for one layer (not all 24 like the
    InterPLM embedder), so memory is 24× smaller and we can run larger batches.
    """

    def __init__(self, layer: int, device: torch.device, dtype: torch.dtype):
        self.layer = layer
        self.device = device
        self.dtype = dtype
        log.info("Loading ProtT5 (%s) on %s, dtype=%s",
                 PROTT5_NAME, device, dtype)
        self.tokenizer = T5Tokenizer.from_pretrained(
            PROTT5_NAME, do_lower_case=False, legacy=True
        )
        self.model = T5EncoderModel.from_pretrained(PROTT5_NAME, torch_dtype=dtype)
        self.model = self.model.to(device).eval()

    @torch.no_grad()
    def encode(self, sequences: List[str]) -> List[np.ndarray]:
        """Run ProtT5 on a list of sequences, return per-sequence (L, D) arrays
        at the configured layer. Strips the trailing EOS token and any padding."""
        # ProtT5 wants whitespace-separated AA tokens with UZOB→X
        processed = [" ".join(re.sub(r"[UZOB]", "X", s)) for s in sequences]
        inputs = self.tokenizer(
            processed, add_special_tokens=True, padding="longest", return_tensors="pt"
        )
        input_ids = inputs["input_ids"].to(self.device)
        attention_mask = inputs["attention_mask"].to(self.device)

        # Single hook on the chosen encoder block
        cache: Dict[str, torch.Tensor] = {}

        def hook(module, inp, out):
            cache["x"] = (out[0] if isinstance(out, tuple) else out).detach()

        h = self.model.encoder.block[self.layer].register_forward_hook(hook)
        try:
            self.model(input_ids=input_ids, attention_mask=attention_mask)
        finally:
            h.remove()

        layer_out = cache["x"]  # (B, T, D) — includes padding + EOS
        # Strip padding & EOS by using each sequence's actual AA length
        outputs: List[np.ndarray] = []
        for i, seq in enumerate(sequences):
            n = len(seq)
            outputs.append(layer_out[i, :n].float().cpu().numpy())
        return outputs


def harvest_one_assay(
    dms_id: str,
    dms_path: Path,
    harvester: ProtT5SingleLayerHarvester,
    output_path: Path,
    batch_size: int,
    max_variants: Optional[int] = None,
) -> Optional[dict]:
    """Returns timing / size info dict, or None if nothing scored."""
    if not dms_path.exists():
        log.error("  DMS file missing: %s", dms_path); return None
    dms = pd.read_csv(dms_path)
    if max_variants and len(dms) > max_variants:
        dms = dms.sample(max_variants, random_state=42).reset_index(drop=True)
    if len(dms) == 0:
        return None

    # WT reconstruction (same logic as score_proteingym_features.py)
    first_mut = parse_mutant_positions(dms["mutant"].iloc[0])
    if not first_mut:
        log.error("  could not parse mutant '%s'", dms["mutant"].iloc[0])
        return None
    wt_seq = list(dms["mutated_sequence"].iloc[0])
    for waa, p, _ in first_mut:
        wt_seq[p - 1] = waa
    wt_seq_str = "".join(wt_seq)
    L_wt = len(wt_seq_str)
    D = D_MODEL

    # Encode WT once
    t0 = time.time()
    wt_emb = harvester.encode([wt_seq_str])[0]  # (L_wt, D)
    emb_wt_mean = wt_emb.mean(axis=0).astype(np.float32)  # (D,)
    t_wt = time.time() - t0

    # Pre-parse all variants
    parsed_per_row: List[List[Tuple[str, int, str]]] = []
    dms_scores: List[float] = []
    mutants: List[str] = []
    for _, r in dms.iterrows():
        ps = parse_mutant_positions(r["mutant"])
        if not ps:
            continue
        parsed_per_row.append(ps)
        dms_scores.append(float(r["DMS_score"]))
        mutants.append(str(r["mutant"]))
    if not parsed_per_row:
        log.warning("  no parsable variants"); return None
    N = len(parsed_per_row)

    # Allocate output arrays
    centroid_pos = np.array(
        [int(round(np.mean([p for _, p, _ in ps]))) for ps in parsed_per_row],
        dtype=np.int32,
    )
    n_mut = np.array([len(ps) for ps in parsed_per_row], dtype=np.int16)
    emb_wt_at_pos = np.zeros((N, D), dtype=np.float32)
    emb_mut_mean = np.zeros((N, D), dtype=np.float32)
    emb_mut_at_pos = np.zeros((N, D), dtype=np.float32)
    for i, ps in enumerate(parsed_per_row):
        # average WT embedding across this variant's mutated positions
        positions = [p - 1 for _, p, _ in ps if 0 <= p - 1 < L_wt]
        if positions:
            emb_wt_at_pos[i] = wt_emb[positions].mean(axis=0)

    # Process variants in batches
    t_var0 = time.time()
    seqs = dms["mutated_sequence"].tolist()
    n_batches = (N + batch_size - 1) // batch_size
    for bi in range(n_batches):
        s, e = bi * batch_size, min((bi + 1) * batch_size, N)
        batch_seqs = seqs[s:e]
        embs = harvester.encode(batch_seqs)  # list of (L_i, D)
        for j, emb in enumerate(embs):
            idx = s + j
            emb_mut_mean[idx] = emb.mean(axis=0)
            positions = [p - 1 for _, p, _ in parsed_per_row[idx] if 0 <= p - 1 < emb.shape[0]]
            if positions:
                emb_mut_at_pos[idx] = emb[positions].mean(axis=0)
        if bi % 50 == 0:
            log.info("    batch %d/%d", bi + 1, n_batches)
    t_var = time.time() - t_var0

    np.savez_compressed(
        output_path,
        mutant=np.array(mutants),
        position=centroid_pos,
        n_mutations=n_mut,
        DMS_score=np.array(dms_scores, dtype=np.float32),
        emb_wt_mean=emb_wt_mean,
        emb_wt_at_pos=emb_wt_at_pos,
        emb_mut_mean=emb_mut_mean,
        emb_mut_at_pos=emb_mut_at_pos,
        layer=np.array(harvester.layer, dtype=np.int16),
        seq_len_wt=np.array(L_wt, dtype=np.int32),
    )
    return {"N": N, "L_wt": L_wt, "t_wt": t_wt, "t_var": t_var}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matches", type=Path, required=True)
    ap.add_argument("--dms_dir", type=Path, required=True)
    ap.add_argument("--output_dir", type=Path, required=True)
    ap.add_argument("--layer", type=int, default=12,
                    help="0-indexed ProtT5 encoder block (default 12 ≈ mid-encoder)")
    ap.add_argument("--batch_size", type=int, default=8)
    ap.add_argument("--max_variants", type=int, default=None,
                    help="Cap variants per assay (smoke test)")
    ap.add_argument("--max_seq_len", type=int, default=2048,
                    help="Skip assays whose target_seq_len exceeds this")
    ap.add_argument("--assays", type=Path, default=None,
                    help="Optional TSV/CSV with DMS_id column to restrict the run")
    ap.add_argument("--device", type=str, default=None,
                    help="cuda / cpu / mps; auto-detect if omitted")
    ap.add_argument("--dtype", type=str, default="auto",
                    choices=["auto", "float16", "bfloat16", "float32"])
    ap.add_argument("--skip_existing", action="store_true")
    args = ap.parse_args()

    assert 0 <= args.layer < N_LAYERS, f"--layer must be in [0, {N_LAYERS})"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    matches = pd.read_csv(args.matches)
    if args.max_seq_len:
        n_before = matches["DMS_id"].nunique()
        matches = matches[matches["target_seq_len"] <= args.max_seq_len].copy()
        n_after = matches["DMS_id"].nunique()
        if n_after < n_before:
            log.info("Length filter dropped %d assays (target_seq_len > %d)",
                     n_before - n_after, args.max_seq_len)

    if args.assays:
        sep = "\t" if args.assays.suffix == ".tsv" else ","
        wanted = set(pd.read_csv(args.assays, sep=sep)["DMS_id"].tolist())
        matches = matches[matches["DMS_id"].isin(wanted)].copy()

    assay_list = list(matches["DMS_id"].drop_duplicates())
    log.info("Will harvest %d assays at layer %d", len(assay_list), args.layer)

    # Device + dtype
    if args.device:
        device = torch.device(args.device)
    else:
        if torch.cuda.is_available(): device = torch.device("cuda")
        elif torch.backends.mps.is_available(): device = torch.device("mps")
        else: device = torch.device("cpu")

    if args.dtype == "auto":
        # fp16 on CUDA, fp32 elsewhere (MPS fp16 historically flaky on T5)
        dtype = torch.float16 if device.type == "cuda" else torch.float32
    else:
        dtype = getattr(torch, args.dtype)
    log.info("Device: %s, dtype: %s", device, dtype)

    harvester = ProtT5SingleLayerHarvester(args.layer, device, dtype)

    total_t0 = time.time()
    total_variants = 0
    for i, dms_id in enumerate(assay_list):
        out_path = args.output_dir / f"{dms_id}.npz"
        if args.skip_existing and out_path.exists():
            log.info("[%d/%d] %s — skip (exists)", i + 1, len(assay_list), dms_id)
            continue
        log.info("[%d/%d] %s", i + 1, len(assay_list), dms_id)
        dms_path = args.dms_dir / f"{dms_id}.csv"
        info = harvest_one_assay(
            dms_id, dms_path, harvester, out_path, args.batch_size, args.max_variants
        )
        if info:
            total_variants += info["N"]
            log.info(
                "  N=%d L_wt=%d  WT=%.1fs  variants=%.1fs (%.3f s/var)",
                info["N"], info["L_wt"], info["t_wt"], info["t_var"],
                info["t_var"] / max(1, info["N"]),
            )

    total_t = time.time() - total_t0
    log.info("Done. %d variants in %.1f min (%.3f s/var)",
             total_variants, total_t / 60, total_t / max(1, total_variants))


if __name__ == "__main__":
    sys.exit(main())
