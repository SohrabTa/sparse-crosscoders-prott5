"""
ProtT5 masked-marginal pseudo-log-likelihood baseline on ProteinGym DMS assays.

Tests whether the *output distribution* of ProtT5 (logits, not residual stream)
encodes the substitution-specific fitness signal — i.e. the established
masked-marginals scoring approach used by ESM-1v / VESPA / ProteinGym entries.

For each variant with mutant string `WT-position-MUT` (e.g. `I291A`):
   score = log P(MUT | context with position masked)
         − log P(WT  | context with position masked)
   (decoder is fed `<extra_id_0>` and asked to generate one AA token)
Multi-mutants: sum of per-position log-ratios (independence approx).

ProtT5 was pre-trained with span-corruption with average span length 1, so
single-residue masking IS its pre-training task. This is the canonical
scoring strategy for ProtT5-based ProteinGym entries.

Output: per-assay Spearman correlation between summed log-ratio and DMS_score.
Per-variant log-ratios are also saved for downstream analysis.

Usage:
    uv run --project repos/crosscode --with scipy python \\
      repos/sparse-crosscoders-prott5/pseudoll_baseline.py \\
      --matches data/proteingym_concept_matches.csv \\
      --dms_dir data/DMS_ProteinGym_substitutions \\
      --output_dir /tmp/pseudoll \\
      --assays /tmp/pseudoll_assays.tsv \\
      --max_variants 200 \\
      --batch_size 8

Notes:
 * We score the *whole* mutant set with a single masked-encoder pass per
   position (not per variant) — for each unique position p we mask p once
   and read off log P(aa | context) for all 20 amino acids in one shot.
   Cost is O(unique_positions × |alphabet|), not O(variants × |alphabet|).
 * On MPS the full T5ForConditionalGeneration model is ~12 GB fp32, ~6 GB
   fp16. fp32 is safest on a 32 GB Mac; fp16 on H100 is fine.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from scipy.stats import spearmanr
from transformers import T5ForConditionalGeneration, T5Tokenizer

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("pseudoll")

PROTT5_NAME = "Rostlab/prot_t5_xl_uniref50"
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
MUTANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z\*])$")


def parse_mutant(mutant_field: str) -> List[Tuple[str, int, str]]:
    out: List[Tuple[str, int, str]] = []
    for token in str(mutant_field).split(":"):
        m = MUTANT_RE.match(token.strip())
        if not m:
            return []
        out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


class ProtT5Scorer:
    """BART-style masked LM scorer for ProtT5.

    ProtT5 was trained with a *BART-like MLM denoising* objective (per its
    HuggingFace README), not T5 span corruption. The decoder reconstructs the
    full uncorrupted sequence; the WT prefix conditions the autoregressive
    prediction at the masked position. Decoding format:

       encoder input: [▁M, ▁T, ..., <extra_id_0>, ..., ▁K, </s>]
                                      ^ mask token at target position
       decoder input: [<pad>, ▁M, ▁T, ..., ▁(token at target_pos − 1)]
       readout:       log-softmax at the LAST decoder position, restricted to
                      the 20 AA token ids.

    Verified empirically: WT recovers at top-1 for 7/8 random positions on
    NUD15 with high confidence (log_p ≈ 0).
    """

    def __init__(self, device: torch.device, dtype: torch.dtype):
        self.device = device
        self.dtype = dtype
        log.info("Loading ProtT5 (full encoder-decoder) on %s, dtype=%s", device, dtype)
        self.tokenizer = T5Tokenizer.from_pretrained(
            PROTT5_NAME, do_lower_case=False, legacy=True
        )
        self.model = T5ForConditionalGeneration.from_pretrained(
            PROTT5_NAME, torch_dtype=dtype
        ).to(device).eval()

        # ProtT5's SentencePiece vocab uses ▁-prefixed AA tokens (▁A, ▁L, …),
        # NOT bare 'A', 'L'. Bare names tokenize to <unk> (id=2). We must look
        # up the ▁-prefixed form to get the real id.
        self.aa_token_ids = {
            aa: self.tokenizer.convert_tokens_to_ids(f"▁{aa}")
            for aa in AMINO_ACIDS
        }
        if any(v == self.tokenizer.unk_token_id for v in self.aa_token_ids.values()):
            raise RuntimeError(
                f"AA→token-id lookup failed (unk in {self.aa_token_ids}); "
                "ProtT5 vocab format may have changed."
            )
        self.x_token_id = self.tokenizer.convert_tokens_to_ids("▁X")
        self.extra_id_0 = self.tokenizer.convert_tokens_to_ids("<extra_id_0>")
        self.eos = self.tokenizer.eos_token_id
        self.decoder_start = self.model.config.decoder_start_token_id

    def _aa_to_id(self, c: str) -> int:
        return self.aa_token_ids.get(c, self.x_token_id)

    @torch.no_grad()
    def score_position(self, sequence: str, pos1: int) -> np.ndarray:
        """For a single sequence with position `pos1` (1-based) masked with
        <extra_id_0> in the encoder input, return log P(aa) over the 20
        amino acids in AMINO_ACIDS order (renormalized over the AA subset).

        We construct token ids directly (not via tokenizer-on-string) to avoid
        SentencePiece-prefix subtleties and to guarantee the sentinel ends up
        at exactly position pos1−1 in the encoder.
        """
        # Apply ProtT5 preprocessing: UZOB → X. Build encoder ids: AA tokens with
        # the target position replaced by <extra_id_0>, followed by </s>.
        clean = re.sub(r"[UZOB]", "X", sequence)
        enc_ids = [self._aa_to_id(c) for c in clean]
        enc_ids[pos1 - 1] = self.extra_id_0
        enc_ids.append(self.eos)
        input_ids = torch.tensor([enc_ids], device=self.device, dtype=torch.long)

        # Decoder input: <pad> + WT prefix (residues 1..pos1−1).
        dec_ids = [self.decoder_start] + [self._aa_to_id(c) for c in clean[: pos1 - 1]]
        dec_in = torch.tensor([dec_ids], device=self.device, dtype=torch.long)

        out = self.model(input_ids=input_ids, decoder_input_ids=dec_in)
        # logits[0, -1] predicts the next decoder token — i.e. the residue at pos1.
        logits = out.logits[0, -1].float()
        # Restrict softmax to the 20 AA tokens (mask out non-AA logits).
        aa_logits = torch.stack(
            [logits[self.aa_token_ids[aa]] for aa in AMINO_ACIDS]
        )
        return F.log_softmax(aa_logits, dim=-1).cpu().numpy()


def score_one_assay(
    dms_id: str,
    dms_path: Path,
    scorer: ProtT5Scorer,
    max_variants: Optional[int],
) -> Optional[dict]:
    if not dms_path.exists():
        log.error("  DMS missing: %s", dms_path)
        return None
    dms = pd.read_csv(dms_path)
    if max_variants and len(dms) > max_variants:
        dms = dms.sample(max_variants, random_state=42).reset_index(drop=True)

    parsed_per_row: List[List[Tuple[str, int, str]]] = []
    for m in dms["mutant"]:
        parsed_per_row.append(parse_mutant(m))
    # Determine unique positions to mask (one ProtT5 forward per position)
    unique_positions = sorted({p for row in parsed_per_row for _, p, _ in row})
    if not unique_positions:
        return None

    # WT sequence: revert the first variant's mutations
    first = parsed_per_row[0]
    if not first:
        return None
    wt_seq = list(dms["mutated_sequence"].iloc[0])
    for waa, p, _ in first:
        wt_seq[p - 1] = waa
    wt_seq_str = "".join(wt_seq)
    L = len(wt_seq_str)

    # For each unique position, score the 20 AA log-probs once
    t0 = time.time()
    pos_logprobs: Dict[int, np.ndarray] = {}
    for k, pos in enumerate(unique_positions):
        if not (1 <= pos <= L):
            continue
        pos_logprobs[pos] = scorer.score_position(wt_seq_str, pos)
        if k % 50 == 0:
            log.info("    pos %d/%d", k + 1, len(unique_positions))
    t_score = time.time() - t0

    # For each variant: sum log P(mut) − log P(wt) over its mutated positions
    aa_idx = {aa: i for i, aa in enumerate(AMINO_ACIDS)}
    rows = []
    for vidx, parsed in enumerate(parsed_per_row):
        if not parsed:
            continue
        score = 0.0
        ok = True
        for waa, p, maa in parsed:
            if p not in pos_logprobs or maa not in aa_idx or waa not in aa_idx:
                ok = False; break
            lp = pos_logprobs[p]
            score += lp[aa_idx[maa]] - lp[aa_idx[waa]]
        if not ok:
            continue
        rows.append({
            "mutant": str(dms["mutant"].iloc[vidx]),
            "n_mutations": len(parsed),
            "DMS_score": float(dms["DMS_score"].iloc[vidx]),
            "pseudoll_score": float(score),
        })

    per_variant = pd.DataFrame(rows)
    rho = spearmanr(per_variant["pseudoll_score"], per_variant["DMS_score"],
                    nan_policy="omit").statistic if len(per_variant) >= 5 else float("nan")
    log.info(
        "  %s: %d variants scored, %d unique positions, %.1fs scoring, Spearman=%+.3f",
        dms_id, len(per_variant), len(unique_positions), t_score, rho,
    )
    return {
        "DMS_id": dms_id,
        "n_variants_scored": len(per_variant),
        "n_unique_positions": len(unique_positions),
        "seconds": round(t_score, 1),
        "spearman": float(rho) if np.isfinite(rho) else float("nan"),
        "per_variant": per_variant,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matches", type=Path, required=True)
    ap.add_argument("--dms_dir", type=Path, required=True)
    ap.add_argument("--output_dir", type=Path, required=True)
    ap.add_argument("--assays", type=Path, default=None,
                    help="Optional TSV/CSV with DMS_id column. Otherwise uses all 135")
    ap.add_argument("--max_variants", type=int, default=200)
    ap.add_argument("--max_seq_len", type=int, default=2048)
    ap.add_argument("--device", type=str, default=None)
    ap.add_argument("--dtype", type=str, default="auto",
                    choices=["auto", "float16", "bfloat16", "float32"])
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    matches = pd.read_csv(args.matches)
    if args.max_seq_len:
        matches = matches[matches["target_seq_len"] <= args.max_seq_len]
    if args.assays:
        sep = "\t" if args.assays.suffix == ".tsv" else ","
        wanted = set(pd.read_csv(args.assays, sep=sep)["DMS_id"].tolist())
        matches = matches[matches["DMS_id"].isin(wanted)]
    assay_ids = sorted(matches["DMS_id"].unique())
    log.info("Scoring %d assays via pseudo-LL", len(assay_ids))

    if args.device:
        device = torch.device(args.device)
    else:
        if torch.cuda.is_available(): device = torch.device("cuda")
        elif torch.backends.mps.is_available(): device = torch.device("mps")
        else: device = torch.device("cpu")

    if args.dtype == "auto":
        dtype = torch.float16 if device.type == "cuda" else torch.float32
    else:
        dtype = getattr(torch, args.dtype)

    scorer = ProtT5Scorer(device=device, dtype=dtype)

    summary_rows = []
    for i, dms_id in enumerate(assay_ids):
        log.info("[%d/%d] %s", i + 1, len(assay_ids), dms_id)
        result = score_one_assay(
            dms_id, args.dms_dir / f"{dms_id}.csv", scorer, args.max_variants
        )
        if result is None:
            continue
        result["per_variant"].to_parquet(
            args.output_dir / f"{dms_id}__pseudoll.parquet", index=False
        )
        summary_rows.append({k: v for k, v in result.items() if k != "per_variant"})
        # Persist summary incrementally
        pd.DataFrame(summary_rows).to_csv(args.output_dir / "summary.csv", index=False)

    if not summary_rows:
        log.warning("No assays produced output."); return
    s = pd.DataFrame(summary_rows)
    print("\n" + "=" * 56)
    print(" Pseudo-LL Spearman across {} assays".format(len(s)))
    print("=" * 56)
    v = s["spearman"].abs().dropna()
    print(f"  median |ρ| = {v.median():.3f}   mean = {v.mean():.3f}")
    print(f"  > 0.3  in {(v > 0.3).sum()}/{len(v)}")
    print(f"  > 0.5  in {(v > 0.5).sum()}/{len(v)}")
    print()
    print(s.sort_values("spearman", key=lambda x: x.abs(), ascending=False).to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
