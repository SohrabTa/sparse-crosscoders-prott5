"""
Per-(assay, feature) Spearman scan over ALL 8192 crosscoder features.

The pooled-metrics run answered "do the InterPLM-paired candidate features
carry fitness signal when read with pool_mean_abs?" → yes. This script answers
the follow-up: **is there more signal in the other 8000 features that InterPLM
didn't pair?** If so, we have features that capture fitness biology beyond the
named-concept space. (Result: yes — unpaired live features reach median best
|ρ| 0.477 per assay vs 0.284 for the paired candidates.)

For each assay (single ProtT5 + crosscoder forward per variant):
  1. Encode WT, then each variant batch.
  2. For each variant, compute pool_mean_abs over all 8192 features (one
     pool-mean scalar per feature per variant).
  3. Per feature: Spearman(per_variant_metric, DMS_score).
  4. Save (assay, feature, metric, spearman) for all 135 × 8192 ≈ 1.1 M rows.

Storage: one float per (assay, feature, metric) → tiny output (~20 MB).
Memory peak per assay: subsampled to ≤max_variants × 8192 × 4 bytes ×
n_metrics → keep ≤ ~13 GB even on huge assays via cap.

Output:
    <output_dir>/full_feature_spearman.csv
      DMS_id, feature, n_variants, fire_rate_mean,
      pool_mean_abs_sp, delta_at_pos_sp

After harvest, join against `heldout_all_top_pairings.csv` to label features
that ARE / ARE NOT InterPLM-paired and look at the unpaired top-rankers.
"""
from __future__ import annotations
import argparse, gc, logging, re, sys, warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
from scipy.stats import ConstantInputWarning, spearmanr

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "repos" / "InterPLM"))
sys.path.insert(0, str(_REPO_ROOT / "repos" / "crosscode"))

from interplm.embedders.prott5 import ProtT5CrosscoderEmbedder  # noqa: E402
from interplm.sae.inference import load_sae  # noqa: E402

warnings.filterwarnings("ignore", category=ConstantInputWarning)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("harvest_full")

MUTANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z\*])$")


def parse_mutant_positions(mutant_field: str) -> List[Tuple[str, int, str]]:
    out: List[Tuple[str, int, str]] = []
    for token in str(mutant_field).split(":"):
        m = MUTANT_RE.match(token.strip())
        if not m:
            return []
        out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


@torch.no_grad()
def encode_full_8192(
    crosscoder, x: torch.Tensor, device: torch.device
) -> torch.Tensor:
    """Return per-residue feature activations (post-BatchTopK) for ALL 8192 features.
    x: [L_seq, 1, 24, 1024]. Returns: (L_seq, 8192) — on GPU."""
    inner = crosscoder.crosscoder
    pre = inner.get_pre_bias_BL(x.to(device=device, dtype=next(inner.parameters()).dtype))
    if inner.b_enc_L is not None:
        pre = pre + inner.b_enc_L
    return inner.activation_fn.forward(pre).float()


def harvest_one_assay(
    dms_id: str,
    dms_dir: Path,
    embedder: ProtT5CrosscoderEmbedder,
    crosscoder,
    n_latents: int,
    device: torch.device,
    batch_size: int,
    max_variants: int,
) -> Optional[pd.DataFrame]:
    dms_path = dms_dir / f"{dms_id}.csv"
    if not dms_path.exists():
        log.error("  DMS missing: %s", dms_path); return None
    dms = pd.read_csv(dms_path)
    if max_variants and len(dms) > max_variants:
        dms = dms.sample(max_variants, random_state=42).reset_index(drop=True)
    if len(dms) == 0:
        return None

    first_mut = parse_mutant_positions(str(dms["mutant"].iloc[0]))
    if not first_mut: return None
    wt_seq = list(dms["mutated_sequence"].iloc[0])
    for waa, p, _ in first_mut:
        wt_seq[p - 1] = waa
    wt_seq_str = "".join(wt_seq)
    L_wt = len(wt_seq_str)
    log.info("  %s: %d variants  L_wt=%d", dms_id, len(dms), L_wt)

    # WT activations once (kept on GPU)
    wt_emb = embedder.extract_embeddings([wt_seq_str], batch_size=1)
    wt_act = encode_full_8192(crosscoder, wt_emb, device)
    # (L_wt, 8192)

    # Per-variant accumulators — float32 on CPU
    n = len(dms)
    pool_mean_abs = np.zeros((n, n_latents), dtype=np.float32)
    delta_at_pos  = np.zeros((n, n_latents), dtype=np.float32)
    fire_per_var  = np.zeros((n, n_latents), dtype=np.float32)
    valid = np.zeros(n, dtype=bool)
    pos_used = np.zeros(n, dtype=np.int32)
    dms_scores = np.zeros(n, dtype=np.float32)

    seqs = dms["mutated_sequence"].tolist()
    n_batches = (n + batch_size - 1) // batch_size
    for bi in range(n_batches):
        s, e = bi * batch_size, min((bi + 1) * batch_size, n)
        batch_seqs = seqs[s:e]
        bundle = embedder.extract_embeddings_with_boundaries(
            batch_seqs, layer=-1, batch_size=batch_size
        )
        emb_all = bundle["embeddings"]
        boundaries = bundle["boundaries"]
        all_act = encode_full_8192(crosscoder, emb_all, device)

        for j, (b_s, b_e) in enumerate(boundaries):
            idx = s + j
            mutant = str(dms["mutant"].iloc[idx])
            parsed = parse_mutant_positions(mutant)
            if not parsed:
                continue
            positions = [p for _, p, _ in parsed if 0 <= p - 1 < L_wt]
            if not positions:
                continue
            mut_a = all_act[b_s:b_e]
            L_cmp = min(mut_a.shape[0], L_wt)
            da = (mut_a[:L_cmp] - wt_act[:L_cmp])
            # Pool-over-sequence stats
            pool_mean_abs[idx] = da.abs().mean(dim=0).cpu().numpy()
            fire_per_var[idx]  = (mut_a[:L_cmp] != 0).float().mean(dim=0).cpu().numpy()
            # At-pos stats (averaged over mutated positions for multi-mutant)
            pos_idx = [p - 1 for p in positions if p - 1 < L_cmp]
            if pos_idx:
                delta_at_pos[idx] = da[pos_idx].mean(dim=0).cpu().numpy()
            valid[idx] = True
            pos_used[idx] = int(round(np.mean(positions)))
            dms_scores[idx] = float(dms["DMS_score"].iloc[idx])
        if bi % 50 == 0:
            log.info("    batch %d/%d", bi + 1, n_batches)
        del all_act, emb_all
        if device.type == "cuda":
            torch.cuda.empty_cache()

    # Subset to valid variants
    valid_idx = np.where(valid)[0]
    if len(valid_idx) < 5:
        return None
    y = dms_scores[valid_idx]
    if np.std(y) < 1e-10:
        return None

    # Per-feature Spearman vs DMS
    log.info("  computing Spearman across %d valid variants × %d features",
             len(valid_idx), n_latents)
    rows = []
    for fi in range(n_latents):
        row = {"DMS_id": dms_id, "feature": int(fi),
               "n_variants": int(len(valid_idx)),
               "fire_rate_mean": float(fire_per_var[valid_idx, fi].mean())}
        for col, arr in [
            ("pool_mean_abs_sp", pool_mean_abs),
            ("delta_at_pos_sp",  delta_at_pos),
        ]:
            x = arr[valid_idx, fi]
            if np.std(x) < 1e-10:
                row[col] = float("nan"); continue
            rho = spearmanr(x, y, nan_policy="omit").statistic
            row[col] = float(rho) if np.isfinite(rho) else float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crosscoder_dir", type=Path, required=True)
    ap.add_argument("--checkpoint", type=str, default="ae_normalized.pt")
    ap.add_argument("--matches", type=Path, required=True)
    ap.add_argument("--dms_dir", type=Path, required=True)
    ap.add_argument("--output_dir", type=Path, required=True)
    ap.add_argument("--batch_size", type=int, default=8)
    ap.add_argument("--max_variants", type=int, default=20000,
                    help="Subsample assays above this variant count (memory cap)")
    ap.add_argument("--max_seq_len", type=int, default=2048)
    ap.add_argument("--assays", type=Path, default=None)
    ap.add_argument("--device", type=str, default=None)
    ap.add_argument("--skip_existing", action="store_true")
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
    log.info("Will scan %d assays × 8192 features", len(assay_ids))

    if args.device:
        device = torch.device(args.device)
    else:
        if torch.cuda.is_available(): device = torch.device("cuda")
        elif torch.backends.mps.is_available(): device = torch.device("mps")
        else: device = torch.device("cpu")
    log.info("Device: %s", device)

    log.info("Loading ProtT5...")
    embedder = ProtT5CrosscoderEmbedder(device=str(device))
    log.info("Loading crosscoder...")
    crosscoder = load_sae(args.crosscoder_dir, device=str(device), model_name=args.checkpoint)
    crosscoder.eval()
    n_latents = crosscoder.dict_size
    log.info("  n_latents = %d", n_latents)

    for i, dms_id in enumerate(assay_ids):
        out_path = args.output_dir / f"{dms_id}__full_feat.parquet"
        if args.skip_existing and out_path.exists():
            log.info("[%d/%d] %s — skip", i + 1, len(assay_ids), dms_id)
            continue
        log.info("[%d/%d] %s", i + 1, len(assay_ids), dms_id)
        try:
            df = harvest_one_assay(
                dms_id, args.dms_dir, embedder, crosscoder, n_latents,
                device, args.batch_size, args.max_variants,
            )
        except Exception as e:
            log.exception("  failed: %s", e); continue
        if df is None or df.empty:
            log.warning("  no rows produced"); continue
        df.to_parquet(out_path, index=False)
        log.info("  wrote %d rows", len(df))
        gc.collect()
        if device.type == "cuda":
            torch.cuda.empty_cache()

    # Final aggregation: concatenate all per-assay parquets into one CSV
    log.info("Aggregating per-assay files...")
    parts = list(args.output_dir.glob("*__full_feat.parquet"))
    if not parts:
        log.warning("No parquets produced"); return
    out = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
    full_csv = args.output_dir / "full_feature_spearman.csv"
    out.to_csv(full_csv, index=False)
    log.info("Wrote final aggregate: %d rows -> %s", len(out), full_csv)


if __name__ == "__main__":
    sys.exit(main())
