"""
Harvest full-sequence pooled SAE metrics for ProteinGym matched assays.

Extends the per-residue-only score_proteingym_features.py to also save
*sequence-pooled* statistics, which the prior analysis suggested may carry
mutation-specific signal where the per-residue Δact does not.

For each (variant, candidate feature) we compute the feature activation
(post-BatchTopK — the only value the crosscoder emits):

  at_pos          activation at the mutated residue (or centroid for
                  multi-mutants)
  delta_at_pos    act_mut[pos] − act_wt[pos]
  pool_mean_abs   mean over the *whole sequence* of |act_mut[i] − act_wt[i]|
  pool_max_abs    max  over the whole sequence of |act_mut[i] − act_wt[i]|
  pool_mean_mut   mean over the whole sequence of act_mut[i]
  pool_max_mut    max  over the whole sequence of act_mut[i]

Output (per assay):
    <output_dir>/<DMS_id>__pooled.parquet
      long-format, one row per (variant, feature) across all matched concepts
      for the assay. Includes DMS_score, position, n_mutations, feature, the
      6 metrics, plus fire_rate aggregated across the whole sequence (fraction
      of positions where the feature is non-zero).

After harvest is complete, run:
    --compute_summary   pass to also write per-(assay,concept,feature)
                        Spearman summaries to summary.csv
or invoke the helper at the bottom of this file directly.

The harvest pipeline mirrors score_proteingym_features.py and reuses the
same per-assay caching trick (one ProtT5 + crosscoder pass per assay across
all matched concepts).

Usage (smoke test):
    uv run --project repos/crosscode --with scipy --with pyarrow python \\
      repos/sparse-crosscoders-prott5/harvest_pooled_metrics.py \\
      --crosscoder_dir model_checkpoints/.../crashed_epoch_0_step_2519836 \\
      --pairings data/crosscoder_eval/pre-auxfix/real/.../heldout_all_top_pairings.csv \\
      --matches data/proteingym_concept_matches.csv \\
      --dms_dir data/external/DMS_ProteinGym_substitutions \\
      --output_dir data/proteingym/pooled_metrics \\
      --max_variants 30 \\
      --assays /tmp/two_assay.tsv

Usage (cluster full run): omit --max_variants and --assays.
"""

from __future__ import annotations

import argparse
import gc
import logging
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import torch
from scipy.stats import ConstantInputWarning, spearmanr

# Pinned schema so the per-batch streaming writer (see harvest_one_assay) gets
# a consistent type layout — without this, pyarrow infers per batch and a column
# that happens to be all-zero in batch 1 vs varied in batch 2 can cause a type
# mismatch crash. Floats are stored as float32 to halve on-disk size.
POOLED_SCHEMA = pa.schema([
    ("DMS_id", pa.string()),
    ("mutant", pa.string()),
    ("position", pa.int32()),
    ("centroid_pos", pa.int32()),
    ("n_mutations", pa.int16()),
    ("DMS_score", pa.float32()),
    ("feature", pa.int32()),
    ("at_pos_wt", pa.float32()),
    ("at_pos_mut", pa.float32()),
    ("delta_at_pos", pa.float32()),
    ("pool_mean_abs", pa.float32()),
    ("pool_max_abs", pa.float32()),
    ("pool_mean_mut", pa.float32()),
    ("pool_max_mut", pa.float32()),
    ("fire_rate_var", pa.float32()),
])

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "repos" / "InterPLM"))
sys.path.insert(0, str(_REPO_ROOT / "repos" / "crosscode"))

from interplm.embedders.prott5 import ProtT5CrosscoderEmbedder  # noqa: E402
from interplm.sae.inference import load_sae  # noqa: E402

warnings.filterwarnings("ignore", category=ConstantInputWarning)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("harvest_pooled")

MUTANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z\*])$")


def parse_mutant_positions(mutant_field: str) -> List[Tuple[str, int, str]]:
    out: List[Tuple[str, int, str]] = []
    for token in str(mutant_field).split(":"):
        m = MUTANT_RE.match(token.strip())
        if not m:
            return []
        out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


def load_candidate_features(
    pairings_path: Path, concept: str
) -> List[Tuple[int, float]]:
    df = pd.read_csv(pairings_path)
    sub = df[df["concept"] == concept].copy()
    if sub.empty:
        return []
    sub = sub.query("tp_per_domain >= 2 or tp >= 2")
    sub = sub.sort_values(
        ["f1_per_domain", "recall_per_domain", "tp"], ascending=False
    ).drop_duplicates("feature", keep="first")
    return list(
        zip(
            sub["feature"].astype(int).tolist(),
            sub["f1_per_domain"].astype(float).tolist(),
        )
    )


@torch.no_grad()
def encode_full_features(
    crosscoder, x: torch.Tensor, feat_idx: torch.Tensor, device: torch.device
) -> np.ndarray:
    """Run the encoder and return per-residue feature activations (post-BatchTopK)
    restricted to the feature subset.

    x: [L_seq, 1, 24, 1024] (the ProtT5 embedder's "batch" dim is amino-acid tokens)
    feat_idx: 1D tensor of feature indices

    Returns:
        act: (L_seq, n_feats) float32 — after activation_fn (BatchTopK at inference)
    """
    inner = crosscoder.crosscoder
    pre = inner.get_pre_bias_BL(x.to(device=device, dtype=next(inner.parameters()).dtype))
    if inner.b_enc_L is not None:
        pre = pre + inner.b_enc_L
    return inner.activation_fn.forward(pre)[:, feat_idx].float().cpu().numpy()


def harvest_one_assay(
    dms_id: str,
    concepts: List[str],
    dms_dir: Path,
    pairings_path: Path,
    embedder: ProtT5CrosscoderEmbedder,
    crosscoder,
    device: torch.device,
    output_dir: Path,
    batch_size: int,
    max_variants: Optional[int],
) -> int:
    """Returns number of (variant, feature) rows written, or 0 if skipped."""
    # Build candidate-feature union across all matched concepts for this assay
    feats_per_concept: Dict[str, List[int]] = {}
    f1_per_feat: Dict[int, float] = {}
    for c in concepts:
        cands = load_candidate_features(pairings_path, c)
        if not cands:
            continue
        feats_per_concept[c] = [f for f, _ in cands]
        for f, s in cands:
            f1_per_feat[f] = max(f1_per_feat.get(f, 0.0), s)
    if not feats_per_concept:
        log.warning("  no candidate features for %s; skipping", dms_id)
        return 0
    union_feats = sorted({f for fs in feats_per_concept.values() for f in fs})
    feat_to_idx = {f: i for i, f in enumerate(union_feats)}
    feat_idx = torch.tensor(union_feats, dtype=torch.long)

    dms_path = dms_dir / f"{dms_id}.csv"
    if not dms_path.exists():
        log.error("  DMS file missing: %s", dms_path)
        return 0
    dms = pd.read_csv(dms_path)
    if max_variants and len(dms) > max_variants:
        dms = dms.sample(max_variants, random_state=42).reset_index(drop=True)
    if len(dms) == 0:
        return 0

    # WT reconstruction: revert first variant's mutated_sequence
    first_mut = parse_mutant_positions(str(dms["mutant"].iloc[0]))
    if not first_mut:
        log.error("  could not parse mutant '%s'", dms["mutant"].iloc[0])
        return 0
    wt_seq = list(dms["mutated_sequence"].iloc[0])
    for waa, p, _ in first_mut:
        wt_seq[p - 1] = waa
    wt_seq_str = "".join(wt_seq)
    log.info(
        "  %s: %d variants, %d concepts, %d union features, L_wt=%d",
        dms_id, len(dms), len(feats_per_concept), len(union_feats), len(wt_seq_str),
    )

    # WT activations (full per-residue, for the union feature subset)
    wt_emb = embedder.extract_embeddings([wt_seq_str], batch_size=1)  # (L_wt, 1, 24, 1024)
    wt_act = encode_full_features(crosscoder, wt_emb, feat_idx, device)
    # (L_wt, n_union)
    L_wt = wt_act.shape[0]

    # Streaming write: open a ParquetWriter on the first non-empty batch and
    # append one row-group per batch. Peak memory is bounded by the largest
    # single batch (~batch_size × n_union × ~22 fields ≈ kilobytes), which lets
    # us safely process assays like PHOT_CHLRE (167 k variants) without OOM.
    out_path = output_dir / f"{dms_id}__pooled.parquet"
    tmp_path = out_path.with_suffix(".parquet.tmp")
    writer: Optional[pq.ParquetWriter] = None
    total_rows = 0
    n_batches = (len(dms) + batch_size - 1) // batch_size
    try:
        for bi in range(n_batches):
            s, e = bi * batch_size, min((bi + 1) * batch_size, len(dms))
            batch_df = dms.iloc[s:e]
            seqs = batch_df["mutated_sequence"].tolist()
            bundle = embedder.extract_embeddings_with_boundaries(
                seqs, layer=-1, batch_size=batch_size
            )
            emb_all = bundle["embeddings"]
            boundaries = bundle["boundaries"]
            all_act = encode_full_features(
                crosscoder, emb_all, feat_idx, device
            )

            batch_rows: List[dict] = []
            for (_, vrow), (b_start, b_end) in zip(batch_df.iterrows(), boundaries):
                parsed = parse_mutant_positions(str(vrow["mutant"]))
                if not parsed:
                    continue
                mut_act = all_act[b_start:b_end]  # (L_var, n_union)
                L_var = mut_act.shape[0]
                L_cmp = min(L_var, L_wt)
                d = mut_act[:L_cmp] - wt_act[:L_cmp]
                m = mut_act[:L_cmp]
                w = wt_act[:L_cmp]
                ad = np.abs(d)
                pool_mean_abs = ad.mean(axis=0)
                pool_max_abs  = ad.max(axis=0)
                pool_mean_mut = m.mean(axis=0)
                pool_max_mut  = m.max(axis=0)
                fire_var = (m != 0).astype(np.float32).mean(axis=0)
                n_mut = len(parsed)
                positions = [p for _, p, _ in parsed if 0 <= p - 1 < L_cmp]
                if not positions:
                    continue
                centroid_pos = int(round(np.mean(positions)))
                mutant_str = str(vrow["mutant"])
                dms_score = float(vrow["DMS_score"])
                for pos in positions:
                    wi = w[pos - 1]; mi = m[pos - 1]; di = d[pos - 1]
                    for feat in union_feats:
                        i = feat_to_idx[feat]
                        batch_rows.append({
                            "DMS_id": dms_id,
                            "mutant": mutant_str,
                            "position": pos,
                            "centroid_pos": centroid_pos,
                            "n_mutations": n_mut,
                            "DMS_score": dms_score,
                            "feature": int(feat),
                            "at_pos_wt": float(wi[i]),
                            "at_pos_mut": float(mi[i]),
                            "delta_at_pos": float(di[i]),
                            "pool_mean_abs": float(pool_mean_abs[i]),
                            "pool_max_abs":  float(pool_max_abs[i]),
                            "pool_mean_mut": float(pool_mean_mut[i]),
                            "pool_max_mut":  float(pool_max_mut[i]),
                            "fire_rate_var": float(fire_var[i]),
                        })

            if batch_rows:
                table = pa.Table.from_pylist(batch_rows, schema=POOLED_SCHEMA)
                if writer is None:
                    writer = pq.ParquetWriter(tmp_path, POOLED_SCHEMA, compression="snappy")
                writer.write_table(table)
                total_rows += len(batch_rows)
                # Explicitly drop to keep peak bounded.
                del batch_rows, table

            if bi % 50 == 0:
                log.info("    batch %d/%d (cumulative rows: %d)", bi + 1, n_batches, total_rows)
    finally:
        if writer is not None:
            writer.close()

    if total_rows == 0:
        if tmp_path.exists():
            tmp_path.unlink()
        return 0
    # Atomic rename so partial writes never get committed under the final name.
    tmp_path.replace(out_path)
    log.info("  wrote %d (variant, feature) rows -> %s", total_rows, out_path.name)
    return total_rows


def compute_summary(output_dir: Path, matches_path: Path, pairings_path: Path) -> None:
    """Per (assay, concept, feature) Spearman of each candidate metric vs DMS_score."""
    matches = pd.read_csv(matches_path)
    pair = pd.read_csv(pairings_path)
    f1map = pair.set_index(["concept", "feature"])["f1_per_domain"].to_dict()

    summary_rows = []
    metric_cols = [
        "at_pos_mut", "delta_at_pos", "pool_mean_abs", "pool_max_abs",
        "pool_mean_mut", "pool_max_mut",
    ]
    # Build per-(assay, feature) → list of concepts mapping
    concept_set: Dict[Tuple[str, int], List[str]] = defaultdict(list)
    for _, r in matches.iterrows():
        # which features are paired to this concept?
        for f, _s in load_candidate_features(pairings_path, r["matched_concept"]):
            concept_set[(r["DMS_id"], int(f))].append(r["matched_concept"])

    for p in sorted(output_dir.glob("*__pooled.parquet")):
        dms_id = p.stem.replace("__pooled", "")
        df = pd.read_parquet(p)
        for feat, g in df.groupby("feature"):
            # For each concept this feature is paired to for this assay
            concepts = concept_set.get((dms_id, int(feat)), [None])
            # Score Spearman per metric column
            y = g["DMS_score"].values
            if len(y) < 5 or np.std(y) < 1e-10:
                continue
            for concept in concepts:
                row = {
                    "DMS_id": dms_id,
                    "concept": concept,
                    "feature": int(feat),
                    "n_variants": int(g["mutant"].nunique()),
                    "n_rows": int(len(g)),
                    "fire_rate_mean": float(g["fire_rate_var"].mean()),
                    "f1_per_domain": f1map.get((concept, int(feat)), float("nan"))
                                       if concept else float("nan"),
                }
                for col in metric_cols:
                    x = g[col].values
                    if np.std(x) < 1e-10:
                        row[col] = float("nan")
                        continue
                    rho = spearmanr(x, y, nan_policy="omit").statistic
                    row[col] = float(rho) if np.isfinite(rho) else float("nan")
                summary_rows.append(row)

    df = pd.DataFrame(summary_rows)
    out = output_dir / "summary.csv"
    df.to_csv(out, index=False)
    log.info("Wrote per-(assay,concept,feature) Spearman summary: %d rows -> %s",
             len(df), out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crosscoder_dir", type=Path, required=True)
    ap.add_argument("--checkpoint", type=str, default="ae_normalized.pt")
    ap.add_argument("--pairings", type=Path, required=True)
    ap.add_argument("--matches", type=Path, required=True)
    ap.add_argument("--dms_dir", type=Path, required=True)
    ap.add_argument("--output_dir", type=Path, required=True)
    ap.add_argument("--batch_size", type=int, default=4)
    ap.add_argument("--max_variants", type=int, default=None)
    ap.add_argument("--max_seq_len", type=int, default=2048)
    ap.add_argument("--assays", type=Path, default=None,
                    help="Optional TSV/CSV with DMS_id column to restrict the run")
    ap.add_argument("--device", type=str, default=None)
    ap.add_argument("--skip_existing", action="store_true")
    ap.add_argument("--compute_summary", action="store_true",
                    help="Also compute per-feature Spearman summary at the end")
    ap.add_argument("--summary_only", action="store_true",
                    help="Skip harvesting; just (re)compute summary from existing parquets")
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.summary_only:
        compute_summary(args.output_dir, args.matches, args.pairings)
        return

    matches = pd.read_csv(args.matches)
    if args.max_seq_len:
        before = matches["DMS_id"].nunique()
        matches = matches[matches["target_seq_len"] <= args.max_seq_len].copy()
        after = matches["DMS_id"].nunique()
        if after < before:
            log.info("Length filter dropped %d assays (target_seq_len > %d)",
                     before - after, args.max_seq_len)
    if args.assays:
        sep = "\t" if args.assays.suffix == ".tsv" else ","
        wanted = set(pd.read_csv(args.assays, sep=sep)["DMS_id"].tolist())
        matches = matches[matches["DMS_id"].isin(wanted)].copy()

    grouped: Dict[str, List[str]] = defaultdict(list)
    for _, r in matches.iterrows():
        grouped[r["DMS_id"]].append(r["matched_concept"])
    log.info(
        "Will harvest %d assays / %d (assay,concept) pairs",
        len(grouped), sum(len(v) for v in grouped.values()),
    )

    if args.device:
        device = torch.device(args.device)
    else:
        if torch.cuda.is_available():
            device = torch.device("cuda")
        elif torch.backends.mps.is_available():
            device = torch.device("mps")
        else:
            device = torch.device("cpu")
    log.info("Device: %s", device)

    log.info("Loading ProtT5...")
    embedder = ProtT5CrosscoderEmbedder(device=str(device))
    log.info("Loading crosscoder: %s / %s", args.crosscoder_dir, args.checkpoint)
    crosscoder = load_sae(args.crosscoder_dir, device=str(device), model_name=args.checkpoint)
    crosscoder.eval()

    total_rows = 0
    for i, (dms_id, concepts) in enumerate(grouped.items()):
        out_path = args.output_dir / f"{dms_id}__pooled.parquet"
        if args.skip_existing and out_path.exists():
            log.info("[%d/%d] %s — skip (exists)", i + 1, len(grouped), dms_id)
            continue
        log.info("[%d/%d] %s", i + 1, len(grouped), dms_id)
        try:
            n = harvest_one_assay(
                dms_id, concepts, args.dms_dir, args.pairings,
                embedder, crosscoder, device,
                args.output_dir, args.batch_size, args.max_variants,
            )
            total_rows += n
        except Exception as e:
            log.exception("  failed on %s: %s", dms_id, e)
        # Defensive cleanup between assays. The cluster's 92 GB OOM on the
        # previous run was almost certainly cumulative growth across 95 assays
        # (rows list could only explain ~8 GB on its own). Cycling the GC and
        # releasing PyTorch's CUDA caching allocator buffers back to the OS
        # bounds per-assay carry-over to near zero.
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        elif torch.backends.mps.is_available():
            torch.mps.empty_cache()
    log.info("Total rows written: %d", total_rows)

    if args.compute_summary:
        compute_summary(args.output_dir, args.matches, args.pairings)


if __name__ == "__main__":
    sys.exit(main())
