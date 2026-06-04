"""
Position-residualized Spearman on the pooled-metrics parquets.

For each (assay, concept, feature) we currently report a *global* Spearman
between a per-variant metric (e.g. pool_mean_abs) and DMS_score. That number
mixes two signals:
  - position importance: at conserved sites *any* substitution is bad,
    so variants at those sites have both extreme metric values and extreme
    DMS scores → spurious correlation
  - substitution-specificity: for variants at the SAME site, does the metric
    track which substitution is more harmful?

This script partials out the per-position mean from both the metric and
DMS_score, then re-correlates the residuals. The within-position component
that survives is the pure substitution-specific signal — analogous to the
position-grouped CV we ran on the raw-residual ridge probe.

For each (assay, concept, feature):
    1. Filter to single-mutant rows (n_mutations == 1) so position is well-defined.
    2. For each variant: residual = value − mean(value | same position).
    3. Spearman of the residuals.

Reports both the raw Spearman (matches existing summary.csv) and the
position-residualized Spearman, so the gap is the leakage estimate.

Usage:
    uv run --project repos/crosscode --with scipy --with pyarrow python \\
      repos/sparse-crosscoders-prott5/position_residualized_spearman.py \\
      --pooled_dir data/proteingym/pooled_metrics \\
      --matches data/proteingym_concept_matches.csv \\
      --pairings data/crosscoder_eval/.../heldout_all_top_pairings.csv \\
      --output data/proteingym/position_residualized_summary.csv
"""
from __future__ import annotations
import argparse, logging, warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import ConstantInputWarning, spearmanr

warnings.filterwarnings("ignore", category=ConstantInputWarning)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("pos_resid")


def load_candidate_features(pairings_path: Path, concept: str) -> List[int]:
    df = pd.read_csv(pairings_path)
    sub = df[df["concept"] == concept].copy()
    if sub.empty:
        return []
    sub = sub.query("tp_per_domain >= 2 or tp >= 2")
    sub = sub.sort_values(
        ["f1_per_domain", "recall_per_domain", "tp"], ascending=False
    ).drop_duplicates("feature", keep="first")
    return sub["feature"].astype(int).tolist()


def residualize(df_sub: pd.DataFrame, value_col: str, pos_col: str = "position") -> np.ndarray:
    """Subtract per-position mean from a per-row column. Returns residuals."""
    pos_mean = df_sub.groupby(pos_col)[value_col].transform("mean")
    return (df_sub[value_col] - pos_mean).values


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pooled_dir", type=Path, required=True)
    ap.add_argument("--matches", type=Path, required=True)
    ap.add_argument("--pairings", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--min_positions", type=int, default=3,
                    help="Skip (assay,concept,feature) where <N unique positions with ≥2 variants")
    ap.add_argument("--min_variants", type=int, default=20,
                    help="Skip (assay,concept,feature) where total single-mutants < N")
    args = ap.parse_args()

    matches = pd.read_csv(args.matches)
    pair = pd.read_csv(args.pairings)
    f1map = pair.set_index(["concept", "feature"])["f1_per_domain"].to_dict()

    # concept → list of paired features
    concept_feats: Dict[str, List[int]] = {}
    for c in matches["matched_concept"].unique():
        concept_feats[c] = load_candidate_features(args.pairings, c)
    # (dms_id, feature) → list of concepts that feature is paired to for this assay
    concept_set: Dict[Tuple[str, int], List[str]] = defaultdict(list)
    for _, r in matches.iterrows():
        for f in concept_feats[r["matched_concept"]]:
            concept_set[(r["DMS_id"], int(f))].append(r["matched_concept"])

    metric_cols = [
        "delta_at_pos", "at_pos_mut", "pool_mean_abs", "pool_max_abs",
        "pool_mean_mut", "pool_max_mut",
    ]

    parquets = sorted(args.pooled_dir.glob("*__pooled.parquet"))
    log.info("Processing %d parquets", len(parquets))

    rows = []
    for i, p in enumerate(parquets):
        dms_id = p.stem.replace("__pooled", "")
        if i % 20 == 0:
            log.info("  [%d/%d] %s", i + 1, len(parquets), dms_id)
        try:
            df = pd.read_parquet(p)
        except Exception as e:
            log.warning("  read failed for %s: %s", dms_id, e); continue
        # Single-mutants only — multi-mutant residualization is ambiguous
        df = df[df["n_mutations"] == 1].copy()
        if len(df) == 0:
            continue

        for feat, g in df.groupby("feature"):
            if len(g) < args.min_variants:
                continue
            # Require ≥min_positions unique positions, each with ≥2 variants
            pos_counts = g["position"].value_counts()
            multi_pos_count = (pos_counts >= 2).sum()
            if multi_pos_count < args.min_positions:
                continue
            # Restrict to positions with ≥2 variants for residualization to be meaningful
            valid_positions = pos_counts.index[pos_counts >= 2]
            g_valid = g[g["position"].isin(valid_positions)].copy()
            if len(g_valid) < args.min_variants:
                continue
            y = g_valid["DMS_score"].values
            y_resid = residualize(g_valid, "DMS_score")
            concepts = concept_set.get((dms_id, int(feat)), [None])
            base_row = {
                "DMS_id": dms_id,
                "feature": int(feat),
                "n_single_mut": int(len(g_valid)),
                "n_unique_pos": int(g_valid["position"].nunique()),
                "n_pos_with_multi": int(multi_pos_count),
            }
            for concept in concepts:
                row = dict(base_row)
                row["concept"] = concept
                row["f1_per_domain"] = (f1map.get((concept, int(feat)), float("nan"))
                                          if concept else float("nan"))
                for col in metric_cols:
                    if col not in g_valid.columns:
                        row[f"{col}_raw"] = float("nan")
                        row[f"{col}_resid"] = float("nan")
                        continue
                    x = g_valid[col].values
                    if np.std(x) < 1e-10 or np.std(y) < 1e-10:
                        row[f"{col}_raw"] = float("nan")
                        row[f"{col}_resid"] = float("nan")
                        continue
                    # Raw Spearman
                    rho_raw = spearmanr(x, y, nan_policy="omit").statistic
                    # Residualized Spearman
                    x_resid = residualize(g_valid, col)
                    if np.std(x_resid) < 1e-10:
                        rho_resid = float("nan")
                    else:
                        rho_resid = spearmanr(x_resid, y_resid, nan_policy="omit").statistic
                    row[f"{col}_raw"] = float(rho_raw) if np.isfinite(rho_raw) else float("nan")
                    row[f"{col}_resid"] = float(rho_resid) if np.isfinite(rho_resid) else float("nan")
                rows.append(row)

    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    log.info("Wrote %d rows -> %s", len(out), args.output)

    # ---- Headline summary ----
    print("\n=== HEADLINE: raw vs position-residualized |Spearman| medians ===")
    print(f"{'metric':<22s}  {'n':>5s}  {'raw_med':>8s}  {'resid_med':>10s}  {'gap':>7s}  "
          f"{'raw>0.3':>8s}  {'resid>0.3':>10s}  {'resid>0.5':>10s}")
    for col in metric_cols:
        v_raw = out[f"{col}_raw"].abs().dropna()
        v_resid = out[f"{col}_resid"].abs().dropna()
        if len(v_resid) == 0: continue
        gap = v_raw.median() - v_resid.median()
        print(f"  {col:<22s}  {len(v_resid):>5d}  {v_raw.median():>8.3f}  "
              f"{v_resid.median():>10.3f}  {gap:>+7.3f}  "
              f"{(v_raw>0.3).sum():>8d}  {(v_resid>0.3).sum():>10d}  {(v_resid>0.5).sum():>10d}")

    print("\n=== Best feature per (assay,concept) ===")
    for col in ["pool_mean_abs"]:
        rawcol, residcol = f"{col}_raw", f"{col}_resid"
        best_raw = out.assign(absx=lambda d: d[rawcol].abs()).sort_values(
            "absx", ascending=False).drop_duplicates(["DMS_id","concept"])
        best_resid = out.assign(absx=lambda d: d[residcol].abs()).sort_values(
            "absx", ascending=False).drop_duplicates(["DMS_id","concept"])
        v_raw = best_raw[rawcol].abs().dropna()
        v_resid = best_resid[residcol].abs().dropna()
        print(f"  {col}: raw best/pair median={v_raw.median():.3f} (n={len(v_raw)})  |  "
              f"resid best/pair median={v_resid.median():.3f} (n={len(v_resid)})  "
              f">0.3:{(v_resid>0.3).sum()}/{len(v_resid)}  >0.5:{(v_resid>0.5).sum()}")


if __name__ == "__main__":
    import sys; sys.exit(main())
