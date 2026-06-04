"""Two literature-grounded ProteinGym analyses on the paired-feature pooled parquets.

A) Per-residue MotifAE-comparable readout (Hu et al. 2025).
   For each (assay, candidate feature) correlate the feature's per-residue WT activation a[p]
   against the per-residue mean mutation effect m[p] (mean DMS_score over substitutions at p,
   positions with >= MIN_SUBS substitutions). Spearman across positions; best feature per assay.
   This measures "does the feature mark the residues where mutations matter" (a position-importance
   readout), the same quantity MotifAE reports (their median 0.41 / plain-SAE 0.33, ESM2).

B) Held-out-split feature selection (InterPLM / InterProt convention).
   The literature guard against best-of-N inflation is a train/test split, not a permutation null.
   Per assay we split variants by POSITION (position-grouped, to avoid the position-conservation
   leak), select the feature with the highest |Spearman(pool_mean_abs, DMS)| on the TRAIN variants,
   and report that feature's |Spearman| on the held-out TEST variants. If the per-assay best is
   selection overfitting, test |rho| collapses toward 0.

Both run on data/proteingym/pooled_metrics/*__pooled.parquet (the InterPLM-paired candidate union
per assay). Single-mutant variants only.

Usage (from the crosscoder root):
  uv run --with pandas --with numpy --with scipy --with pyarrow \
    python repos/sparse-crosscoders-prott5/per_residue_and_holdout.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path("/Users/sohrab.tawana/private/crosscoder")
POOLED = ROOT / "data/proteingym/pooled_metrics"
OUT = ROOT / "data/proteingym/per_residue_holdout_summary.csv"
MIN_SUBS = 10      # MotifAE: residues with >=10 assayed substitutions
MIN_POS = 5        # require >=5 usable positions for a per-residue correlation
MIN_FOLD = 15      # require >=15 variants in each split fold
SEED = 0


def main() -> None:
    rng = np.random.default_rng(SEED)
    rows = []
    for p in sorted(POOLED.glob("*__pooled.parquet")):
        dms_id = p.stem.replace("__pooled", "")
        df = pd.read_parquet(p)
        df = df[df["n_mutations"] == 1].copy()
        if len(df) == 0:
            continue

        # ---- A) per-residue MotifAE-comparable ----
        per_res_best = np.nan
        # mean DMS per position (over substitutions); positions with >= MIN_SUBS variants
        var = df.drop_duplicates(["mutant"])  # one row per variant for the position->effect map
        pos_effect = var.groupby("position")["DMS_score"].agg(["mean", "size"])
        good_pos = pos_effect.index[pos_effect["size"] >= MIN_SUBS]
        if len(good_pos) >= MIN_POS:
            m = pos_effect.loc[good_pos, "mean"]
            best = 0.0
            for feat, g in df.groupby("feature"):
                a = g.groupby("position")["at_pos_wt"].first().reindex(good_pos)
                ok = a.notna() & m.notna()
                if ok.sum() < MIN_POS or a[ok].std() < 1e-9:
                    continue
                rho = spearmanr(a[ok], m[ok]).statistic
                if np.isfinite(rho):
                    best = max(best, abs(rho))
            per_res_best = best if best > 0 else np.nan

        # ---- B) position-grouped held-out split on pool_mean_abs ----
        train_sel, test_sel = np.nan, np.nan
        positions = np.array(sorted(df["position"].unique()))
        if len(positions) >= 4:
            rng.shuffle(positions)
            cut = len(positions) // 2
            train_pos, test_pos = set(positions[:cut]), set(positions[cut:])
            tr = df[df["position"].isin(train_pos)]
            te = df[df["position"].isin(test_pos)]
            if tr["mutant"].nunique() >= MIN_FOLD and te["mutant"].nunique() >= MIN_FOLD:
                best_feat, best_tr = None, -1.0
                for feat, g in tr.groupby("feature"):
                    if g["pool_mean_abs"].std() < 1e-9 or g["DMS_score"].std() < 1e-9:
                        continue
                    r = abs(spearmanr(g["pool_mean_abs"], g["DMS_score"]).statistic)
                    if np.isfinite(r) and r > best_tr:
                        best_tr, best_feat = r, feat
                if best_feat is not None:
                    gte = te[te["feature"] == best_feat]
                    if gte["pool_mean_abs"].std() > 1e-9 and gte["DMS_score"].std() > 1e-9:
                        rte = spearmanr(gte["pool_mean_abs"], gte["DMS_score"]).statistic
                        train_sel, test_sel = best_tr, abs(rte) if np.isfinite(rte) else np.nan

        rows.append({"DMS_id": dms_id, "per_residue_best": per_res_best,
                     "holdout_train_sel": train_sel, "holdout_test_sel": test_sel})

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT} ({len(out)} assays)\n")
    print("A) per-residue MotifAE-comparable (best paired feature per assay):")
    v = out["per_residue_best"].dropna()
    print(f"   median |rho| = {v.median():.3f}  (n={len(v)} assays, >0.3={int((v>0.3).sum())}, >0.5={int((v>0.5).sum())})")
    print(f"   [MotifAE per-residue best, ESM2: 0.41; plain SAE: 0.33]\n")
    print("B) position-grouped held-out split on pool_mean_abs:")
    tr = out["holdout_train_sel"].dropna(); te = out["holdout_test_sel"].dropna()
    print(f"   train-selected |rho| (in-sample, inflated) median = {tr.median():.3f}  (n={len(tr)})")
    print(f"   held-out TEST |rho| of train-selected feature median = {te.median():.3f}  (>0.3={int((te>0.3).sum())}, >0.5={int((te>0.5).sum())})")


if __name__ == "__main__":
    main()
