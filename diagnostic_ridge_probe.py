"""
Ridge-probe diagnostic battery for the ProteinGym × Crosscoder evaluation.

Consumes:
  - data/proteingym/raw_prott5_layer12/<DMS_id>.npz   (output of harvest_raw_prott5.py)
  - data/proteingym/scoring/<DMS_id>__<concept>__per_variant.parquet (existing)
  - data/proteingym_concept_matches.csv

Produces per-assay (and per-(assay,concept)) ridge-probe Spearman scores for:

  A) raw_mean_pool      : ridge( emb_mut_mean → DMS )
  B) raw_pos            : ridge( emb_mut_at_pos → DMS )
  C) raw_pos_delta      : ridge( emb_mut_at_pos − emb_wt_at_pos → DMS )
  D) raw_pool_delta     : ridge( emb_mut_mean − emb_wt_mean → DMS )
  E) sae_candidates     : ridge( act_mut[candidate features] → DMS )   per (assay,concept)
  F) sae_candidates_Δ   : ridge( delta_act[candidate features] → DMS ) per (assay,concept)

Each probe is evaluated with 5-fold CV using held-out Spearman.

The headline question: does the raw-ProtT5 ridge probe (A/B/C/D) get a
*meaningfully higher* Spearman than our existing Δactivation correlation on
the same assays? If yes, there is real signal we are leaving on the table
and the SAE / metric is the bottleneck. If no, even raw ProtT5 doesn't
encode DMS at this layer and we need a different instrument entirely
(pseudo-log-likelihood, supervised probes like VESPA, etc.).

Usage:
    uv run --project repos/crosscode --with scikit-learn python \\
        repos/sparse-crosscoders-prott5/diagnostic_ridge_probe.py \\
        --raw_dir data/proteingym/raw_prott5_layer12 \\
        --scoring_dir data/proteingym/scoring \\
        --matches data/proteingym_concept_matches.csv \\
        --output data/proteingym/diagnostic_summary.csv
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# scikit-learn for ridge + KFold
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("diag")


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    x, y = x[m], y[m]
    if len(x) < 3:
        return float("nan")
    sx = pd.Series(x); sy = pd.Series(y)
    if sx.std() == 0 or sy.std() == 0:
        return float("nan")
    return float(np.corrcoef(sx.rank(), sy.rank())[0, 1])


def cv_ridge_spearman(
    X: np.ndarray, y: np.ndarray, n_splits: int = 5, seed: int = 42
) -> dict:
    """Returns dict with held-out Spearman (mean/median across folds) +
    diagnostics. X: (N, D), y: (N,). Uses sklearn RidgeCV with sensible alphas."""
    N = len(y)
    if N < 20:
        return {"n": N, "spearman_cv": float("nan"), "ok": False, "reason": "n<20"}
    if X.shape[1] == 0:
        return {"n": N, "spearman_cv": float("nan"), "ok": False, "reason": "D=0"}
    # if y is constant or near-constant, bail
    if np.std(y) < 1e-10:
        return {"n": N, "spearman_cv": float("nan"), "ok": False, "reason": "y_const"}

    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    preds = np.zeros_like(y, dtype=float)
    alphas = (1.0, 10.0, 100.0, 1000.0, 10000.0)
    for fold_idx, (tr, te) in enumerate(kf.split(X)):
        # standardize features on train, apply to test
        mu = X[tr].mean(axis=0); sigma = X[tr].std(axis=0)
        sigma = np.where(sigma < 1e-8, 1.0, sigma)
        Xtr = (X[tr] - mu) / sigma
        Xte = (X[te] - mu) / sigma
        try:
            model = RidgeCV(alphas=alphas, cv=3 if len(tr) > 60 else None).fit(Xtr, y[tr])
        except Exception:
            preds[te] = y[tr].mean()
            continue
        preds[te] = model.predict(Xte)
    rho = spearman(preds, y)
    return {"n": N, "spearman_cv": rho, "ok": np.isfinite(rho), "reason": ""}


def diag_raw_one_assay(npz_path: Path) -> dict:
    """Run diagnostics A/B/C/D on the raw ProtT5 npz for one assay."""
    z = np.load(npz_path, allow_pickle=True)
    y = z["DMS_score"]
    emb_wt_mean = z["emb_wt_mean"]                 # (D,)
    emb_wt_at_pos = z["emb_wt_at_pos"]             # (N, D)
    emb_mut_mean = z["emb_mut_mean"]               # (N, D)
    emb_mut_at_pos = z["emb_mut_at_pos"]           # (N, D)
    N = len(y)

    out = {"DMS_id": npz_path.stem, "N_variants": N,
           "L_wt": int(z["seq_len_wt"]), "layer": int(z["layer"])}

    feats = {
        "A_raw_mean_pool":   emb_mut_mean,
        "B_raw_pos":         emb_mut_at_pos,
        "C_raw_pos_delta":   emb_mut_at_pos - emb_wt_at_pos,
        "D_raw_pool_delta":  emb_mut_mean   - emb_wt_mean[None, :],
    }
    for name, X in feats.items():
        r = cv_ridge_spearman(X, y)
        out[name] = r["spearman_cv"]
        out[f"{name}_n"] = r["n"]
    return out


def load_scoring_parquet(parquet_path: Path) -> pd.DataFrame:
    df = pd.read_parquet(parquet_path)
    return df


def diag_sae_one_pair(
    parquet_path: Path, dms_id: str, concept: str
) -> dict:
    """Run diagnostics E/F on the SAE candidate-feature activations for one
    (assay, concept) pair, using the existing per_variant.parquet."""
    df = load_scoring_parquet(parquet_path)
    if df.empty:
        return {"DMS_id": dms_id, "concept": concept, "ok": False, "reason": "empty"}

    # Reshape to wide: rows=variants, cols=features. Use mean over multi-mutant rows.
    wide_mut = df.pivot_table(
        index="mutant", columns="feature", values="act_mut", aggfunc="mean"
    ).fillna(0.0)
    wide_delta = df.pivot_table(
        index="mutant", columns="feature", values="delta_act", aggfunc="mean"
    ).fillna(0.0)
    # DMS_score is constant per mutant — take first
    y = df.drop_duplicates("mutant").set_index("mutant")["DMS_score"].reindex(wide_mut.index).values
    X_mut = wide_mut.values
    X_del = wide_delta.values
    n_feat = X_mut.shape[1]

    out = {"DMS_id": dms_id, "concept": concept,
           "n_candidate_features": int(n_feat), "N_variants": len(y)}
    rE = cv_ridge_spearman(X_mut, y); out["E_sae_act"]   = rE["spearman_cv"]
    rF = cv_ridge_spearman(X_del, y); out["F_sae_delta"] = rF["spearman_cv"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw_dir", type=Path, required=True,
                    help="Dir of raw ProtT5 .npz files (output of harvest_raw_prott5.py)")
    ap.add_argument("--scoring_dir", type=Path, required=True,
                    help="Dir of existing per_variant.parquet scoring outputs")
    ap.add_argument("--matches", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    matches = pd.read_csv(args.matches)
    # only consider assays where we have BOTH raw npz and scoring parquets
    raw_files = {p.stem: p for p in args.raw_dir.glob("*.npz")}
    log.info("Found %d raw_prott5 .npz files", len(raw_files))

    # ---- Diagnostic A/B/C/D : raw-ProtT5 ridge probe per assay ----
    raw_rows = []
    for dms_id, npz in sorted(raw_files.items()):
        try:
            row = diag_raw_one_assay(npz)
            raw_rows.append(row)
        except Exception as e:
            log.warning("  %s raw failed: %s", dms_id, e)
    raw_df = pd.DataFrame(raw_rows)
    log.info("Raw probe done: %d assays scored", len(raw_df))

    # ---- Diagnostic E/F : SAE candidate-feature ridge probe per (assay,concept) ----
    sae_rows = []
    for _, m in matches.iterrows():
        dms_id = m["DMS_id"]; concept = m["matched_concept"]
        if dms_id not in raw_files:
            continue
        safe_concept = concept.replace("/", "_").replace(" ", "_")
        parquet = args.scoring_dir / f"{dms_id}__{safe_concept}__per_variant.parquet"
        if not parquet.exists():
            continue
        try:
            row = diag_sae_one_pair(parquet, dms_id, concept)
            sae_rows.append(row)
        except Exception as e:
            log.warning("  %s|%s sae failed: %s", dms_id, concept, e)
    sae_df = pd.DataFrame(sae_rows)
    log.info("SAE probe done: %d (assay,concept) pairs scored", len(sae_df))

    # ---- Merge to a single per-assay summary table ----
    if not sae_df.empty:
        sae_best = sae_df.assign(
            E_abs=lambda d: d["E_sae_act"].abs(),
            F_abs=lambda d: d["F_sae_delta"].abs(),
        )
        best_per_assay = sae_best.sort_values("E_abs", ascending=False).drop_duplicates("DMS_id")
        best_per_assay = best_per_assay[["DMS_id", "concept", "E_sae_act", "F_sae_delta",
                                          "n_candidate_features", "N_variants"]]
        best_per_assay = best_per_assay.rename(columns={
            "concept": "best_concept_for_sae",
            "E_sae_act": "E_sae_act_best",
            "F_sae_delta": "F_sae_delta_best",
            "n_candidate_features": "sae_n_features",
            "N_variants": "sae_N_variants",
        })
    else:
        best_per_assay = pd.DataFrame(columns=["DMS_id"])

    merged = raw_df.merge(best_per_assay, on="DMS_id", how="left")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, index=False)
    sae_df.to_csv(args.output.with_suffix("").with_name(args.output.stem + "_sae_all.csv"), index=False)
    log.info("Wrote per-assay summary -> %s", args.output)

    # ---- Headline aggregate ----
    print("\n" + "=" * 72)
    print(" HEADLINE — Held-out 5-fold Spearman across assays")
    print("=" * 72)
    cols = ["A_raw_mean_pool", "B_raw_pos", "C_raw_pos_delta", "D_raw_pool_delta",
            "E_sae_act_best", "F_sae_delta_best"]
    rows = []
    for c in cols:
        if c in merged.columns:
            v = merged[c].dropna()
            v_abs = v.abs() if c.startswith(("E_", "F_")) else v  # SAE we take |.|
            rows.append({
                "probe": c,
                "n_assays": int(v.notna().sum()),
                "median": float(v_abs.median()) if len(v_abs) else float("nan"),
                "mean":   float(v_abs.mean())   if len(v_abs) else float("nan"),
                "p25":    float(v_abs.quantile(.25)) if len(v_abs) else float("nan"),
                "p75":    float(v_abs.quantile(.75)) if len(v_abs) else float("nan"),
                "frac_gt_0.30": float((v_abs > .30).mean()) if len(v_abs) else float("nan"),
                "frac_gt_0.50": float((v_abs > .50).mean()) if len(v_abs) else float("nan"),
            })
    head = pd.DataFrame(rows)
    print(head.to_string(index=False))

    # ---- Diff: ceiling gap ----
    if not merged.empty and "B_raw_pos" in merged and "E_sae_act_best" in merged:
        d = merged.dropna(subset=["B_raw_pos", "E_sae_act_best"]).copy()
        d["gap_B_vs_E"] = d["B_raw_pos"].abs() - d["E_sae_act_best"].abs()
        d["gap_C_vs_F"] = d["C_raw_pos_delta"].abs() - d["F_sae_delta_best"].abs()
        print("\nCEILING GAP (raw_pos − best_SAE_per_assay)")
        for c in ["gap_B_vs_E", "gap_C_vs_F"]:
            v = d[c].dropna()
            if not len(v): continue
            print(f"  {c}: n={len(v)}  median={v.median():+.3f}  mean={v.mean():+.3f}  "
                  f"p25={v.quantile(.25):+.3f}  p75={v.quantile(.75):+.3f}  "
                  f"frac>0={(v>0).mean()*100:.1f}%")


if __name__ == "__main__":
    import sys; sys.exit(main())
