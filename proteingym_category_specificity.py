#!/usr/bin/env python3
"""ProteinGym category-specificity test (experiment 03).

Claim under test: a feature paired (InterPLM) to a concept tracks the DMS *functional category*
that matches the concept's biology more than it tracks the other categories. The position/AA
prior that confounds the absolute ProteinGym correlation (eta=0.575) is category-blind, so a
matched-vs-mismatched gap is signal the prior cannot manufacture.

Design
------
- Unit = one (assay, concept) pairing. Metric = MEAN |pool_mean_abs| over the concept's features
  for that pairing (selection-free: no best-of-N, so the matched/mismatched contrast is unbiased).
  `pool_mean_abs` is the canonical pooled Spearman (signed, [-1,1]); we use |.|.
- A-priori concept -> category map (textbook biology, fixed before looking at correlations; see
  CONCEPT_CATEGORY). Region_Disordered and the zinc-finger / J-domain concepts are excluded as
  having no unambiguous functional category. The two glycosylation labels are merged.
- Restrict to assays whose coarse_selection_type is one of the 4 specific categories
  {Stability, Activity, Expression, Binding}. OrganismalFitness is a composite growth phenotype
  (not concept-specific) and is excluded.
- matched := (assay.category == concept.category).
- Tests: (1) Mann-Whitney U on |rho| matched vs mismatched (one-sided: matched > mismatched);
  (2) label-permutation null on the mean gap = mean(matched) - mean(mismatched): randomly
  reassign each concept to one of the 4 categories, recompute the gap, N times; p = frac >= obs.

Null choice: the baseline crosscoder is NOT a usable control here (it has ~1 paired concept and no
pooled_metrics_baseline/), so specificity has no baseline pairings to test. The label-permutation
is the principled self-contained null. (A baseline pooled harvest remains a documented open item.)

Inputs (read-only):
  - data/proteingym/pooled_metrics/summary.csv   (our scored (DMS_id, concept, feature) rows)
  - data/external/DMS_substitutions.csv          (ProteinGym reference; coarse_selection_type)
Outputs:
  - data/proteingym/category_specificity_pairings.csv   (per-pairing matched/mismatched + |rho|)
  - data/proteingym/category_specificity_summary.csv     (per-concept + headline numbers)
  - prints the full result to stdout

Seed: 0 (only used for the permutation null). Reproducible.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

SEED = 0
N_PERM = 10000
SPECIFIC = ["Stability", "Activity", "Expression", "Binding"]

# A-priori concept -> matched functional category (fixed before inspecting any correlation).
CONCEPT_CATEGORY = {
    "Disulfide bond": "Stability",
    "Coiled coil": "Stability",
    "Glycosylation": "Expression",            # both glyco labels merged below
    "Domain_SH3": "Binding",
    "Domain_PDZ": "Binding",
    "Domain_NR LBD": "Binding",
    "Domain_C-type lectin": "Binding",
    "Domain_Protein kinase": "Activity",
    "Domain_Nudix hydrolase": "Activity",
    "Domain_N-acetyltransferase": "Activity",
    "Domain_Helicase C-terminal": "Activity",
    "Domain_Helicase ATP-binding": "Activity",
    "Motif_Nudix box": "Activity",
    "Modified residue_Phosphohistidine": "Activity",
}
# Excluded concepts (no unambiguous functional category): Region_Disordered, Zinc finger_any,
# Zinc finger_RING-type, Domain_J.

ROOT = Path(__file__).resolve().parents[2]
if not (ROOT / "data").exists():
    ROOT = Path.cwd()
SUMMARY = ROOT / "data" / "proteingym" / "pooled_metrics" / "summary.csv"
REF = ROOT / "data" / "external" / "DMS_substitutions.csv"
OUT_PAIR = ROOT / "data" / "proteingym" / "category_specificity_pairings.csv"
OUT_SUM = ROOT / "data" / "proteingym" / "category_specificity_summary.csv"


def main():
    ref = pd.read_csv(REF, usecols=["DMS_id", "coarse_selection_type"])
    ref["coarse_selection_type"] = ref["coarse_selection_type"].astype(str).str.strip()
    cat_by_dms = dict(zip(ref["DMS_id"], ref["coarse_selection_type"]))

    df = pd.read_csv(SUMMARY, usecols=["DMS_id", "concept", "feature", "pool_mean_abs"])
    df["concept"] = df["concept"].astype(str).str.strip()
    # merge the two near-duplicate glycosylation labels into one concept
    df.loc[df["concept"].str.startswith("Glycosylation"), "concept"] = "Glycosylation"

    df["category"] = df["DMS_id"].map(cat_by_dms)
    df["abs_rho"] = df["pool_mean_abs"].abs()
    df = df.dropna(subset=["abs_rho"])  # drop features with undefined Spearman (constant activation)

    # keep only mapped concepts and the 4 specific categories
    df = df[df["concept"].isin(CONCEPT_CATEGORY)]
    df = df[df["category"].isin(SPECIFIC)]

    # one row per (assay, concept) pairing, under two aggregations of the concept's features:
    #   mean  = selection-free (no best-of-N);  best = max |rho| (best-of-N, symmetric across groups)
    pair = (df.groupby(["DMS_id", "concept", "category"])
              .agg(abs_rho_mean=("abs_rho", "mean"),
                   abs_rho_best=("abs_rho", "max"),
                   n_features=("feature", "nunique"))
              .reset_index())
    pair["concept_category"] = pair["concept"].map(CONCEPT_CATEGORY)
    pair["matched"] = pair["category"] == pair["concept_category"]
    pair.to_csv(OUT_PAIR, index=False)

    def run(metric_col):
        m = pair.loc[pair["matched"], metric_col].to_numpy()
        mm = pair.loc[~pair["matched"], metric_col].to_numpy()
        obs_gap = m.mean() - mm.mean()
        u_stat, u_p = (mannwhitneyu(m, mm, alternative="greater") if len(m) and len(mm)
                       else (float("nan"), float("nan")))
        rng = np.random.default_rng(SEED)
        concepts = pair["concept"].unique()
        cat_arr = pair["category"].to_numpy()
        rho_arr = pair[metric_col].to_numpy()
        concept_idx = pair["concept"].to_numpy()
        perm_gaps = np.empty(N_PERM)
        for i in range(N_PERM):
            assign = {c: rng.choice(SPECIFIC) for c in concepts}
            matched = np.array([cat_arr[j] == assign[concept_idx[j]] for j in range(len(rho_arr))])
            perm_gaps[i] = (0.0 if matched.all() or (~matched).all()
                            else rho_arr[matched].mean() - rho_arr[~matched].mean())
        perm_p = (1 + np.sum(perm_gaps >= obs_gap)) / (N_PERM + 1)
        return dict(n_m=len(m), n_mm=len(mm), mean_m=m.mean(), mean_mm=mm.mean(),
                    gap=obs_gap, u=u_stat, u_p=u_p, perm_p=perm_p,
                    perm_mean=perm_gaps.mean(), perm_p95=np.percentile(perm_gaps, 95))

    res_mean = run("abs_rho_mean")
    res_best = run("abs_rho_best")

    # per-concept descriptive: matched-cell mean vs that concept's mismatched mean
    rows = []
    for c in sorted(pair["concept"].unique()):
        sub = pair[pair["concept"] == c]
        cm = sub.loc[sub["matched"], "abs_rho_mean"]
        cmm = sub.loc[~sub["matched"], "abs_rho_mean"]
        rows.append({
            "concept": c,
            "matched_category": CONCEPT_CATEGORY[c],
            "n_matched_assays": len(cm),
            "n_mismatched_assays": len(cmm),
            "matched_mean_absrho": round(cm.mean(), 4) if len(cm) else np.nan,
            "mismatched_mean_absrho": round(cmm.mean(), 4) if len(cmm) else np.nan,
        })
    per_concept = pd.DataFrame(rows)
    per_concept.to_csv(OUT_SUM, index=False)

    pd.set_option("display.width", 170, "display.max_columns", 20)
    print("=== ProteinGym category-specificity test ===")
    print(f"seed={SEED}  N_perm={N_PERM}  metric=mean|pool_mean_abs| per (assay,concept) pairing\n")
    print("Per-concept matched vs mismatched mean |rho|:")
    print(per_concept.to_string(index=False))
    for name, r in [("mean over features (selection-free)", res_mean),
                    ("best feature per pairing (symmetric best-of-N)", res_best)]:
        print(f"\n--- aggregation: {name} ---")
        print(f"  matched n={r['n_m']} mean|rho|={r['mean_m']:.4f}   |   "
              f"mismatched n={r['n_mm']} mean|rho|={r['mean_mm']:.4f}")
        print(f"  gap (matched - mismatched) = {r['gap']:+.4f}")
        print(f"  Mann-Whitney (matched > mismatched): U={r['u']:.1f}, one-sided p={r['u_p']:.4f}")
        print(f"  label-permutation null: p={r['perm_p']:.4f} "
              f"(null mean={r['perm_mean']:+.4f}, p95={r['perm_p95']:+.4f})")
    print(f"\nWrote {OUT_PAIR.name}, {OUT_SUM.name}")


if __name__ == "__main__":
    main()
