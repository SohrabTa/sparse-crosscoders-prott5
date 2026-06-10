"""
Qualitative analysis of the ProteinGym (feature, assay) pooled-metric correlations.

Goal (experiment 03, design "(b)"): instead of asking "does the activation readout
validate the labels" (it does not — see 03-proteingym.md), characterise *where* the
(feature, assay) signal is strong vs weak, to understand the metric's strengths/flaws
and to scope whether the concept-localized fitness test (design "(1)") is worth running.

We slice the per-(assay, concept, feature) pooled-metric Spearman by:
  - concept structural locality  (point/site vs domain vs disordered/extended)
  - assay functional category     (Stability / Activity / Binding / Expression / OrganismalFitness)
  - concept↔category match         (does the feature's concept sit on the assay's functional axis)
  - pairing strength               (f1_per_domain)   -- known weak +0.27 trend
  - firing sparsity                (fire_rate_mean)  -- concept specialists fire sparsely

INPUTS  (read-only):
  data/proteingym/pooled_metrics/summary.csv   (job 5664029; 2678 rows, per (DMS_id, concept, feature))
      cols used: DMS_id, concept, feature, n_variants, fire_rate_mean, f1_per_domain,
                 pool_mean_abs, pool_mean_abs_resid (position-residualized = substitution-specific)
  data/proteingym/concept_matches.csv          (match_proteingym_concepts.py)
      cols used: DMS_id, matched_concept, coarse_selection_type, target_seq_len,
                 n_concept_residues, frac_variants_in_concept

OUTPUTS (written under data/proteingym/qualitative/):
  tidy.csv                         -- merged per-(assay,concept,feature) table + derived columns
  by_locality.csv                  -- |rho| distribution by concept locality class
  by_category.csv                  -- |rho| distribution by assay functional category
  by_match.csv                     -- matched vs mismatched concept/category
  drivers.csv                      -- Spearman of |rho| against f1_per_domain & fire_rate, overall + per locality
  top_bottom.csv                   -- the 25 strongest and 25 weakest (assay,concept,feature) tuples + properties

Aggregation mirrors the experiment doc: signal = |pool_mean_abs| (and the residualized
substitution-specific |pool_mean_abs_resid|). Bootstrap 95% CIs on group medians use SEED.

Disposition: COMMIT (backs a reported qualitative result in 03-proteingym.md).

Usage:
    uv run --with pandas --with numpy --with scipy python \
        repos/sparse-crosscoders-prott5/qualitative_proteingym_analysis.py
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SEED = 42
N_BOOT = 2000

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("qual")

ROOT = Path(__file__).resolve().parents[2]            # crosscoder root
PG = ROOT / "data" / "proteingym"
OUT = PG / "qualitative"
OUT.mkdir(exist_ok=True)


def locality_class(concept: str) -> str:
    """Structural locality of a Swiss-Prot concept: how spatially concentrated its
    residues are. Point/site concepts have a 3D locus (where design (1) should work);
    disordered/extended concepts do not."""
    c = concept.lower()
    if c.startswith(("disulfide", "glycosylation", "modified residue", "zinc finger", "motif")):
        return "point_site"
    if c.startswith("domain"):
        return "domain"
    if c.startswith("region") or c.startswith("coiled coil"):
        return "disordered_extended"
    return "other"


def boot_median_ci(x: np.ndarray, rng: np.random.Generator) -> tuple[float, float, float]:
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    if len(x) == 0:
        return (np.nan, np.nan, np.nan)
    med = float(np.median(x))
    if len(x) < 3:
        return (med, np.nan, np.nan)
    boots = [np.median(rng.choice(x, size=len(x), replace=True)) for _ in range(N_BOOT)]
    return (med, float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)))


def main() -> None:
    rng = np.random.default_rng(SEED)
    log.info("SEED=%d  N_BOOT=%d", SEED, N_BOOT)

    summ = pd.read_csv(PG / "pooled_metrics" / "summary.csv")
    matches = pd.read_csv(PG / "concept_matches.csv")
    log.info("summary rows=%d  matches rows=%d", len(summ), len(matches))

    meta = matches[
        ["DMS_id", "matched_concept", "coarse_selection_type",
         "target_seq_len", "n_concept_residues", "frac_variants_in_concept"]
    ].rename(columns={"matched_concept": "concept", "coarse_selection_type": "assay_category"})
    meta = meta.drop_duplicates(["DMS_id", "concept"])

    df = summ.merge(meta, on=["DMS_id", "concept"], how="left")
    df["abs_rho"] = df["pool_mean_abs"].abs()
    df["abs_rho_resid"] = df["pool_mean_abs_resid"].abs()
    df["locality"] = df["concept"].map(locality_class)

    # concept's "home" functional axis -> does it match this assay's category?
    concept_axis = {
        "Disulfide bond": "Stability", "Coiled coil": "Stability",
        "Zinc finger_RING-type": "Binding", "Zinc finger_any": "Binding",
        "Glycosylation_any": "Expression",
        "Glycosylation_N-linked (GlcNAc...) asparagine": "Expression",
        "Domain_Protein kinase": "Activity", "Domain_Nudix hydrolase": "Activity",
        "Domain_NR LBD": "Binding", "Domain_SH3": "Binding", "Domain_PDZ": "Binding",
        "Motif_Nudix box": "Activity", "Modified residue_Phosphohistidine": "Activity",
    }
    df["concept_axis"] = df["concept"].map(concept_axis)
    df["match"] = np.where(
        df["concept_axis"].isna(), "unknown",
        np.where(df["concept_axis"] == df["assay_category"], "matched", "mismatched"),
    )

    df.to_csv(OUT / "tidy.csv", index=False)
    log.info("wrote tidy.csv (%d rows)", len(df))

    def grouped(col_group: str) -> pd.DataFrame:
        rows = []
        for key, g in df.groupby(col_group):
            med, lo, hi = boot_median_ci(g["abs_rho"].values, rng)
            medr, _, _ = boot_median_ci(g["abs_rho_resid"].values, rng)
            rows.append({
                col_group: key, "n_rows": len(g),
                "n_assays": g["DMS_id"].nunique(),
                "median_abs_rho": med, "ci_lo": lo, "ci_hi": hi,
                "median_abs_rho_resid": medr,
                "frac_gt_0.3": float((g["abs_rho"] > 0.3).mean()),
                "frac_gt_0.5": float((g["abs_rho"] > 0.5).mean()),
            })
        return pd.DataFrame(rows).sort_values("median_abs_rho", ascending=False)

    by_loc = grouped("locality")
    by_cat = grouped("assay_category")
    by_match = grouped("match")
    by_loc.to_csv(OUT / "by_locality.csv", index=False)
    by_cat.to_csv(OUT / "by_category.csv", index=False)
    by_match.to_csv(OUT / "by_match.csv", index=False)

    # drivers: does pairing strength / sparsity predict the signal?
    drv = []
    for scope, g in [("ALL", df)] + [(f"loc={k}", v) for k, v in df.groupby("locality")]:
        for driver in ["f1_per_domain", "fire_rate_mean"]:
            sub = g[[driver, "abs_rho"]].dropna()
            if len(sub) >= 10:
                rho, p = spearmanr(sub[driver], sub["abs_rho"])
            else:
                rho, p = (np.nan, np.nan)
            drv.append({"scope": scope, "driver": driver, "n": len(sub),
                        "spearman_vs_abs_rho": rho, "p": p})
    pd.DataFrame(drv).to_csv(OUT / "drivers.csv", index=False)

    cols = ["DMS_id", "concept", "feature", "locality", "assay_category", "match",
            "abs_rho", "abs_rho_resid", "f1_per_domain", "fire_rate_mean",
            "n_variants", "n_concept_residues", "frac_variants_in_concept", "target_seq_len"]
    ranked = df.dropna(subset=["abs_rho"]).sort_values("abs_rho", ascending=False)
    top_bottom = pd.concat([ranked.head(25).assign(end="top"),
                            ranked.tail(25).assign(end="bottom")])
    top_bottom[cols + ["end"]].to_csv(OUT / "top_bottom.csv", index=False)

    # --- console digest ---
    log.info("\n=== by locality (|pool_mean_abs|) ===\n%s", by_loc.to_string(index=False))
    log.info("\n=== by assay category ===\n%s", by_cat.to_string(index=False))
    log.info("\n=== matched vs mismatched concept axis ===\n%s", by_match.to_string(index=False))
    log.info("\n=== drivers (Spearman of |rho| vs pairing strength / sparsity) ===\n%s",
             pd.DataFrame(drv).to_string(index=False))


if __name__ == "__main__":
    main()
