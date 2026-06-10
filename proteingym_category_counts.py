#!/usr/bin/env python3
"""Feasibility table for the ProteinGym category-specificity test.

Question it supports (experiment 03): does a feature paired to a concept track the
DMS *functional category* that matches its biology (e.g. a disulfide-bond feature -> Stability
assays) more than the other categories? Before running that test we need to know it has the
statistical power: how many *distinct assays* fall in each (concept x coarse_selection_type) cell.

Inputs (read-only):
  - data/DMS_substitutions.csv                         (ProteinGym reference; has coarse_selection_type)
  - data/proteingym/pooled_metrics/summary.csv         (our (DMS_id, concept, feature) scored rows)

Output:
  - data/proteingym/category_counts.csv                (concept x category matrix, distinct-assay counts)
  - prints the same matrix + per-concept matched-category coverage to stdout

No randomness; no seed needed. One-off feasibility check.
"""
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]  # crosscoder root (…/repos/sparse-crosscoders-prott5 -> up 2)
if not (ROOT / "data").exists():
    # fall back: allow running from crosscoder root directly
    ROOT = Path.cwd()

REF = ROOT / "data" / "DMS_substitutions.csv"
SUMMARY = ROOT / "data" / "proteingym" / "pooled_metrics" / "summary.csv"
OUT = ROOT / "data" / "proteingym" / "category_counts.csv"

ref = pd.read_csv(REF, usecols=["DMS_id", "coarse_selection_type"])
ref["coarse_selection_type"] = ref["coarse_selection_type"].astype(str).str.strip()
cat_by_dms = dict(zip(ref["DMS_id"], ref["coarse_selection_type"]))

summ = pd.read_csv(SUMMARY, usecols=["DMS_id", "concept"])
summ["concept"] = summ["concept"].astype(str).str.strip()
summ = summ.drop_duplicates(["DMS_id", "concept"])  # one row per (assay, concept) pairing
summ["category"] = summ["DMS_id"].map(cat_by_dms)

missing = summ[summ["category"].isna()]["DMS_id"].unique()
if len(missing):
    print(f"WARN: {len(missing)} scored DMS_ids not found in reference: {sorted(missing)[:5]}…",
          file=sys.stderr)
summ = summ.dropna(subset=["category"])

# distinct assays per (concept x category)
mat = (summ.groupby(["concept", "category"])["DMS_id"]
       .nunique().unstack(fill_value=0))
# order columns by total assay count in the benchmark
col_order = ["OrganismalFitness", "Stability", "Activity", "Expression", "Binding"]
mat = mat.reindex(columns=[c for c in col_order if c in mat.columns], fill_value=0)
mat["TOTAL_assays"] = summ.groupby("concept")["DMS_id"].nunique()
mat = mat.sort_values("TOTAL_assays", ascending=False)

OUT.parent.mkdir(parents=True, exist_ok=True)
mat.to_csv(OUT)

pd.set_option("display.width", 160, "display.max_columns", 20)
print(f"Distinct assays per (concept x coarse_selection_type)  ->  {OUT}\n")
print(mat.to_string())
print(f"\nTotal distinct assays covered: {summ['DMS_id'].nunique()}")
print(f"Total (assay, concept) pairings: {len(summ)}")
