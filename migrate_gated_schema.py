"""One-time schema migration: gated-only ProteinGym data.

Every crosscoder feature activation we use is gated (post-BatchTopK), so the "dense pre-gate"
(`d_*`) readout was removed as invalid, and the now-redundant `g_` prefix / "gated" qualifier is
dropped from column names. This script rewrites the committed ProteinGym data in place:

  - drop every `d_*` column (dense pre-gate readout)
  - strip the `g_` prefix from every gated column   (g_pool_mean_abs -> pool_mean_abs)
  - rename fire_rate_gated_{mean,var} -> fire_rate_{mean,var}

Files migrated (CSVs + per-assay parquets):
  data/proteingym/pooled_metrics/summary.csv               (+ residualized columns merged in)
  data/proteingym/pooled_metrics/*__pooled.parquet
  data/proteingym/position_residualized_summary.csv
  data/proteingym/full_feature_spearman/full_feature_spearman.csv
  data/proteingym/full_feature_spearman/*__full_feat.parquet

It also merges the residualized gated columns (`*_resid`) from
position_residualized_summary.csv into pooled_metrics/summary.csv, keyed on
(DMS_id, concept, feature), so the residualized robustness signal lives alongside the raw signal.

(The leaderboard comparison.csv / pseudoll artifacts were since removed entirely — the ProteinGym
analysis no longer compares to fitness predictors.)

Idempotent: a file already lacking `g_`/`d_` columns is left unchanged.

Usage (from the crosscoder root):
    uv run --with pandas --with pyarrow python repos/sparse-crosscoders-prott5/migrate_gated_schema.py
"""
from __future__ import annotations

import glob
from pathlib import Path

import pandas as pd

ROOT = Path("/Users/sohrab.tawana/private/crosscoder")
PG = ROOT / "data/proteingym"

# residualized metric columns to merge into summary.csv (post-rename names)
RESID_COLS = [
    "at_pos_mut_resid", "delta_at_pos_resid", "pool_mean_abs_resid",
    "pool_max_abs_resid", "pool_mean_mut_resid", "pool_max_mut_resid",
]


def rename_col(c: str) -> str:
    c = c.replace("fire_rate_gated", "fire_rate")
    if c.startswith("g_"):
        return c[2:]
    return c


def migrate_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Drop d_* columns, strip g_ prefix, rename fire_rate_gated*."""
    keep = [c for c in df.columns if not c.startswith("d_")]
    df = df[keep].copy()
    df.columns = [rename_col(c) for c in df.columns]
    return df


def already_migrated(df: pd.DataFrame) -> bool:
    return not any(c.startswith(("g_", "d_")) or "fire_rate_gated" in c for c in df.columns)


def migrate_csv(path: Path) -> None:
    df = pd.read_csv(path)
    if already_migrated(df):
        print(f"  skip (already migrated): {path.name}")
        return
    out = migrate_frame(df)
    out.to_csv(path, index=False)
    print(f"  migrated CSV: {path.name}  ({df.shape[1]} -> {out.shape[1]} cols)")


def migrate_parquets(pattern: str) -> int:
    files = sorted(glob.glob(pattern))
    n = 0
    for f in files:
        df = pd.read_parquet(f)
        if already_migrated(df):
            continue
        migrate_frame(df).to_parquet(f, index=False)
        n += 1
    print(f"  migrated {n}/{len(files)} parquets matching {Path(pattern).name}")
    return n


def merge_residualized_into_summary() -> None:
    summ_path = PG / "pooled_metrics/summary.csv"
    resid_path = PG / "position_residualized_summary.csv"
    summ = pd.read_csv(summ_path)
    resid = pd.read_csv(resid_path)
    if any(c in summ.columns for c in RESID_COLS):
        print("  skip resid-merge (already present in summary.csv)")
        return
    keep = ["DMS_id", "concept", "feature"] + [c for c in RESID_COLS if c in resid.columns]
    merged = summ.merge(resid[keep], on=["DMS_id", "concept", "feature"], how="left")
    merged.to_csv(summ_path, index=False)
    print(f"  merged residualized columns into summary.csv "
          f"({len([c for c in keep if c not in ('DMS_id','concept','feature')])} cols added)")


def main() -> None:
    print("Migrating CSVs:")
    migrate_csv(PG / "position_residualized_summary.csv")   # rename first so RESID_COLS exist
    migrate_csv(PG / "pooled_metrics/summary.csv")
    migrate_csv(PG / "full_feature_spearman/full_feature_spearman.csv")
    print("Migrating parquets:")
    migrate_parquets(str(PG / "pooled_metrics/*__pooled.parquet"))
    migrate_parquets(str(PG / "full_feature_spearman/*__full_feat.parquet"))
    print("Merging residualized columns into summary.csv:")
    merge_residualized_into_summary()
    print("done.")


if __name__ == "__main__":
    main()
