"""Build the EVAL-INDEPENDENT assay list (--matches) for harvest_full_feature_spearman.py.

The full-feature ProteinGym scan (§5.6 supplementary: "do ANY of the 8192 features
track fitness across the benchmark?") must NOT depend on the InterPLM pairings —
otherwise the assay set would be conditioned on which features got annotated. So we
build the assay list straight from the ProteinGym reference table, covering ALL 217
substitution assays.

Inputs
  --reference   data/external/DMS_substitutions.csv  (ProteinGym reference; DMS_id + seq_len)
Outputs
  --out         CSV with columns: DMS_id, target_seq_len   (one row per assay)

harvest_full_feature_spearman.py reads this via --matches and applies its own
--max_seq_len filter (default 2048, which drops the 4 assays > 2048 aa). We emit
all 217 here so the drop is visible and controllable downstream, not baked in.

Disposition: COMMITTED reproducibility script (it produces a reported artifact —
the assay list backing the §5.6 scan). Deterministic; no randomness.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", type=Path,
                    default=Path("data/external/DMS_substitutions.csv"))
    ap.add_argument("--out", type=Path,
                    default=Path("data/proteingym_full_feature_matches.csv"))
    ap.add_argument("--max_seq_len", type=int, default=None,
                    help="optional: also print how many assays a downstream cap drops")
    args = ap.parse_args()

    ref = pd.read_csv(args.reference)
    missing = [c for c in ("DMS_id", "seq_len") if c not in ref.columns]
    if missing:
        raise SystemExit(f"reference missing columns {missing}")

    matches = (
        ref[["DMS_id", "seq_len"]]
        .rename(columns={"seq_len": "target_seq_len"})
        .drop_duplicates("DMS_id")
        .sort_values("DMS_id")
        .reset_index(drop=True)
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    matches.to_csv(args.out, index=False)

    print(f"Wrote {len(matches)} assays -> {args.out}")
    print(f"  target_seq_len: min={matches.target_seq_len.min()} "
          f"median={int(matches.target_seq_len.median())} "
          f"max={matches.target_seq_len.max()}")
    cap = args.max_seq_len or 2048
    n_over = int((matches.target_seq_len > cap).sum())
    print(f"  assays > {cap} aa (dropped by harvest --max_seq_len {cap}): {n_over} "
          f"-> {len(matches) - n_over} scanned")


if __name__ == "__main__":
    main()
