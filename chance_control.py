"""Chance-level controls for the ProteinGym feature-fitness correlations.

A per-assay best-feature |Spearman| is only interpretable against a null: maxing |ρ| over many
directions inflates it, so the question is whether the InterPLM-paired (or any live) feature beats
what you'd get from arbitrary directions. Per assay we compute |ρ| (pool_mean_abs vs DMS_score) for:

  paired_best          best of the InterPLM-paired features
  unpaired_live_best   best of the unpaired LIVE features
  random_live_best_K   best of K random LIVE features, K = #paired-live (the label-vs-random control,
                       averaged over R draws) -> "does the label pick better than a random subset?"
  dead_best            best of the training-DEAD features (max activation 0 on the eval set)
                       -> null-direction floor
  live_median          median over all live features -> a typical feature

then reports the median across the 135 assays.

Reads:
  data/proteingym/full_feature_spearman/full_feature_spearman.csv
  <ckpt>/uniprotkb_modern_score45_67k/max_activations_per_feature.pt   (dead mask: max==0)
  data/crosscoder_eval/pre-auxfix/real/.../heldout_all_top_pairings.csv (paired feature ids)
Writes:
  data/proteingym/chance_control_summary.csv

Usage (from the crosscoder root):
  uv run --with pandas --with numpy --with torch python repos/sparse-crosscoders-prott5/chance_control.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path("/Users/sohrab.tawana/private/crosscoder")
PG = ROOT / "data/proteingym"
R_DRAWS = 200
SEED = 0

# default (pre-auxfix / bugged) inputs — overridable via CLI so the same analysis runs on auxfix
DEF_FF = PG / "full_feature_spearman/full_feature_spearman.csv"
DEF_MAXACT = (ROOT / "model_checkpoints/crosscoder_l8192_k32_bs512_full_2026-03-12_06-03-41/"
              "crashed_epoch_0_step_2519836/uniprotkb_modern_score45_67k/max_activations_per_feature.pt")
DEF_PAIRINGS = (ROOT / "data/crosscoder_eval/pre-auxfix/real/uniprotkb_modern_score45_67k/"
                "test_counts/heldout_all_top_pairings.csv")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ff", type=Path, default=DEF_FF, help="full_feature_spearman.csv")
    ap.add_argument("--maxact", type=Path, default=DEF_MAXACT, help="max_activations_per_feature.pt (dead mask)")
    ap.add_argument("--pairings", type=Path, default=DEF_PAIRINGS, help="heldout_all_top_pairings.csv")
    ap.add_argument("--baseline_ff", type=Path, default=PG / "full_feature_spearman_baseline/full_feature_spearman.csv")
    ap.add_argument("--metric", default="pool_mean_abs_sp", choices=["pool_mean_abs_sp", "delta_at_pos_sp"])
    ap.add_argument("--out", type=Path, default=PG / "chance_control_summary.csv")
    args = ap.parse_args()

    ff = pd.read_csv(args.ff)
    ff["rho"] = ff[args.metric].abs().fillna(0.0)

    maxact = torch.load(args.maxact, map_location="cpu").numpy()
    dead = set(np.where(maxact == 0)[0].tolist())
    paired = set(pd.read_csv(args.pairings)["feature"].unique())
    all_feats = set(range(len(maxact)))
    live = all_feats - dead
    paired_live = sorted(paired & live)
    unpaired_live = sorted(live - paired)
    K = len(paired_live)
    rng = np.random.default_rng(SEED)
    live_arr = np.array(sorted(live))
    print(f"features: {len(all_feats)} total | {len(live)} live | {len(dead)} dead | "
          f"{len(paired)} paired ({K} paired-live) | metric={args.metric}")

    rows = []
    for dms, g in ff.groupby("DMS_id"):
        rho = g.set_index("feature")["rho"]
        def best(feats):
            s = rho.reindex(feats).dropna()
            return float(s.max()) if len(s) else float("nan")
        # label-vs-random: best of K random live features, averaged over draws
        rand_best = np.mean([rho.reindex(rng.choice(live_arr, size=K, replace=False)).max()
                             for _ in range(R_DRAWS)])
        rows.append({
            "DMS_id": dms,
            "paired_best": best(paired_live),
            "unpaired_live_best": best(unpaired_live),
            "random_live_best_K": float(rand_best),
            "dead_best": best(sorted(dead)),
            "live_median": float(rho.reindex(sorted(live)).median()),
        })

    df = pd.DataFrame(rows)

    # Random-init-model null (SAEBench/InterPLM convention): best feature per assay from the
    # crosscoder trained on a randomly-initialized ProtT5. Auto-included once the cluster job
    # (submit_baseline_proteingym.sh) has written this CSV.
    cols = ["paired_best", "unpaired_live_best", "random_live_best_K", "dead_best", "live_median"]
    if args.baseline_ff.exists():
        bff = pd.read_csv(args.baseline_ff)
        bff["rho"] = bff[args.metric].abs().fillna(0.0)
        bbest = bff.groupby("DMS_id")["rho"].max().rename("baseline_model_best")
        df = df.merge(bbest, on="DMS_id", how="left")
        cols.append("baseline_model_best")
    else:
        print(f"(baseline-model null not present yet -> {args.baseline_ff})\n")

    df.to_csv(args.out, index=False)
    print(f"wrote {args.out} ({len(df)} assays; K={K} paired-live, R={R_DRAWS} draws)\n")
    print(f"{'quantity':<22s} {'median |rho|':>12s}   {'>0.3':>5s} {'>0.5':>5s}")
    for c in cols:
        v = df[c].dropna()
        print(f"{c:<22s} {v.median():>12.3f}   {int((v>0.3).sum()):>5d} {int((v>0.5).sum()):>5d}")


if __name__ == "__main__":
    main()
