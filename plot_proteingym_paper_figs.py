"""ProteinGym paper figures — label-validity framing (no leaderboard, no pseudo-LL).

ProteinGym here is not a fitness benchmark we compete on; it is ground-truth fitness data to test
whether InterPLM-annotated crosscoder features track the effect of mutations. The readout is the
crosscoder feature activation (post-BatchTopK). Comparisons to dedicated fitness predictors
(pseudo-LL / SaProt / VESPA) were removed because their method differs from ours.

Three figures (all from CSVs already on disk; no GPU):
  1. f1_vs_fitness   — InterPLM pairing strength predicts the per-variant fitness signal (+0.27)
  2. chance_control  — concept-matched / best-paired vs best-of-N random vs unpaired vs dead floor,
                       plus (when present) the random-init-MODEL null (baseline_model_best)
  3. holdout         — held-out-split test |rho| vs in-sample selected (InterPLM/InterProt guard)
(The per-residue / MotifAE comparison was dropped from the paper — it is confounded and does not
separate trained from random-init; see documentation/experiments/03-proteingym.md.)

Usage:
  # pre-auxfix defaults (unchanged):
  uv run --with pandas --with numpy --with scipy --with matplotlib python \
    repos/sparse-crosscoders-prott5/plot_proteingym_paper_figs.py
  # auxfix chance-control figure with the random-init baseline-model null (only cc CSV needed):
  uv run ... plot_proteingym_paper_figs.py \
    --cc data/proteingym/chance_control_summary_auxfix.csv --suffix _auxfix
fig_01 / fig_03 are skipped automatically if their (pooled / holdout) CSVs are absent — so the auxfix
run, which currently only has the chance-control CSV, renders just the chance-control panel.
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path("/Users/sohrab.tawana/private/crosscoder")
PG = ROOT / "data/proteingym"
OUTDIR = ROOT / "data/figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

OURS = "#1f77b4"
NULL = "#b0b0b0"


def fig_f1_vs_fitness(summ, suffix=""):
    d = summ.dropna(subset=["f1_per_domain", "pool_mean_abs"]).copy()
    d["y"] = d["pool_mean_abs"].abs()
    rho, p = spearmanr(d["f1_per_domain"], d["y"])
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.scatter(d["f1_per_domain"], d["y"], s=10, alpha=0.25, color=OURS, edgecolors="none")
    bins = np.linspace(0, 1, 11)
    idx = np.digitize(d["f1_per_domain"], bins) - 1
    xs, ys = [], []
    for b in range(len(bins) - 1):
        m = idx == b
        if m.sum() >= 5:
            xs.append((bins[b] + bins[b + 1]) / 2); ys.append(d["y"][m].median())
    ax.plot(xs, ys, color="black", lw=2.2, marker="o", ms=4, label="binned median")
    ax.set_xlabel("InterPLM pairing strength  (per-domain F1)")
    ax.set_ylabel("ProteinGym signal  |ρ|  (pool_mean_abs)")
    ax.set_title(f"Stronger annotations → stronger fitness signal  (ρ = {rho:+.2f}, p = {p:.0e})")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(OUTDIR / f"pg_01_f1_vs_fitness{suffix}.png", dpi=200); plt.close(fig)
    return rho, p


def fig_chance_control(cc, summ=None, suffix=""):
    # Bars built from the chance-control summary; the concept-matched bar needs the pooled summary
    # and is included only when that CSV is available. baseline_model_best (the crosscoder trained
    # on random-init ProtT5) is included only when present in the CSV.
    bars = [("dead / typical\nfeature (floor)", float(cc["dead_best"].median()), NULL)]
    if summ is not None:
        concept_matched = (summ.assign(a=summ["pool_mean_abs"].abs())
                           .sort_values("a", ascending=False).drop_duplicates("DMS_id")["a"].median())
        bars.append(("concept-matched\npaired", float(concept_matched), OURS))
    bars.append(("best paired\n(InterPLM)", float(cc["paired_best"].median()), OURS))
    bars.append(("best of K\nRANDOM live", float(cc["random_live_best_K"].median()), NULL))
    if "baseline_model_best" in cc.columns:
        bars.append(("best of\nRANDOM-INIT model", float(cc["baseline_model_best"].median()), NULL))
    bars.append(("best unpaired\nlive", float(cc["unpaired_live_best"].median()), NULL))

    labels = [b[0] for b in bars]; vals = [b[1] for b in bars]; cols = [b[2] for b in bars]
    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    rects = ax.bar(labels, vals, color=cols, width=0.7, edgecolor="black", lw=0.3)
    for r, v in zip(rects, vals):
        ax.text(r.get_x() + r.get_width() / 2, v + 0.008, f"{v:.2f}", ha="center", fontsize=10)
    ax.set_ylabel("Median best-feature-per-assay |ρ|")
    ax.set_ylim(0, max(vals) * 1.2)
    ax.set_title("Chance control: random feature subsets — and a random-init crosscoder — beat the\n"
                 "annotated ones (best-of-N selection inflates |ρ|; InterPLM labels do not enrich for fitness)",
                 loc="left", pad=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(OUTDIR / f"pg_02_chance_control{suffix}.png", dpi=200); plt.close(fig)


def fig_holdout(ho, suffix=""):
    tr = ho["holdout_train_sel"].dropna().median()
    te = ho["holdout_test_sel"].dropna().median()
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    rects = ax.bar(["selected on TRAIN\n(in-sample)", "same feature on\nHELD-OUT TEST"],
                   [tr, te], color=[NULL, OURS], width=0.6, edgecolor="black", lw=0.3)
    for r, v in zip(rects, [tr, te]):
        ax.text(r.get_x() + r.get_width() / 2, v + 0.006, f"{v:.2f}", ha="center", fontsize=11)
    ax.set_ylabel("Median |ρ|  (pool_mean_abs)")
    ax.set_ylim(0, max(tr, te) * 1.25)
    ax.set_title("Position-grouped held-out split (InterPLM/InterProt guard)\n"
                 "selection inflation drops the signal to the honest held-out value", loc="left", pad=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(OUTDIR / f"pg_03_holdout{suffix}.png", dpi=200); plt.close(fig)
    return tr, te


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cc", type=Path, default=PG / "chance_control_summary.csv")
    ap.add_argument("--summ", type=Path, default=PG / "pooled_metrics/summary.csv",
                    help="pooled summary; fig_01 + the concept-matched bar are skipped if absent")
    ap.add_argument("--ho", type=Path, default=PG / "per_residue_holdout_summary.csv",
                    help="held-out summary; fig_03 is skipped if absent")
    ap.add_argument("--suffix", default="", help="appended to output filenames, e.g. _auxfix")
    args = ap.parse_args()

    cc = pd.read_csv(args.cc)
    summ = pd.read_csv(args.summ) if args.summ.exists() else None

    if summ is not None:
        rho, p = fig_f1_vs_fitness(summ, args.suffix)
        print(f"fig_01 f1_per_domain vs signal: rho={rho:+.3f} p={p:.1e}")
    else:
        print(f"(skip fig_01 / concept-matched bar — no pooled summary at {args.summ})")

    fig_chance_control(cc, summ, args.suffix)
    print(f"fig_02 chance control: " + ", ".join(
        f"{c}={cc[c].median():.3f}" for c in
        ["paired_best", "random_live_best_K", "baseline_model_best", "unpaired_live_best", "dead_best"]
        if c in cc.columns))

    if args.ho.exists():
        tr, te = fig_holdout(pd.read_csv(args.ho), args.suffix)
        print(f"fig_03 holdout: train-selected {tr:.3f} -> held-out test {te:.3f}")
    else:
        print(f"(skip fig_03 — no holdout summary at {args.ho})")

    print("wrote figures to", OUTDIR)


if __name__ == "__main__":
    main()
