"""ProteinGym paper figures, regenerated on the full 135-assay set.

Four figures, each a single legible comparison:
  1. headline      — distribution of our best-feature-per-assay |rho| over 135 assays
  2. vs_ceiling    — our median |rho| vs the supervised/SOTA ceiling (VESPA/VespaG/SaProt)
  3. f1_vs_fitness — InterPLM pairing strength predicts ProteinGym signal (scatter)
  4. residualized  — raw vs position-residualized |rho| (signal is substitution-specific)
"""

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
CEIL = "#888888"


def fig_headline(comp):
    ours = comp["crosscoder_d_pool_mean"].abs()
    med = ours.median()
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.hist(ours, bins=20, range=(0, 1), color=OURS, alpha=0.85)
    ax.axvline(med, color="black", lw=2, ls="--")
    ax.text(med + 0.02, ax.get_ylim()[1] * 0.9, f"median {med:.2f}", fontsize=11)
    ax.set_xlabel("Spearman |ρ|  (best feature per assay)")
    ax.set_ylabel("Number of DMS assays")
    ax.set_title(f"Feature activation vs. DMS fitness across {len(ours)} ProteinGym assays")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(OUTDIR / "pg_01_headline.png", dpi=200); plt.close(fig)
    return med


def fig_vs_ceiling(comp):
    models = [("Ours\n(unsupervised)", "crosscoder_d_pool_mean", OURS),
              ("VESPA", "VESPA", CEIL),
              ("VespaG", "VespaG", CEIL),
              ("SaProt 650M\n(SOTA)", "SaProt (650M)", CEIL)]
    meds = [comp[col].abs().median() for _, col, _ in models]
    labels = [m[0] for m in models]
    cols = [m[2] for m in models]
    fig, ax = plt.subplots(figsize=(7, 4.2))
    bars = ax.bar(labels, meds, color=cols, width=0.6)
    for b, m in zip(bars, meds):
        ax.text(b.get_x() + b.get_width() / 2, m + 0.01, f"{m:.2f}", ha="center", fontsize=11)
    ax.set_ylabel("Median Spearman |ρ|  (135 assays)")
    ax.set_ylim(0, max(meds) * 1.18)
    ax.set_title("Unsupervised crosscoder features vs. the supervised / SOTA ceiling")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(OUTDIR / "pg_02_vs_ceiling.png", dpi=200); plt.close(fig)
    return dict(zip(labels, meds))


def fig_f1_vs_fitness(summ):
    d = summ.dropna(subset=["f1_per_domain", "d_pool_mean_abs"]).copy()
    d["y"] = d["d_pool_mean_abs"].abs()
    rho, p = spearmanr(d["f1_per_domain"], d["y"])
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.scatter(d["f1_per_domain"], d["y"], s=10, alpha=0.25, color=OURS, edgecolors="none")
    # binned mean trend for legibility
    bins = np.linspace(0, 1, 11)
    idx = np.digitize(d["f1_per_domain"], bins) - 1
    xs, ys = [], []
    for b in range(len(bins) - 1):
        m = idx == b
        if m.sum() >= 5:
            xs.append((bins[b] + bins[b + 1]) / 2); ys.append(d["y"][m].median())
    ax.plot(xs, ys, color="black", lw=2.2, marker="o", ms=4, label="binned median")
    ax.set_xlabel("InterPLM pairing strength  (per-domain F1)")
    ax.set_ylabel("ProteinGym signal  |ρ|")
    ax.set_title(f"Stronger annotations → stronger fitness signal  (ρ = {rho:+.2f}, p = {p:.0e})")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(OUTDIR / "pg_03_f1_vs_fitness.png", dpi=200); plt.close(fig)
    return rho, p


def fig_residualized(resid):
    # best feature per (assay, concept) by *residualized* |d_pool_mean_abs|
    # (the most substitution-specific feature; matches the hand-in's 0.289 headline)
    d = resid.dropna(subset=["d_pool_mean_abs_raw", "d_pool_mean_abs_resid"]).copy()
    d["raw"] = d["d_pool_mean_abs_raw"].abs()
    d["res"] = d["d_pool_mean_abs_resid"].abs()
    best = d.loc[d.groupby(["DMS_id", "concept"])["res"].idxmax()]
    mres = best["res"].median()
    frac03 = (best["res"] > 0.3).mean()
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    ax.scatter(best["raw"], best["res"], s=14, alpha=0.4, color=OURS, edgecolors="none")
    ax.plot([0, 1], [0, 1], color="black", lw=1, ls="--")
    ax.axhline(0.3, color="gray", lw=1, ls=":")
    ax.set_xlabel("Raw  |ρ|")
    ax.set_ylabel("Position-residualized  |ρ|")
    ax.set_title(f"Signal is substitution-specific\nmedian residualized |ρ| = {mres:.2f}  "
                 f"({frac03:.0%} of pairs > 0.3)")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(OUTDIR / "pg_04_residualized.png", dpi=200); plt.close(fig)
    return best["raw"].median(), mres


def main():
    comp = pd.read_csv(PG / "vespa_comparison/comparison.csv")
    summ = pd.read_csv(PG / "pooled_metrics/summary.csv")
    resid = pd.read_csv(PG / "position_residualized_summary.csv")

    med = fig_headline(comp)
    ceil = fig_vs_ceiling(comp)
    rho, p = fig_f1_vs_fitness(summ)
    mr, mres = fig_residualized(resid)

    print("wrote pg_01..pg_04 to", OUTDIR)
    print(f"headline median |rho| (best feat/assay) = {med:.3f}  over {len(comp)} assays")
    print("medians:", {k: round(v, 3) for k, v in ceil.items()})
    print(f"f1_per_domain vs signal: rho={rho:+.3f} p={p:.1e}")
    print(f"residualized: raw {mr:.3f} -> resid {mres:.3f}")


if __name__ == "__main__":
    main()
