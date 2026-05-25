"""
Generate validation + exploration plots from ProteinGym scoring output.

Reads:
  --summary       <output_dir>/summary.csv
  --matches       data/proteingym_concept_matches.csv  (for selection_type joins)
  --scoring_dir   <output_dir>/                         (per-variant parquets,
                                                         only used for --deep_dive)

Writes a set of PNGs to --plots_dir.

Usage:
    uv run --with pandas --with scipy --with matplotlib --with seaborn --with pyarrow \\
      python repos/sparse-crosscoders-prott5/analyze_proteingym_scoring.py \\
        --summary data/proteingym/scoring/summary.csv \\
        --matches data/proteingym_concept_matches.csv \\
        --scoring_dir data/proteingym/scoring \\
        --plots_dir data/proteingym/plots

Specific plots can also be triggered individually — see --only.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr


# ---------------------------------------------------------------------------
# data loading + joining
# ---------------------------------------------------------------------------

def load(summary_csv: Path, matches_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(summary_csv)
    matches = pd.read_csv(matches_csv)[[
        "DMS_id", "matched_concept", "coarse_selection_type",
        "selection_type", "frac_variants_in_concept", "target_seq_len",
    ]]
    df = df.merge(
        matches, left_on=["DMS_id", "concept"],
        right_on=["DMS_id", "matched_concept"], how="left",
    ).drop(columns=["matched_concept"])
    # Convenience columns
    df["abs_sp"] = df["spearman_abs"].abs()
    df["abs_sgn"] = df["spearman_signed"].abs()
    df["is_live"] = df["mean_abs_delta"] > 0
    return df


# ---------------------------------------------------------------------------
# Headline validation plot — Plot 1
# ---------------------------------------------------------------------------

def plot_headline(df: pd.DataFrame, out: Path):
    """Scatter: f1_per_domain (x) vs |spearman_abs| (y).
    Color by coarse_selection_type. Annotated with overall Spearman."""
    live = df[df["is_live"]]
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(
        data=live, x="f1_per_domain", y="abs_sp",
        hue="coarse_selection_type", alpha=0.5, s=20, ax=ax,
    )
    # Linear trendline
    x = live["f1_per_domain"].values
    y = live["abs_sp"].values
    coef = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, coef[0] * xs + coef[1], "k--", lw=1.2, alpha=0.6, label="fit")
    rho, p = spearmanr(x, y)
    ax.set_title(
        f"InterPLM pairing strength vs ProteinGym fitness signal\n"
        f"Spearman = {rho:+.3f},  p = {p:.1e},  n = {len(live):,} live (feature, assay, concept)"
    )
    ax.set_xlabel("InterPLM f1_per_domain (concept-feature pairing strength)")
    ax.set_ylabel("|Spearman| of |Δactivation| vs |ΔDMS|")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Validation broken down by selection class — Plot 2
# ---------------------------------------------------------------------------

def plot_by_selection_class(df: pd.DataFrame, out: Path):
    """Box/strip plots: |spearman_abs| distribution per coarse_selection_type,
    plus the within-class meta-Spearman as text."""
    live = df[df["is_live"]].copy()
    classes = live["coarse_selection_type"].dropna().unique()

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    # Left: distribution per class (best feature per (assay,concept))
    best = (live.sort_values("abs_sp", ascending=False)
                .drop_duplicates(["DMS_id", "concept"], keep="first"))
    sns.violinplot(
        data=best, x="coarse_selection_type", y="abs_sp",
        ax=axes[0], cut=0, inner="point",
    )
    axes[0].set_title("Best |Spearman| per (assay, concept), by selection class")
    axes[0].set_xlabel(""); axes[0].set_ylabel("|Spearman_abs| of best feature")
    axes[0].tick_params(axis="x", rotation=20)

    # Right: meta-Spearman per class
    rows = []
    for c in sorted(classes):
        sub = live[live["coarse_selection_type"] == c]
        rho, p = spearmanr(sub["f1_per_domain"], sub["abs_sp"])
        rows.append((c, len(sub), rho, p))
    summ = pd.DataFrame(rows, columns=["class", "n", "rho", "p"])
    summ = summ.sort_values("rho", ascending=True)
    colors = ["#cccccc" if p > 0.05 else "#1f77b4" for p in summ["p"]]
    axes[1].barh(summ["class"], summ["rho"], color=colors)
    for i, (_, r) in enumerate(summ.iterrows()):
        axes[1].text(r["rho"] + 0.005, i,
                     f"ρ={r['rho']:+.2f}, n={r['n']}, p={r['p']:.0e}",
                     va="center", fontsize=9)
    axes[1].set_title("Validity: Spearman(f1_per_domain, |spearman_abs|) per class")
    axes[1].set_xlabel("Meta-Spearman ρ (grey = p > 0.05)")
    axes[1].axvline(0, color="k", lw=0.5)
    axes[1].set_xlim(min(summ["rho"].min(), 0) - 0.05, summ["rho"].max() + 0.20)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Top hits — Plot 3
# ---------------------------------------------------------------------------

def plot_top_hits(df: pd.DataFrame, out: Path, n: int = 25):
    """Horizontal bar chart of the top N (assay, concept, feature) by |spearman_abs|."""
    live = df[df["is_live"]].copy()
    top = live.nlargest(n, "abs_sp").copy()
    top["label"] = top.apply(
        lambda r: f"{r['DMS_id'][:30]} | {r['concept'][:25]} | f/{int(r['feature'])}",
        axis=1,
    )
    fig, ax = plt.subplots(figsize=(10, max(5, n * 0.3)))
    palette = dict(zip(top["coarse_selection_type"].unique(),
                       sns.color_palette("tab10", top["coarse_selection_type"].nunique())))
    colors = [palette.get(c, "grey") for c in top["coarse_selection_type"]]
    ax.barh(top["label"][::-1], top["abs_sp"][::-1], color=colors[::-1])
    for i, (_, r) in enumerate(top[::-1].iterrows()):
        ax.text(r["abs_sp"] + 0.003, i, f"f1pd={r['f1_per_domain']:.2f}",
                va="center", fontsize=8)
    # Legend
    handles = [plt.Rectangle((0,0),1,1, color=c) for c in palette.values()]
    ax.legend(handles, palette.keys(), loc="lower right", fontsize=8)
    ax.set_xlabel("|Spearman_abs|")
    ax.set_title(f"Top {n} (assay, concept, feature) by ProteinGym fitness correlation")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Per-concept distribution — Plot 4
# ---------------------------------------------------------------------------

def plot_per_concept(df: pd.DataFrame, out: Path):
    """Boxplot of best-per-(assay,concept) |spearman_abs| grouped by concept."""
    live = df[df["is_live"]]
    best = (live.sort_values("abs_sp", ascending=False)
                .drop_duplicates(["DMS_id", "concept"], keep="first"))
    order = (best.groupby("concept")["abs_sp"].median()
                  .sort_values(ascending=False).index.tolist())
    n_per = best.groupby("concept").size()
    keep = [c for c in order if n_per[c] >= 2]  # only concepts with ≥2 assays
    sub = best[best["concept"].isin(keep)]
    fig, ax = plt.subplots(figsize=(10, max(5, 0.4 * len(keep))))
    sns.boxplot(data=sub, y="concept", x="abs_sp", order=keep, ax=ax, color="lightblue")
    sns.stripplot(data=sub, y="concept", x="abs_sp", order=keep, ax=ax, color="black", size=3, alpha=0.6)
    for i, c in enumerate(keep):
        ax.text(sub[sub["concept"]==c]["abs_sp"].max() + 0.01, i,
                f"n={n_per[c]}", va="center", fontsize=8)
    ax.set_xlabel("Best |Spearman_abs| per (assay, concept)")
    ax.set_ylabel("")
    ax.set_title("ProteinGym fitness-correlation strength by concept (concepts seen in ≥2 assays)")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Within-pair validity test histogram — Plot 5
# ---------------------------------------------------------------------------

def plot_within_pair(df: pd.DataFrame, out: Path):
    """For (assay, concept) pairs with ≥3 live features, compute the within-pair
    Spearman of f1_per_domain vs |spearman_abs|. Plot as histogram."""
    live = df[df["is_live"]]
    rhos = []
    for (a, c), g in live.groupby(["DMS_id", "concept"]):
        if len(g) < 3:
            continue
        rho, _ = spearmanr(g["f1_per_domain"], g["abs_sp"])
        if not np.isnan(rho):
            rhos.append(rho)
    rhos = np.array(rhos)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(rhos, bins=25, color="#1f77b4", edgecolor="white")
    ax.axvline(0, color="k", lw=1)
    ax.axvline(np.median(rhos), color="red", lw=1.5,
               label=f"median = {np.median(rhos):+.2f}")
    pct_pos = (rhos > 0).mean() * 100
    ax.set_title(
        f"Within-(assay,concept) feature ranking: f1_per_domain vs |spearman_abs|\n"
        f"n = {len(rhos)} pairs with ≥3 live features  |  "
        f"{pct_pos:.0f}% positive  |  median = {np.median(rhos):+.2f}"
    )
    ax.set_xlabel("Within-pair Spearman ρ")
    ax.set_ylabel("Count")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Deep-dive scatter for one (assay, concept, feature) — Plot 6
# ---------------------------------------------------------------------------

def plot_deep_dive(
    scoring_dir: Path, dms_id: str, concept: str, feature: int, out: Path,
):
    """Scatter of Δactivation vs DMS_score, with regression and marginals."""
    safe = concept.replace("/", "_").replace(" ", "_")
    parq = scoring_dir / f"{dms_id}__{safe}__per_variant.parquet"
    if not parq.exists():
        print(f"  WARN: parquet missing: {parq}")
        return
    per_var = pd.read_parquet(parq)
    pv = per_var[per_var["feature"] == feature].copy()
    if pv.empty:
        print(f"  WARN: no rows for feature {feature} in {parq.name}")
        return
    rho, p = spearmanr(pv["delta_act"], pv["DMS_score"])
    rho_a, p_a = spearmanr(pv["delta_act"].abs(),
                           (pv["DMS_score"] - pv["DMS_score"].median()).abs())
    g = sns.jointplot(
        data=pv, x="delta_act", y="DMS_score", kind="scatter",
        height=6, alpha=0.4, s=8, marginal_kws={"bins": 40},
    )
    g.fig.suptitle(
        f"{dms_id}  |  {concept}  |  f/{feature}\n"
        f"signed ρ = {rho:+.3f} (p={p:.1e})    |Δ|↔|ΔDMS| ρ = {rho_a:+.3f} (p={p_a:.1e})    "
        f"n = {len(pv)}",
        y=1.02, fontsize=10,
    )
    g.ax_joint.axhline(pv["DMS_score"].median(), color="grey", lw=0.5, ls="--")
    g.ax_joint.axvline(0, color="grey", lw=0.5, ls="--")
    g.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(g.fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--matches", type=Path, required=True)
    ap.add_argument("--scoring_dir", type=Path, default=None,
                    help="Per-variant parquet dir (needed for deep-dive)")
    ap.add_argument("--plots_dir", type=Path, required=True)
    ap.add_argument("--only", type=str, default=None,
                    help="Comma-list of plot names to run "
                         "(headline,by_class,top_hits,per_concept,within_pair,deep_dive). "
                         "Default: all.")
    ap.add_argument("--deep_dive_top_n", type=int, default=6,
                    help="Number of top hits to deep-dive into")
    args = ap.parse_args()

    args.plots_dir.mkdir(parents=True, exist_ok=True)
    df = load(args.summary, args.matches)
    print(f"Loaded {len(df):,} summary rows; {df['is_live'].sum():,} live")

    plots_to_run = set(args.only.split(",")) if args.only else {
        "headline", "by_class", "top_hits", "per_concept", "within_pair", "deep_dive",
    }

    if "headline" in plots_to_run:
        plot_headline(df, args.plots_dir / "01_headline.png"); print("  wrote 01_headline.png")
    if "by_class" in plots_to_run:
        plot_by_selection_class(df, args.plots_dir / "02_by_class.png"); print("  wrote 02_by_class.png")
    if "top_hits" in plots_to_run:
        plot_top_hits(df, args.plots_dir / "03_top_hits.png"); print("  wrote 03_top_hits.png")
    if "per_concept" in plots_to_run:
        plot_per_concept(df, args.plots_dir / "04_per_concept.png"); print("  wrote 04_per_concept.png")
    if "within_pair" in plots_to_run:
        plot_within_pair(df, args.plots_dir / "05_within_pair.png"); print("  wrote 05_within_pair.png")
    if "deep_dive" in plots_to_run:
        if not args.scoring_dir:
            print("  SKIP deep_dive: pass --scoring_dir to enable")
        else:
            top = (df[df["is_live"]]
                   .nlargest(args.deep_dive_top_n, "abs_sp")
                   .reset_index(drop=True))
            for _, r in top.iterrows():
                fname = (
                    f"06_deepdive__{r['DMS_id']}__"
                    f"{r['concept'].replace('/','_').replace(' ','_')}__f{int(r['feature'])}.png"
                )
                plot_deep_dive(args.scoring_dir, r["DMS_id"], r["concept"],
                               int(r["feature"]), args.plots_dir / fname)
                print(f"  wrote {fname}")

    print(f"\nDone — plots in {args.plots_dir}")


if __name__ == "__main__":
    main()
