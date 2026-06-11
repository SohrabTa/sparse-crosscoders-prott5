"""Per-layer decoder-norm profiles of concept-paired crosscoder features,
colored by the biological concept *category* the feature is associated with.

Replicates the original sparse-crosscoder per-layer-norm figure, but recolors
the lines by biology (Domain / Region / Disulfide bond / ...) instead of peak layer.

x-axis: ProtT5 encoder layer (1..24)
y-axis: feature decoder norm at each layer, per-feature max-normalized to [0,1]
color : concept category of the feature's best (per-domain-F1) pairing
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import torch

ROOT = Path("/Users/sohrab.tawana/private/crosscoder")
CKPT = ROOT / "model_checkpoints/crosscoder_l8192_k32_bs512_full_2026-03-12_06-03-41/crashed_epoch_0_step_2519836/model.pt"
PAIRINGS = ROOT / "data/crosscoder_eval/pre-auxfix/real/uniprotkb_modern_score45_67k/test_counts/heldout_all_top_pairings.csv"
OUTDIR = ROOT / "data/figures"
OUT = OUTDIR / "layerwise_feature_anatomy.png"
OUT_FACET = OUTDIR / "layerwise_feature_anatomy_faceted.png"
OUTDIR.mkdir(parents=True, exist_ok=True)


def decoder_norms_per_layer(state, unfold_scale=True):
    """[n_latents, n_layers] scale-invariant decoder norms."""
    w_dec = state["_W_dec_LXoDo"][:, 0, :, :]          # [L=8192, P=24, D=1024]
    norms = w_dec.norm(dim=-1)                          # [8192, 24]
    if unfold_scale and bool(state["is_folded"]):
        scale_P = state["folded_scaling_factors_out_Xo"][0]   # [24]
        norms = norms * scale_P[None, :]
    return norms


def best_pairing_per_feature(df):
    """One category per feature: the pairing with the highest per-domain F1."""
    idx = df.groupby("feature")["f1_per_domain"].idxmax()
    best = df.loc[idx, ["feature", "concept"]].copy()
    best["category"] = best["concept"].str.split("_").str[0]
    return best


def main():
    state = torch.load(CKPT, map_location="cpu", weights_only=False)
    norms = decoder_norms_per_layer(state)             # [8192, 24]
    norms = norms / norms.max(dim=1, keepdim=True).values.clamp_min(1e-12)
    norms = norms.numpy()
    n_layers = norms.shape[1]

    pairings = pd.read_csv(PAIRINGS)
    best = best_pairing_per_feature(pairings)

    counts = best["category"].value_counts()
    order = counts.index.tolist()                       # largest first
    cmap = plt.get_cmap("tab10")
    colors = {cat: cmap(i % 10) for i, cat in enumerate(order)}
    x = list(range(1, n_layers + 1))

    def feats(cat):
        return [int(f) for f in best.loc[best["category"] == cat, "feature"]]

    # ---- Figure 1: combined panel — faint individual lines + bold per-category mean ----
    fig, ax = plt.subplots(figsize=(9, 5.2))
    for cat in order:
        ids = feats(cat)
        for fid in ids:
            ax.plot(x, norms[fid], color=colors[cat], alpha=0.18, linewidth=0.8)
        if len(ids) >= 3:                               # mean only where it's meaningful
            ax.plot(x, norms[ids].mean(0), color=colors[cat], linewidth=2.6)
    handles = [plt.Line2D([0], [0], color=colors[c], lw=2.6,
                          label=f"{c} (n={counts[c]})") for c in order]
    ax.legend(handles=handles, title="Concept category", fontsize=9, title_fontsize=10,
              loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    ax.set_xlabel("ProtT5 encoder layer")
    ax.set_ylabel("Decoder norm (per-feature normalized)")
    ax.set_xlim(1, n_layers); ax.set_ylim(0, 1.05); ax.set_xticks([1, 6, 12, 18, 24])
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 2: small multiples — one panel per category, shared axes ----
    ncol = 4
    nrow = (len(order) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(13, 3.1 * nrow), sharex=True, sharey=True)
    axes = axes.ravel()
    for i, cat in enumerate(order):
        a = axes[i]
        ids = feats(cat)
        for fid in ids:
            a.plot(x, norms[fid], color=colors[cat], alpha=0.5, linewidth=1.0)
        if len(ids) >= 3:
            a.plot(x, norms[ids].mean(0), color="black", linewidth=2.2)
        a.set_title(f"{cat}  (n={counts[cat]})", fontsize=11)
        a.set_xlim(1, n_layers); a.set_ylim(0, 1.05); a.set_xticks([1, 12, 24])
        a.spines[["top", "right"]].set_visible(False)
    for j in range(len(order), len(axes)):
        axes[j].axis("off")
    fig.supxlabel("ProtT5 encoder layer", fontsize=12)
    fig.supylabel("Decoder norm (per-feature normalized)", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_FACET, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"wrote {OUT}\nwrote {OUT_FACET}")
    print(f"{len(best)} concept-paired features across {len(order)} categories")
    print(counts.to_string())


if __name__ == "__main__":
    main()
