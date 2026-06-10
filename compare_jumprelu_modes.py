"""Ablation: global vs decoder_norm JumpReLU threshold for the BatchTopK->JumpReLU
conversion (experiment 05).

The headline conversion uses a single GLOBAL scalar theta (Bussmann 2024 Eq. 7),
faithful to how this model was trained (BatchTopK selects on the raw pre-activation
f, no decoder-norm scaling). This script asks the robustness question: does the
per-latent decoder_norm variant (Minder et al. 2025, theta_v / s_j) materially
change inference behavior? If the two gates barely differ against BatchTopK, the
global headline is robust.

It estimates BOTH thresholds on the same calibration tokens (reusing
convert_batchtopk_to_jumprelu.estimate_threshold), then on the same held-out
batches compares three gates against the train-time BatchTopK active set:
  - mean L0 (nonzero latents/token; target ~= k)
  - per-token active-set Jaccard vs BatchTopK
  - flip fraction ((token,latent) activations that differ from BatchTopK)
  - FVU (reconstruction error of the gated latents, vs the variance baseline)

Run (local M1, ~100k tokens):
  python compare_jumprelu_modes.py --crosscoder_dir <auxfix BatchTopK ckpt> \
    --checkpoint model.pt --fasta calib_seqs.txt --k 32 --batch_tokens 512 \
    --out_json decoder_norm_ablation.json

NOTE on tokens: like verify_jumprelu_conversion.py, thresholds are estimated and
evaluated on the same token batches. theta (global) and theta_v (decoder_norm) are
batch-aggregate order statistics, not per-token fits, so this is not meaningful
overfitting; the global-vs-decoder_norm comparison is relative and uses identical
tokens for both. Read-only w.r.t. the checkpoint.

Disposition: COMMITTED reproducibility script for the experiment-05 ablation. The
headline conversion stays GLOBAL regardless of outcome; this only reports how much
(if anything) the per-latent variant would change.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import torch

_HERE = Path(__file__).resolve()
_REPO_ROOT = _HERE.parents[2]
sys.path.insert(0, str(_REPO_ROOT / "repos" / "InterPLM"))
sys.path.insert(0, str(_REPO_ROOT / "repos" / "crosscode"))
sys.path.insert(0, str(_HERE.parent))

from convert_batchtopk_to_jumprelu import (  # noqa: E402
    embed_fasta_to_shards,
    estimate_threshold,
    iter_activation_batches,
    load_batchtopk_crosscoder,
    pre_activation_BL,
)


@torch.no_grad()
def gate_stats_vs_batchtopk(inner, batches, k, thr_L_global, thr_L_decnorm):
    """Compare global and decoder_norm JumpReLU gates against BatchTopK on the
    same tokens. Returns per-gate dicts of mean_L0, jaccard, flip_frac, fvu."""
    acc = {
        "batchtopk": {"l0": [], "fvu_num": 0.0, "fvu_den": 0.0},
        "global": {"l0": [], "jacc": [], "flip": 0, "ever": 0, "fvu_num": 0.0, "fvu_den": 0.0},
        "decoder_norm": {"l0": [], "jacc": [], "flip": 0, "ever": 0, "fvu_num": 0.0, "fvu_den": 0.0},
    }

    def fvu_accumulate(d, x_BMPD, latents_BL):
        recon = inner.decode_BMPD(latents_BL)
        num = ((x_BMPD - recon) ** 2).sum().item()
        den = ((x_BMPD - x_BMPD.mean(dim=0, keepdim=True)) ** 2).sum().item()
        d["fvu_num"] += num
        d["fvu_den"] += den

    for x in batches:
        pre = pre_activation_BL(inner, x)  # [B, L]
        B = pre.shape[0]
        # BatchTopK reference: keep the top k*B raw pre-activations across the batch.
        kth = pre.reshape(-1).topk(k * B, sorted=True).values[-1]
        bt_mask = pre >= kth
        bt_lat = pre * bt_mask
        acc["batchtopk"]["l0"].append(bt_mask.sum(-1).float().mean().item())
        fvu_accumulate(acc["batchtopk"], x, bt_lat)

        for name, thr_L in (("global", thr_L_global), ("decoder_norm", thr_L_decnorm)):
            thr = thr_L.to(pre.device)
            g_mask = pre > thr
            g_lat = pre * g_mask
            d = acc[name]
            d["l0"].append(g_mask.sum(-1).float().mean().item())
            inter = (bt_mask & g_mask).sum(-1).float()
            union = (bt_mask | g_mask).sum(-1).float().clamp(min=1)
            d["jacc"].append((inter / union).mean().item())
            d["flip"] += int((bt_mask ^ g_mask).sum())
            d["ever"] += int((bt_mask | g_mask).sum())
            fvu_accumulate(d, x, g_lat)

    mean = lambda xs: sum(xs) / len(xs) if xs else float("nan")
    out = {}
    for name, d in acc.items():
        rec = {"mean_L0": mean(d["l0"]),
               "fvu": d["fvu_num"] / d["fvu_den"] if d["fvu_den"] else float("nan")}
        if name != "batchtopk":
            rec["jaccard_vs_batchtopk"] = mean(d["jacc"])
            rec["flip_frac_vs_batchtopk"] = d["flip"] / max(d["ever"], 1)
        out[name] = rec
    return out


@torch.no_grad()
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--crosscoder_dir", type=Path, required=True)
    ap.add_argument("--checkpoint", type=str, default="model.pt")
    ap.add_argument("--activations", type=Path, default=None)
    ap.add_argument("--fasta", type=Path, default=None)
    ap.add_argument("--k", type=int, default=32)
    ap.add_argument("--batch_tokens", type=int, default=512)
    ap.add_argument("--max_batches", type=int, default=200)
    ap.add_argument("--device", type=str, default=None)
    ap.add_argument("--out_json", type=Path, default=None)
    args = ap.parse_args()

    if args.device:
        device = args.device
    elif torch.cuda.is_available():
        device = "cuda"
    elif torch.backends.mps.is_available():
        device = "mps"
    else:
        device = "cpu"

    if args.activations is None:
        if args.fasta is None:
            ap.error("provide --activations or --fasta")
        args.activations = embed_fasta_to_shards(
            args.fasta, args.crosscoder_dir / "_ablation_acts", device,
            max_tokens=args.batch_tokens * (args.max_batches + 1),
        )

    _, inner = load_batchtopk_crosscoder(args.crosscoder_dir, args.checkpoint, device)
    dtype = next(inner.parameters()).dtype

    def batches():
        for i, b in enumerate(iter_activation_batches(args.activations, args.batch_tokens, device, dtype)):
            if i >= args.max_batches:
                break
            yield b

    log_thr_global, theta_g, nb = estimate_threshold(inner, batches(), args.k, "global")
    log_thr_decnorm, theta_v, _ = estimate_threshold(inner, batches(), args.k, "decoder_norm")
    thr_global = log_thr_global.exp()
    thr_decnorm = log_thr_decnorm.exp()

    stats = gate_stats_vs_batchtopk(inner, batches(), args.k, thr_global, thr_decnorm)

    rep = {
        "theta_global": theta_g, "theta_v_decoder_norm": theta_v,
        "n_batches": nb, "n_tokens": nb * args.batch_tokens,
        "k": args.k, "batch_tokens": args.batch_tokens,
        "decoder_norm_threshold_spread": {
            "min": float(thr_decnorm.min()), "median": float(thr_decnorm.median()),
            "max": float(thr_decnorm.max()),
        },
        "global_threshold": float(thr_global.flatten()[0]),
        "gates": stats,
    }

    print(f"\ntheta(global)={theta_g:.4f}   theta_v(decoder_norm)={theta_v:.4f}   "
          f"batches={nb}  tokens={nb*args.batch_tokens}")
    print(f"global threshold (scalar): {rep['global_threshold']:.4f}")
    print(f"decoder_norm raw-f thresholds: min={thr_decnorm.min():.4g} "
          f"median={thr_decnorm.median():.4g} max={thr_decnorm.max():.4g}")
    print(f"\n{'gate':<14}{'meanL0':>9}{'Jaccard':>10}{'flip_frac':>11}{'FVU':>9}")
    for name in ("batchtopk", "global", "decoder_norm"):
        g = stats[name]
        j = g.get("jaccard_vs_batchtopk", float("nan"))
        fl = g.get("flip_frac_vs_batchtopk", float("nan"))
        print(f"{name:<14}{g['mean_L0']:>9.2f}{j:>10.3f}{fl:>11.3f}{g['fvu']:>9.4f}")
    print("\nIf global and decoder_norm rows are close (Jaccard/flip/FVU), the "
          "global headline conversion is robust to the per-latent variant.")

    if args.out_json:
        args.out_json.write_text(json.dumps(rep, indent=2))
        print(f"wrote {args.out_json}")


if __name__ == "__main__":
    main()
