"""Verify a BatchTopK -> global-JumpReLU conversion (diagnostic for experiment 05).

Reuses the estimator in convert_batchtopk_to_jumprelu.py and reports the three
numbers that decide whether the conversion is safe to push through the downstream
InterPLM / ProteinGym evals:

  1. theta convergence  - running mean of the per-batch (k*B)-th largest pre-activation,
                          to confirm the estimate has stabilized (it is an order
                          statistic averaged over batches; Bussmann 2024 Eq. 7 /
                          Minder 2025 Eq. 14 define it as an expectation over batches
                          and prescribe NO token count -- convergence is the criterion).
  2. mean L0            - average nonzero latents per token under the gate; want ~= k.
  3. active-set divergence vs BatchTopK on the SAME tokens:
                          per-token Jaccard of the active sets, and the fraction of
                          (token, latent) activations that flip. `global` should track
                          BatchTopK closely by design (faithful to the train-time rule).

INPUTS
  --crosscoder_dir   trained BatchTopK crosscoder dir (config.yaml + model.pt)
  --checkpoint       checkpoint filename (default model.pt)
  --activations      dir of [N_tokens, M=1, P=24, D=1024] residual-stream .pt tensors
                     (same layout convert_batchtopk_to_jumprelu.py consumes; produced
                     by ProtT5CrosscoderEmbedder.extract_embeddings). OPTIONAL if --fasta.
  --fasta            FASTA/one-seq-per-line file to embed on the fly via
                     ProtT5CrosscoderEmbedder instead of reading --activations.
  --k, --batch_tokens, --max_batches   match training (k=32, batch_tokens=512).

OUTPUTS
  stdout report + optional --out_json with {theta, mean_L0, jaccard, flip_frac,
  n_batches, n_tokens, theta_trace}. Read-only w.r.t. the checkpoint.

Run example (local M1, ~100k tokens):
  python verify_jumprelu_conversion.py \
    --crosscoder_dir <auxfix ckpt> --fasta seqs.txt --k 32 --batch_tokens 512

Disposition: COMMITTED reproducibility script for experiment 05 (the conversion
decision gate). Documented in documentation/experiments/05-auxk-fix-rerun.md.
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
    iter_activation_batches,
    load_batchtopk_crosscoder,
    pre_activation_BL,
)


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
            args.fasta, args.crosscoder_dir / "_verify_acts", device,
            max_tokens=args.batch_tokens * (args.max_batches + 1),
        )

    _, inner = load_batchtopk_crosscoder(args.crosscoder_dir, args.checkpoint, device)
    dtype = next(inner.parameters()).dtype
    k, B = args.k, args.batch_tokens

    def batches():
        for i, b in enumerate(iter_activation_batches(args.activations, B, device, dtype)):
            if i >= args.max_batches:
                break
            yield b

    # Pass 1: theta + running-mean convergence trace.
    per_batch, theta_trace = [], []
    for x in batches():
        p = pre_activation_BL(inner, x).reshape(-1)
        per_batch.append(p.topk(k * B, sorted=True).values[-1].item())
        theta_trace.append(sum(per_batch) / len(per_batch))
    if not per_batch:
        raise RuntimeError("No usable batches; check --activations/--batch_tokens.")
    theta = sum(per_batch) / len(per_batch)
    n_batches = len(per_batch)

    # Pass 2: mean L0 + active-set divergence vs BatchTopK at the final theta.
    jacc, l0_bt, l0_jr, flips, ever = [], [], [], 0, 0
    for x in batches():
        p = pre_activation_BL(inner, x)
        kth = p.reshape(-1).topk(k * B, sorted=True).values[-1]
        bt, jr = p >= kth, p > theta
        inter = (bt & jr).sum(-1).float()
        union = (bt | jr).sum(-1).float().clamp(min=1)
        jacc.append((inter / union).mean().item())
        l0_bt.append(bt.sum(-1).float().mean().item())
        l0_jr.append(jr.sum(-1).float().mean().item())
        flips += int((bt ^ jr).sum())
        ever += int((bt | jr).sum())

    mean = lambda xs: sum(xs) / len(xs)
    rep = {
        "theta": theta, "n_batches": n_batches, "n_tokens": n_batches * B,
        "mean_L0_batchtopk": mean(l0_bt), "mean_L0_jumprelu": mean(l0_jr),
        "active_set_jaccard": mean(jacc), "flip_frac": flips / max(ever, 1),
        "theta_trace": theta_trace, "k": k, "batch_tokens": B,
    }
    print(f"theta={theta:.4f}  batches={n_batches}  tokens={n_batches*B}")
    print(f"theta convergence (running mean): "
          f"b1={theta_trace[0]:.3f}  b{min(10,n_batches)}={theta_trace[min(9,n_batches-1)]:.3f}  "
          f"final={theta_trace[-1]:.3f}")
    print(f"mean L0   BatchTopK={mean(l0_bt):.2f}   JumpReLU={mean(l0_jr):.2f}   (target k={k})")
    print(f"per-token active-set Jaccard (BatchTopK vs JumpReLU): {mean(jacc):.3f}")
    print(f"activations that flip / ever-active: {flips}/{ever} = {flips/max(ever,1):.3f}")
    if args.out_json:
        args.out_json.write_text(json.dumps(rep, indent=2))
        print(f"wrote {args.out_json}")


if __name__ == "__main__":
    main()
