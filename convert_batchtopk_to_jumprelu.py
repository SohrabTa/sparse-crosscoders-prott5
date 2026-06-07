"""
Convert a BatchTopK-trained crosscoder to a JumpReLU activation for inference.

WHY THIS EXISTS
---------------
We train with BatchTopK (crosscode `BatchTopkActivation`, k_per_example=32): the
top-(k*B) pre-activations across the *whole* token batch survive. That makes a
token's active set depend on the rest of its batch (see insights.md 2026-06-04),
which is wrong for single-sequence inference (InterPLM eval, ProteinGym).

The standard fix (Bussmann, Leask, Nanda 2024, "BatchTopK Sparse Autoencoders",
arXiv:2412.06410) is to replace BatchTopK with a JumpReLU gate at inference:
estimate a threshold theta from data once, then keep latent j iff its
pre-activation clears theta. No retraining; the encoder/decoder/bias are unchanged.

THRESHOLD MODES (this is the part the literature disagrees on — see
references/notes/jumprelu-threshold-conversion.md):

  global       Single scalar theta, the canonical BatchTopK->JumpReLU recipe
               (Bussmann Eq. 7). theta = mean over batches of the minimum
               selected pre-activation = mean over batches of the (k*B)-th
               largest pre-activation. Broadcast to all latents:
                   log_threshold_L[j] = log(theta)   for every j.
               Faithful to how THIS model was trained (BatchTopK selects on the
               raw pre-activation f, with no decoder-norm scaling). DEFAULT.

  decoder_norm Effectively per-latent, the crosscoder recipe of Minder et al.
               2025 ("Overcoming Sparsity Artifacts in Crosscoders"). Scale the
               pre-activation by the summed decoder norm
                   s_j = sum_{m,p} || W_dec[j, m, p, :] ||_2
               estimate theta on v = f * s, then fold the scaling back into a
               per-latent threshold on the raw f:
                   log_threshold_L[j] = log(theta_v / s_j).
               Latents with large decoder norms gate at a lower raw pre-act.
               CAVEAT: training selected on raw f, not on v, so applying the
               scaling only at inference is an approximation, not an exact match
               to the train-time active set. A faithful Minder model would also
               *select* on v during training. Use this to experiment with the
               per-latent variant; prefer `global` for the headline conversion.

Genuine per-feature *learned* thresholds (Rajamanoharan et al. 2024 arXiv:2407.14435;
Ge et al. 2025) are NOT recoverable here: BatchTopK imposes one global cutoff per
batch, so post-hoc you can only estimate that scalar (optionally reshaped by
decoder norm). A true per-feature threshold needs native JumpReLU training
(crosscode jan_update path), which is a different run, not a conversion.

Both modes write into the SAME per-latent parameter that crosscode's JumpReLU
already has (`AnthropicSTEJumpReLUActivation.log_threshold_L`, size = n_latents),
so the existing activation class is reused unmodified — `global` just ties all
entries to one value.

INPUTS
------
  --crosscoder_dir   trained BatchTopK crosscoder dir (config.yaml + model.pt)
  --checkpoint       checkpoint filename inside that dir (e.g. model.pt)
  --activations      directory of saved residual-stream activation tensors, each
                     a torch .pt of shape [N_tokens, M, P, D] matching what the
                     crosscoder consumes (M=1, P=24, D=1024 for ProtT5). This is the
                     layout ProtT5CrosscoderEmbedder.extract_embeddings produces (NOT
                     extract_activations.py, which saves (layers, B, S, D)).
  --fasta            ALTERNATIVELY pass a FASTA / one-seq-per-line file and let this
                     script embed it on the fly (ProtT5CrosscoderEmbedder), so the
                     whole conversion is one self-contained command. Memory-bounded
                     via sharding (see embed_fasta_to_shards). Use --activations OR
                     --fasta, not both.

OUTPUTS
-------
  <out_dir>/jumprelu_threshold.pt   dict: {log_threshold_L: FloatTensor[n_latents],
                                    theta: float, mode: str, k: int,
                                    batch_tokens: int, n_batches: int,
                                    mean_L0_check: float}
  <out_dir>/converted/              (only with --write_converted) a self-contained
                                    JumpReLU checkpoint dir: config.yaml with the
                                    activation_fn block rewritten to
                                    AnthropicSTEJumpReLUActivation and model.pt with
                                    log_threshold_L injected. Loads via the existing
                                    score_proteingym_features.py `keep` path with no
                                    further changes.

The sidecar is also consumed directly by score_proteingym_features.py
--inference_activation jumprelu --jumprelu_threshold <path>.

STATUS: written 2026-06-05 BEFORE the auxfix checkpoint
(crosscoder_l8192_k32_bs512_full_auxfix, experiment 05) landed, so the threshold
math and checkpoint surgery are exercised, but the --activations loader contract
must be confirmed against extract_activations.py output before the first real run.
Documented in documentation/experiments/05-auxk-fix-rerun.md.
"""

from __future__ import annotations

import argparse
import logging
import re
import shutil
import sys
from pathlib import Path
from typing import Iterator

import torch

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "repos" / "InterPLM"))
sys.path.insert(0, str(_REPO_ROOT / "repos" / "crosscode"))

import yaml  # noqa: E402
import crosscode.log  # noqa: E402,F401

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("bt2jr")
log.setLevel(logging.INFO)

JUMPRELU_CLASSNAME = "AnthropicSTEJumpReLUActivation"


def load_batchtopk_crosscoder(crosscoder_dir: Path, checkpoint: str, device: str):
    """Load the trained crosscoder with its native BatchTopK activation."""
    from interplm.sae.inference import load_sae

    wrapper = load_sae(crosscoder_dir, device=device, model_name=checkpoint)
    inner = wrapper.crosscoder
    from crosscode.models.activations.topk import BatchTopkActivation

    if not isinstance(inner.activation_fn, BatchTopkActivation):
        log.warning(
            "Expected BatchTopkActivation, found %s. Estimation still runs on the "
            "raw pre-activations, but double-check this is the model you meant.",
            type(inner.activation_fn).__name__,
        )
    return wrapper, inner


def embed_fasta_to_shards(
    fasta: Path,
    out_dir: Path,
    device: str,
    model_name: str = "Rostlab/prot_t5_xl_uniref50",
    seqs_per_shard: int = 64,
    max_seq_len: int = 512,
    max_tokens: int | None = None,
    store_dtype: torch.dtype = torch.float16,
) -> Path:
    """Embed a FASTA / one-seq-per-line file into a dir of [N_tokens, 1, 24, 1024]
    shards (the layout iter_activation_batches reads), via the SAME embedder the
    evals use (ProtT5CrosscoderEmbedder). Makes the conversion a single command:
    sequences in, threshold out.

    Memory-bounded: embeds `seqs_per_shard` sequences at a time and writes each
    group as its own shard, so peak RAM is one group's activations + the model,
    not the whole calibration set. Stops once `max_tokens` tokens are written.
    Sequences are capped at `max_seq_len` (default 512) to match the train-time
    sequence-length cap (uniref50_length_512) and bound per-shard memory. Shards
    are stored as float16 by default (halves disk; cast back to model dtype on load).
    """
    from interplm.embedders.prott5 import ProtT5CrosscoderEmbedder

    text = fasta.read_text()
    if text.lstrip().startswith(">"):
        parts = re.split(r"^>.*$", text, flags=re.MULTILINE)
        seqs = ["".join(s.split()) for s in parts if s.strip()]
    else:
        seqs = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if max_seq_len:
        seqs = [s[:max_seq_len] for s in seqs]
    if not seqs:
        raise RuntimeError(f"No sequences parsed from {fasta}")

    out_dir.mkdir(parents=True, exist_ok=True)
    log.info("Embedding up to %d seqs (cap %d aa) via %s -> %s", len(seqs), max_seq_len, model_name, out_dir)
    embedder = ProtT5CrosscoderEmbedder(model_name=model_name, device=device)

    written, shard_i = 0, 0
    for start in range(0, len(seqs), seqs_per_shard):
        group = seqs[start : start + seqs_per_shard]
        acts = embedder.extract_embeddings(group, batch_size=4)  # (N, 1, 24, 1024)
        torch.save(acts.to(store_dtype).contiguous(), out_dir / f"shard_{shard_i:04d}.pt")
        written += int(acts.shape[0])
        shard_i += 1
        log.info("  shard %d: +%d tokens (total %d)", shard_i, acts.shape[0], written)
        if max_tokens is not None and written >= max_tokens:
            break
    if written == 0:
        raise RuntimeError(f"No tokens embedded from {fasta}")
    log.info("Embedded %d tokens into %d shards under %s", written, shard_i, out_dir)
    return out_dir


def iter_activation_batches(
    activations_dir: Path, batch_tokens: int, device: str, dtype: torch.dtype
) -> Iterator[torch.Tensor]:
    """Yield [batch_tokens, M, P, D] activation chunks from saved .pt tensors.

    Contract: each file under activations_dir is a torch tensor of shape
    [N_tokens, M, P, D] (the residual-stream activations the crosscoder consumes,
    same layout as crosscoder_dictionary.encode's input). We concatenate across
    files and re-chunk to batch_tokens so the BatchTopK cutoff statistics are
    estimated over token-batches of the SAME size used at train time (512).

    INTEGRATION POINT: if extract_activations.py saves a different layout (e.g.
    per-sequence, or [M, P, N, D]), adapt the load/reshape here. To stream from
    the ProtT5 embedder instead of a store, replace this generator with one that
    runs ProtT5CrosscoderEmbedder over a FASTA and yields the same [B, M, P, D].
    """
    files = sorted(p for p in activations_dir.iterdir() if p.suffix == ".pt")
    if not files:
        raise FileNotFoundError(f"No .pt activation tensors under {activations_dir}")

    buf: list[torch.Tensor] = []
    n_buffered = 0
    for f in files:
        t = torch.load(f, map_location="cpu", weights_only=True)
        if t.ndim != 4:
            raise ValueError(
                f"{f}: expected [N, M, P, D], got shape {tuple(t.shape)}. "
                "Adapt iter_activation_batches to the real layout."
            )
        buf.append(t)
        n_buffered += t.shape[0]
        while n_buffered >= batch_tokens:
            cat = torch.cat(buf, dim=0)
            yield cat[:batch_tokens].to(device=device, dtype=dtype)
            rest = cat[batch_tokens:]
            buf = [rest] if rest.shape[0] else []
            n_buffered = rest.shape[0]
    # trailing partial batch is dropped: BatchTopK's cutoff depends on the batch
    # size, so a short final batch would bias theta downward.


@torch.no_grad()
def pre_activation_BL(inner, x_BMPD: torch.Tensor) -> torch.Tensor:
    """Latent pre-activation f = (encoder x) + b_enc, the value BatchTopK selects on."""
    pre = inner.get_pre_bias_BL(x_BMPD)
    if inner.b_enc_L is not None:
        pre = pre + inner.b_enc_L
    return pre


def decoder_norm_scale_L(inner) -> torch.Tensor:
    """s_j = sum over (model, hookpoint) of the per-block decoder norm (Minder v(.))."""
    # W_dec_LMPD: [n_latents, n_models, n_hookpoints, d_model]
    return inner.W_dec_LMPD.detach().norm(dim=-1).sum(dim=(1, 2))  # -> [n_latents]


@torch.no_grad()
def estimate_threshold(inner, batches: Iterator[torch.Tensor], k: int, mode: str):
    """Return (log_threshold_L, theta, n_batches) for the requested mode.

    Per-batch threshold = the (k*B)-th largest value of the (scaled) pre-activation
    matrix flattened over [tokens x latents] = the minimum selected value, exactly
    what BatchTopK keeps. theta = mean of those over batches.
    """
    s_L = decoder_norm_scale_L(inner) if mode == "decoder_norm" else None
    if s_L is not None:
        s_L = s_L.clamp(min=1e-6)  # dead latents have ~0 decoder norm -> huge threshold (stays off)

    per_batch_min: list[torch.Tensor] = []
    for x_BMPD in batches:
        pre_BL = pre_activation_BL(inner, x_BMPD)
        scored_BL = pre_BL * s_L if s_L is not None else pre_BL
        B = scored_BL.shape[0]
        n_keep = k * B
        flat = scored_BL.reshape(-1)
        if n_keep >= flat.numel():
            continue
        # the n_keep-th largest value == minimum of the kept set
        kth_largest = flat.topk(n_keep, sorted=True).values[-1]
        per_batch_min.append(kth_largest.float().cpu())

    if not per_batch_min:
        raise RuntimeError("No usable batches; check --activations and --batch_tokens.")
    theta = torch.stack(per_batch_min).mean().item()
    n_batches = len(per_batch_min)

    n_latents = inner.n_latents
    if mode == "global":
        log_threshold_L = torch.full((n_latents,), float(torch.log(torch.tensor(theta))))
    elif mode == "decoder_norm":
        # keep raw f iff f * s_j > theta  <=>  f > theta / s_j
        thr_L = (theta / s_L.cpu()).clamp(min=1e-12)
        log_threshold_L = thr_L.log()
    else:
        raise ValueError(f"unknown mode {mode!r}")
    return log_threshold_L.float(), theta, n_batches


@torch.no_grad()
def mean_l0_under_jumprelu(inner, log_threshold_L: torch.Tensor, batch: torch.Tensor) -> float:
    """Sanity check: average nonzero latents per token under the new gate. Want ~= k."""
    pre_BL = pre_activation_BL(inner, batch)
    thr_L = log_threshold_L.exp().to(pre_BL.device)
    gated = (pre_BL > thr_L) * pre_BL
    return (gated != 0).float().sum(dim=-1).mean().item()


def write_converted_checkpoint(
    crosscoder_dir: Path, checkpoint: str, out_dir: Path,
    log_threshold_L: torch.Tensor, bandwidth: float, device: str,
):
    """Write a self-contained JumpReLU checkpoint dir (config rewritten, param injected).

    Loadable by score_proteingym_features.py with --inference_activation keep.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = yaml.unsafe_load((crosscoder_dir / "config.yaml").read_text())
    n_latents = cfg["n_latents"]
    cfg["activation_fn"] = {
        "classname": JUMPRELU_CLASSNAME,
        "cfg": {"size": n_latents, "bandwidth": bandwidth},
    }
    (out_dir / "config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))

    # Derive the EXACT state_dict key the JumpReLU threshold will have by scaffolding
    # the target model — do NOT hardcode it. crosscode registers the activation as the
    # private attr `_activation_fn` (like `_W_enc`/`_W_dec`), so the key is
    # `_activation_fn.log_threshold_L`, not `activation_fn.log_threshold_L`. Hardcoding
    # the wrong name makes load_state_dict(strict=False) silently drop it, leaving the
    # threshold at its init exp(1)=2.72 -> wrong activations (the load-time gotcha).
    from crosscode.models.acausal_crosscoder import ModelHookpointAcausalCrosscoder

    target = ModelHookpointAcausalCrosscoder._scaffold_from_cfg(cfg)
    thr_names = [n for n, _ in target.named_parameters() if n.endswith("log_threshold_L")]
    if len(thr_names) != 1:
        raise RuntimeError(f"expected exactly one log_threshold_L param in scaffold, got {thr_names}")
    base_key = thr_names[0]  # e.g. "_activation_fn.log_threshold_L"

    state_dict = torch.load(crosscoder_dir / checkpoint, map_location=device, weights_only=True)
    wrapped = any(key.startswith("crosscoder.") for key in state_dict)
    param_key = f"crosscoder.{base_key}" if wrapped else base_key
    # Drop any stale/misnamed threshold keys before injecting the correct one.
    for k in [k for k in state_dict if k.endswith("log_threshold_L")]:
        del state_dict[k]
    state_dict[param_key] = log_threshold_L.to(device)
    torch.save(state_dict, out_dir / checkpoint)
    log.info("Wrote converted JumpReLU checkpoint to %s (param key %s)", out_dir, param_key)

    # Round-trip verification: reload via load_sae and confirm the threshold actually
    # took (the in-memory mean_L0 check does NOT exercise the save/load path).
    from interplm.sae.inference import load_sae

    reloaded = load_sae(out_dir, device=device, model_name=checkpoint).crosscoder
    got = reloaded.activation_fn.log_threshold_L.exp()
    want = log_threshold_L.to(got.device).exp()
    if not torch.allclose(got, want, rtol=1e-3, atol=1e-4):
        raise RuntimeError(
            f"Round-trip check FAILED: reloaded threshold exp mean {got.mean():.4f} != "
            f"injected {want.mean():.4f}. The state_dict key '{param_key}' did not bind."
        )
    log.info("Round-trip OK: reloaded threshold exp mean=%.4f matches injected", got.mean().item())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--crosscoder_dir", type=Path, required=True)
    ap.add_argument("--checkpoint", type=str, default="model.pt")
    ap.add_argument("--activations", type=Path, default=None,
                    help="dir of [N, M, P, D] residual-stream activation .pt tensors "
                         "(omit if --fasta is given)")
    ap.add_argument("--fasta", type=Path, default=None,
                    help="FASTA / one-seq-per-line file to embed on the fly (single self-contained "
                         "run: sequences in, threshold out). Mutually exclusive with --activations.")
    ap.add_argument("--model_name", type=str, default="Rostlab/prot_t5_xl_uniref50",
                    help="ProtT5 model for --fasta embedding (match training)")
    ap.add_argument("--seqs_per_shard", type=int, default=64,
                    help="sequences embedded per shard when using --fasta (memory knob)")
    ap.add_argument("--max_seq_len", type=int, default=512,
                    help="cap each sequence length when using --fasta (match train length cap)")
    ap.add_argument("--mode", choices=["global", "decoder_norm"], default="global")
    ap.add_argument("--k", type=int, default=32, help="k_per_example used at training")
    ap.add_argument("--batch_tokens", type=int, default=512,
                    help="token-batch size for the BatchTopK cutoff; match training batch_size")
    ap.add_argument("--max_batches", type=int, default=200,
                    help="cap batches used for estimation (theta converges quickly)")
    ap.add_argument("--bandwidth", type=float, default=2.0,
                    help="JumpReLU STE bandwidth stored in the config (unused at pure inference)")
    ap.add_argument("--out_dir", type=Path, default=None,
                    help="defaults to <crosscoder_dir>")
    ap.add_argument("--write_converted", action="store_true",
                    help="also write a self-contained converted JumpReLU checkpoint dir")
    ap.add_argument("--device", type=str, default=None, help="cuda|mps|cpu (auto if unset)")
    args = ap.parse_args()

    if args.device:
        device = args.device
    elif torch.cuda.is_available():
        device = "cuda"
    elif torch.backends.mps.is_available():
        device = "mps"
    else:
        device = "cpu"
    out_dir = args.out_dir or args.crosscoder_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.fasta is not None:
        if args.activations is not None:
            ap.error("pass either --activations or --fasta, not both")
        # +1 batch of headroom so estimation always has max_batches full batches available
        target_tokens = args.batch_tokens * (args.max_batches + 1)
        args.activations = embed_fasta_to_shards(
            args.fasta, out_dir / "_calib_acts", device,
            model_name=args.model_name, seqs_per_shard=args.seqs_per_shard,
            max_seq_len=args.max_seq_len, max_tokens=target_tokens,
        )
    elif args.activations is None:
        ap.error("provide --activations or --fasta")

    wrapper, inner = load_batchtopk_crosscoder(args.crosscoder_dir, args.checkpoint, device)
    dtype = next(inner.parameters()).dtype

    def capped_batches():
        for i, b in enumerate(iter_activation_batches(args.activations, args.batch_tokens, device, dtype)):
            if i >= args.max_batches:
                break
            yield b

    # estimation consumes the iterator; keep one batch for the L0 sanity check
    sanity_batch = next(iter_activation_batches(args.activations, args.batch_tokens, device, dtype))
    log_threshold_L, theta, n_batches = estimate_threshold(inner, capped_batches(), args.k, args.mode)
    mean_l0 = mean_l0_under_jumprelu(inner, log_threshold_L, sanity_batch)

    log.info("mode=%s  theta=%.6g  n_batches=%d  mean_L0=%.2f (target k=%d)",
             args.mode, theta, n_batches, mean_l0, args.k)
    if not (0.5 * args.k <= mean_l0 <= 2.0 * args.k):
        log.warning("mean_L0 %.2f is far from k=%d — threshold estimate looks off; "
                    "check --batch_tokens and the activation layout.", mean_l0, args.k)

    sidecar = out_dir / "jumprelu_threshold.pt"
    torch.save({
        "log_threshold_L": log_threshold_L.cpu(),
        "theta": theta, "mode": args.mode, "k": args.k,
        "batch_tokens": args.batch_tokens, "n_batches": n_batches,
        "bandwidth": args.bandwidth, "mean_L0_check": mean_l0,
    }, sidecar)
    log.info("Wrote threshold sidecar to %s", sidecar)

    if args.write_converted:
        write_converted_checkpoint(
            args.crosscoder_dir, args.checkpoint, out_dir / "converted",
            log_threshold_L, args.bandwidth, device,
        )
        shutil.copy2(sidecar, out_dir / "converted" / "jumprelu_threshold.pt")


if __name__ == "__main__":
    main()
