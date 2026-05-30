"""
Diagnostic: compare alternative readouts for the ProteinGym crosscoder eval.

The production eval (score_proteingym_features.py) uses ONE readout — signed
Δ of the gated (post-topk) activation at the single mutated residue — and it
came out weak (median |Spearman| ≈ 0.03). Audits showed the cause is *coverage
collapse*: the median candidate feature fires at the mutated position in only
~1% of variants, so the gated Δ is a structural zero almost everywhere.

This script recomputes, in a single forward pass per sequence, every readout we
want to compare, on a small assay set, so we can see which one actually carries
signal before committing to a full re-run.

Readouts computed per (assay, concept, candidate-feature), correlated with DMS:

  GATED (post-topk, the production family)
    g_delta_signed   a_mut[p] - a_wt[p]            (PRODUCTION metric)
    g_abs            a_mut[p]                       (H-E1.1: absolute activation)
    g_pool_max       max_r |a_mut[r] - a_wt[r]|     (H-E1.2: sequence-pooled Δ)
    g_pool_mean      mean_r |a_mut[r] - a_wt[r]|

  DENSE (pre-topk, read the candidate feature's own pre-activation densely;
         this restores coverage while still looking at ONE named direction — the
         surgical version of H-T1, NOT a full 8192-dense readout)
    d_delta_signed   z_mut[p] - z_wt[p]
    d_abs            z_mut[p]
    d_pool_max       max_r |z_mut[r] - z_wt[r]|
    d_pool_mean      mean_r |z_mut[r] - z_wt[r]|

  RECONSTRUCTION (feature-independent; H-E1.3)
    recon_delta      e_mut[p] - e_wt[p]   where e = ||x - x_hat|| per residue
    recon_pool_max   max_r |e_mut[r] - e_wt[r]|

  RAW ProtT5 BASELINE (feature-independent; H-T2 — the crucial control)
    raw_delta        ||x_mut[p] - x_wt[p]||   (summed over 24 layers)
    raw_pool_max     max_r ||x_mut[r] - x_wt[r]||

For every readout we report the better-aligned of:
    rho_dir = Spearman(metric, DMS_score)
    rho_mag = Spearman(|metric|, |DMS_score - median|)
and surface |rho| so directional and magnitude framings are comparable.

The raw-ProtT5 baseline is the control that decides interpretation: if a dense
or pooled readout beats raw, the crosscoder concentrates signal usefully; if it
merely matches raw, we are effectively probing the PLM through a learned basis.

Usage (cluster, inside the same env as the production run):
    python diagnose_proteingym_metrics.py \\
        --crosscoder_dir <ckpt> --checkpoint ae_normalized.pt \\
        --pairings <heldout_all_top_pairings.csv> \\
        --matches <concept_matches.csv> \\
        --dms_dir <DMS_ProteinGym_substitutions> \\
        --output_dir <out> --max_variants 3000
"""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
from scipy.stats import ConstantInputWarning, spearmanr

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "repos" / "InterPLM"))
sys.path.insert(0, str(_REPO_ROOT / "repos" / "crosscode"))

from interplm.embedders.prott5 import ProtT5CrosscoderEmbedder  # noqa: E402
from interplm.sae.inference import load_sae  # noqa: E402

# Reuse the exact candidate-feature selection + mutant parsing from production
from score_proteingym_features import (  # noqa: E402
    load_candidate_features,
    parse_mutant_positions,
)

warnings.filterwarnings("ignore", category=ConstantInputWarning)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("diagnose_pg")


# Small, information-dense default set: a known strong hit, a dense-recoverable
# one, a gated-silent functional domain, a disordered case, plus two controls.
DEFAULT_ASSAYS = [
    ("PPARG_HUMAN_Majithia_2016", "Domain_NR LBD"),          # strong under abs/dense
    ("DNJA1_HUMAN_Tsuboyama_2023_2LO1", "Domain_J"),          # gated~0, abs=-0.42
    ("MK01_HUMAN_Brenan_2016", "Domain_Protein kinase"),      # gated-silent (93% zero)
    ("NUD15_HUMAN_Suiter_2020", "Domain_Nudix hydrolase"),    # functional, mid
    ("BBC1_YEAST_Tsuboyama_2023_1TG0", "Domain_SH3"),         # stability, SH3
    ("TAT_HV1BR_Fernandes_2016", "Region_Disordered"),        # disordered baseline
]


@torch.no_grad()
def encode_full(
    embedder: ProtT5CrosscoderEmbedder,
    crosscoder,
    sequence: str,
    feat_idx: List[int],
    device: torch.device,
) -> Dict[str, np.ndarray]:
    """One forward pass; return all per-residue quantities we need.

    Returns dict with:
      gated[L, F]   post-topk activation for candidate features
      dense[L, F]   relu(pre-bias + b_enc) for candidate features (pre-gate)
      recon[L]      per-residue reconstruction error ||x - x_hat||
      raw[L, 24]    per-residue, per-layer L2 norm of x (for raw-ProtT5 baseline)
    """
    cc = crosscoder.crosscoder  # underlying BaseCrosscoder
    dtype = next(cc.parameters()).dtype

    # x: [L, M=1, P=24, D=1024]
    x = embedder.extract_embeddings([sequence], batch_size=1).to(device=device, dtype=dtype)

    # Dense pre-gate activations (linear projection + encoder bias)
    pre = cc.get_pre_bias_BL(x)
    if cc.b_enc_L is not None:
        pre = pre + cc.b_enc_L
    dense_all = torch.relu(pre)                      # [L, 8192]  pre-gate, ReLU'd
    gated_all = cc.activation_fn.forward(pre)        # [L, 8192]  post-topk

    # Reconstruction from the GATED latents (the only thing the decoder ever sees)
    x_hat = cc.decode_BMPD(gated_all)                # [L, 1, 24, 1024]
    recon = torch.linalg.vector_norm(
        (x - x_hat).reshape(x.shape[0], -1), dim=1
    )                                                # [L]

    # Raw ProtT5 per-residue, per-layer norm — kept as x itself for Δ later.
    # We return the raw embedding tensor (moved to cpu) for difference math.
    return {
        "gated": gated_all[:, feat_idx].float().cpu().numpy(),
        "dense": dense_all[:, feat_idx].float().cpu().numpy(),
        "recon": recon.float().cpu().numpy(),
        "raw": x.squeeze(1).float().cpu().numpy(),   # [L, 24, 1024]
    }


def _best_abs_spearman(metric: np.ndarray, dms: np.ndarray) -> Tuple[float, float, float]:
    """Return (best_abs_rho, rho_dir, rho_mag)."""
    if len(metric) < 5:
        return np.nan, np.nan, np.nan
    med = np.median(dms)
    try:
        rho_dir = spearmanr(metric, dms, nan_policy="omit").statistic
    except Exception:
        rho_dir = np.nan
    try:
        rho_mag = spearmanr(np.abs(metric), np.abs(dms - med), nan_policy="omit").statistic
    except Exception:
        rho_mag = np.nan
    cand = [r for r in (rho_dir, rho_mag) if r is not None and not np.isnan(r)]
    best = max(cand, key=abs) if cand else np.nan
    return best, rho_dir, rho_mag


def diagnose_assay(
    dms_id: str,
    concept: str,
    dms_dir: Path,
    pairings_path: Path,
    embedder: ProtT5CrosscoderEmbedder,
    crosscoder,
    device: torch.device,
    max_variants: Optional[int],
) -> Tuple[List[dict], Optional[dict]]:
    candidates = load_candidate_features(pairings_path, concept)
    if not candidates:
        log.warning("  no candidates for %s", concept)
        return [], None
    feat_ids = [f for f, _ in candidates]
    feat_f1pd = {f: f1 for f, f1 in candidates}
    feat_pos = {f: i for i, f in enumerate(feat_ids)}

    dms_path = dms_dir / f"{dms_id}.csv"
    if not dms_path.exists():
        log.error("  missing %s", dms_path)
        return [], None
    dms = pd.read_csv(dms_path)
    if max_variants and len(dms) > max_variants:
        dms = dms.sample(max_variants, random_state=42).reset_index(drop=True)

    # WT sequence by reverting the first variant's mutation(s)
    first = parse_mutant_positions(str(dms["mutant"].iloc[0]))
    if not first:
        return [], None
    wt = list(dms["mutated_sequence"].iloc[0])
    for waa, p, _ in first:
        wt[p - 1] = waa
    wt_str = "".join(wt)

    log.info("  WT len=%d, %d candidate feats, %d variants", len(wt_str), len(feat_ids), len(dms))
    wt_enc = encode_full(embedder, crosscoder, wt_str, feat_ids, device)

    # Accumulators: one list per readout, per feature (feature-dependent) +
    # feature-independent (recon, raw) per variant.
    n = len(dms)
    F = len(feat_ids)
    # feature-dependent metrics: arrays [n_variants, F]
    acc = {k: np.full((n, F), np.nan) for k in
           ["g_delta_signed", "g_abs", "g_pool_max", "g_pool_mean",
            "d_delta_signed", "d_abs", "d_pool_max", "d_pool_mean"]}
    # feature-independent: arrays [n_variants]
    acc_fi = {k: np.full(n, np.nan) for k in
              ["recon_delta", "recon_pool_max", "raw_delta", "raw_pool_max"]}
    dms_scores = np.full(n, np.nan)

    Lwt = wt_enc["gated"].shape[0]
    for vi, (_, vrow) in enumerate(dms.iterrows()):
        parsed = parse_mutant_positions(str(vrow["mutant"]))
        if not parsed:
            continue
        p1 = parsed[0][1]  # first mutated position (1-based)
        if p1 - 1 >= Lwt:
            continue
        mut_enc = encode_full(embedder, crosscoder, vrow["mutated_sequence"], feat_ids, device)
        Lm = mut_enc["gated"].shape[0]
        L = min(Lwt, Lm)
        if p1 - 1 >= L:
            continue
        dms_scores[vi] = float(vrow["DMS_score"])

        # ---- feature-dependent ----
        g_wt, g_mut = wt_enc["gated"][:L], mut_enc["gated"][:L]
        d_wt, d_mut = wt_enc["dense"][:L], mut_enc["dense"][:L]
        gd = g_mut - g_wt           # [L, F]
        dd = d_mut - d_wt
        acc["g_delta_signed"][vi] = gd[p1 - 1]
        acc["g_abs"][vi] = g_mut[p1 - 1]
        acc["g_pool_max"][vi] = np.abs(gd).max(axis=0)
        acc["g_pool_mean"][vi] = np.abs(gd).mean(axis=0)
        acc["d_delta_signed"][vi] = dd[p1 - 1]
        acc["d_abs"][vi] = d_mut[p1 - 1]
        acc["d_pool_max"][vi] = np.abs(dd).max(axis=0)
        acc["d_pool_mean"][vi] = np.abs(dd).mean(axis=0)

        # ---- reconstruction (feature-independent) ----
        e_wt, e_mut = wt_enc["recon"][:L], mut_enc["recon"][:L]
        ed = e_mut - e_wt
        acc_fi["recon_delta"][vi] = ed[p1 - 1]
        acc_fi["recon_pool_max"][vi] = np.abs(ed).max()

        # ---- raw ProtT5 (feature-independent) ----
        # per-residue L2 over 24 layers x 1024 dims
        rdiff = (mut_enc["raw"][:L] - wt_enc["raw"][:L]).reshape(L, -1)
        rnorm = np.linalg.norm(rdiff, axis=1)   # [L]
        acc_fi["raw_delta"][vi] = rnorm[p1 - 1]
        acc_fi["raw_pool_max"][vi] = rnorm.max()

        if vi % 500 == 0:
            log.info("    variant %d/%d", vi + 1, n)

    valid = ~np.isnan(dms_scores)
    dms_v = dms_scores[valid]

    # Per-feature rows
    rows = []
    for f in feat_ids:
        fi = feat_pos[f]
        row = {"DMS_id": dms_id, "concept": concept, "feature": f,
               "f1_per_domain": feat_f1pd[f], "n_variants": int(valid.sum()),
               "fire_rate_gated": float((acc["g_abs"][valid, fi] > 1e-6).mean()),
               "fire_rate_dense": float((acc["d_abs"][valid, fi] > 1e-6).mean())}
        for key, arr in acc.items():
            best, rdir, rmag = _best_abs_spearman(arr[valid, fi], dms_v)
            row[key] = best
        rows.append(row)

    # Feature-independent (one row per assay)
    fi_row = {"DMS_id": dms_id, "concept": concept, "n_variants": int(valid.sum())}
    for key, arr in acc_fi.items():
        best, rdir, rmag = _best_abs_spearman(arr[valid], dms_v)
        fi_row[key] = best
    return rows, fi_row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crosscoder_dir", type=Path, required=True)
    ap.add_argument("--checkpoint", type=str, default="ae_normalized.pt")
    ap.add_argument("--pairings", type=Path, required=True)
    ap.add_argument("--matches", type=Path, required=True)
    ap.add_argument("--dms_dir", type=Path, required=True)
    ap.add_argument("--output_dir", type=Path, required=True)
    ap.add_argument("--assays", type=Path, default=None,
                    help="Optional TSV/CSV (DMS_id,concept). Default: built-in 6-assay set.")
    ap.add_argument("--max_variants", type=int, default=3000)
    ap.add_argument("--device", type=str, default=None)
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.assays:
        adf = pd.read_csv(args.assays, sep="\t" if args.assays.suffix == ".tsv" else ",")
        assays = list(zip(adf["DMS_id"], adf["concept"]))
    else:
        assays = DEFAULT_ASSAYS

    if args.device:
        device = torch.device(args.device)
    elif torch.cuda.is_available():
        device = torch.device("cuda")
    elif torch.backends.mps.is_available():
        device = torch.device("mps")
    else:
        device = torch.device("cpu")
    log.info("Device: %s | %d assays", device, len(assays))

    embedder = ProtT5CrosscoderEmbedder(device=str(device))
    crosscoder = load_sae(args.crosscoder_dir, device=str(device), model_name=args.checkpoint)
    crosscoder.eval()

    all_feat_rows, all_fi_rows = [], []
    for dms_id, concept in assays:
        log.info("=== %s | %s ===", dms_id, concept)
        rows, fi = diagnose_assay(
            dms_id, concept, args.dms_dir, args.pairings,
            embedder, crosscoder, device, args.max_variants,
        )
        all_feat_rows.extend(rows)
        if fi:
            all_fi_rows.append(fi)
        # incremental save
        if all_feat_rows:
            pd.DataFrame(all_feat_rows).to_csv(args.output_dir / "metric_comparison_per_feature.csv", index=False)
        if all_fi_rows:
            pd.DataFrame(all_fi_rows).to_csv(args.output_dir / "metric_comparison_baselines.csv", index=False)

    feat = pd.DataFrame(all_feat_rows)
    fi = pd.DataFrame(all_fi_rows)
    if feat.empty:
        log.warning("no rows"); return

    readouts_g = ["g_delta_signed", "g_abs", "g_pool_max", "g_pool_mean"]
    readouts_d = ["d_delta_signed", "d_abs", "d_pool_max", "d_pool_mean"]

    # Best feature per (assay, concept) under each readout
    log.info("\n=== BEST |Spearman| per (assay, concept), by readout ===")
    log.info("(production metric = g_delta_signed)")
    hdr = f"{'assay':40s} {'concept':22s}" + "".join(f"{r:>14s}" for r in readouts_g + readouts_d)
    log.info(hdr)
    for (dms_id, concept), g in feat.groupby(["DMS_id", "concept"]):
        cells = "".join(f"{g[r].abs().max():>14.3f}" for r in readouts_g + readouts_d)
        log.info(f"{dms_id[:39]:40s} {concept[:21]:22s}{cells}")

    log.info("\n=== FEATURE-INDEPENDENT BASELINES (best |Spearman| per assay) ===")
    log.info(f"{'assay':40s} {'concept':22s}{'recon_delta':>14s}{'recon_pool':>14s}{'raw_delta':>14s}{'raw_pool':>14s}")
    for _, r in fi.iterrows():
        log.info(f"{r['DMS_id'][:39]:40s} {r['concept'][:21]:22s}"
                 f"{r['recon_delta']:>14.3f}{r['recon_pool_max']:>14.3f}"
                 f"{r['raw_delta']:>14.3f}{r['raw_pool_max']:>14.3f}")

    # Aggregate: median of (best-feature |rho|) across assays, per readout
    log.info("\n=== AGGREGATE: median of (best-feature |Spearman|) across assays ===")
    for r in readouts_g + readouts_d:
        med = feat.groupby(["DMS_id", "concept"])[r].apply(lambda s: s.abs().max()).median()
        log.info(f"  {r:18s}  median best |rho| = {med:.3f}")
    for r in ["recon_delta", "recon_pool_max", "raw_delta", "raw_pool_max"]:
        log.info(f"  {r:18s}  median |rho|       = {fi[r].abs().median():.3f}")

    log.info("\nWrote metric_comparison_per_feature.csv and metric_comparison_baselines.csv to %s", args.output_dir)


if __name__ == "__main__":
    sys.exit(main())
