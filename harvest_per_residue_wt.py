"""Per-residue MotifAE-comparable readout — WT-only, all 8192 features.

Replicates MotifAE's (Hu et al. 2025) feature-vs-DMS setup as closely as we can, to give an
apples-to-apples comparison instead of our concept-matched-only number:

  - per RESIDUE (not per variant): correlate a feature's WT activation a[p] against the per-residue
    mean mutation effect m[p] (mean DMS_score over substitutions at p, positions with >= MIN_SUBS);
  - over ALL 8192 features (best-of-N per assay), the same best-of-many MotifAE does — not just the
    handful of concept-matched candidates;
  - with MotifAE's per-protein activity filter (feature active on >= ACTIVE_FRAC of residues);
  - stratified live / dead / paired so the dead features give the chance floor (if the live best is
    not above the dead best, the number is just max-over-N selection noise).

Cheap: ONE ProtT5 + crosscoder forward per protein (WT only; no mutant forwards) → (L_wt, 8192)
activations. ~135 forwards. Runs locally on an M1 (MPS) — the readout needs no per-variant scan.

Reads:
  --matches    data/proteingym/concept_matches.csv
  --dms_dir    data/DMS_ProteinGym_substitutions/
  --crosscoder_dir / --checkpoint   the trained crosscoder (or the baseline, for the null)
  --maxact     <ckpt>/.../max_activations_per_feature.pt   (dead mask, for the analysis)
  --pairings   heldout_all_top_pairings.csv                (paired feature ids, for the analysis)
Writes:
  <output_dir>/<DMS_id>__per_res_wt.parquet  (8192 rows: feature, spearman, frac_active, n_pos)
  <output_dir>/per_residue_wt_spearman.csv   (aggregate)
  data/proteingym/per_residue_wt_summary.csv (stratified best-of-N medians + printed table)

Usage (local M1, from the crosscoder root):
  uv run --project repos/crosscode --with scipy --with pyarrow python \
    repos/sparse-crosscoders-prott5/harvest_per_residue_wt.py \
      --crosscoder_dir model_checkpoints/crosscoder_l8192_k32_bs512_full_2026-03-12_06-03-41/crashed_epoch_0_step_2519836 \
      --matches data/proteingym/concept_matches.csv \
      --dms_dir data/DMS_ProteinGym_substitutions \
      --maxact model_checkpoints/crosscoder_l8192_k32_bs512_full_2026-03-12_06-03-41/crashed_epoch_0_step_2519836/uniprotkb_modern_score45_67k/max_activations_per_feature.pt \
      --pairings data/crosscoder_eval/uniprotkb_modern_score45_67k/test_counts/heldout_all_top_pairings.csv \
      --output_dir data/proteingym/per_residue_wt --skip_existing
"""
from __future__ import annotations

import argparse, gc, logging, re, sys, warnings
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import torch
from scipy.stats import ConstantInputWarning, spearmanr

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "repos" / "InterPLM"))
sys.path.insert(0, str(_REPO_ROOT / "repos" / "crosscode"))
from interplm.embedders.prott5 import ProtT5CrosscoderEmbedder  # noqa: E402
from interplm.sae.inference import load_sae  # noqa: E402

warnings.filterwarnings("ignore", category=ConstantInputWarning)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("per_res_wt")

MUTANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z\*])$")
MIN_SUBS = 10      # MotifAE: residues with >= 10 assayed substitutions
MIN_POS = 5        # require >= 5 usable positions for a per-residue correlation
ACTIVE_FRAC = 0.30  # MotifAE: features active on >= 30% of a protein's residues


def parse_mutant_positions(mutant_field: str) -> List[Tuple[str, int, str]]:
    out = []
    for token in str(mutant_field).split(":"):
        m = MUTANT_RE.match(token.strip())
        if not m:
            return []
        out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


@torch.no_grad()
def encode_wt(crosscoder, x, device) -> np.ndarray:
    """Per-residue post-BatchTopK activations for ALL features. x: [L,1,24,1024] -> (L, 8192)."""
    inner = crosscoder.crosscoder
    pre = inner.get_pre_bias_BL(x.to(device=device, dtype=next(inner.parameters()).dtype))
    if inner.b_enc_L is not None:
        pre = pre + inner.b_enc_L
    return inner.activation_fn.forward(pre).float().cpu().numpy()


def harvest_one(dms_id, dms_dir, embedder, crosscoder, n_latents, device):
    dms_path = dms_dir / f"{dms_id}.csv"
    if not dms_path.exists():
        log.error("  DMS missing: %s", dms_path); return None
    dms = pd.read_csv(dms_path)
    dms = dms[dms["mutant"].apply(lambda s: len(parse_mutant_positions(s)) == 1)].copy()
    if len(dms) == 0:
        return None
    # WT reconstruction from the first variant
    first = parse_mutant_positions(str(dms["mutant"].iloc[0]))
    wt = list(dms["mutated_sequence"].iloc[0])
    for waa, p, _ in first:
        wt[p - 1] = waa
    wt_str = "".join(wt)
    L = len(wt_str)

    # per-residue mean mutation effect m[p] (1-indexed position -> mean DMS over substitutions)
    dms["pos"] = dms["mutant"].apply(lambda s: parse_mutant_positions(s)[0][1])
    eff = dms.groupby("pos")["DMS_score"].agg(["mean", "size"])
    good = eff.index[(eff["size"] >= MIN_SUBS) & (eff.index >= 1) & (eff.index <= L)]
    if len(good) < MIN_POS:
        return None
    m = eff.loc[good, "mean"].values
    pos0 = (good.values - 1)  # 0-indexed residue positions

    wt_act = encode_wt(crosscoder, embedder.extract_embeddings([wt_str], batch_size=1), device)  # (L, 8192)
    A = wt_act[pos0, :]            # (n_pos, 8192) activations at the scored positions
    frac_active = (wt_act > 0).mean(axis=0)  # over ALL residues (MotifAE >=30% filter basis)

    # Vectorized Spearman = Pearson on ranks (average ties), m fixed across all features.
    rm = pd.Series(m).rank().to_numpy(dtype=np.float64)
    rm_c = rm - rm.mean()
    rm_norm = np.sqrt((rm_c ** 2).sum())
    rA = pd.DataFrame(A.astype(np.float64)).rank().to_numpy()   # (n_pos, 8192), avg ties per column
    rA_c = rA - rA.mean(axis=0, keepdims=True)
    denom = np.sqrt((rA_c ** 2).sum(axis=0)) * rm_norm          # 0 for constant features / constant m
    with np.errstate(invalid="ignore", divide="ignore"):
        rhos = (rA_c.T @ rm_c) / denom                          # (8192,)
    rhos = np.where(np.isfinite(rhos), rhos, np.nan).astype(np.float32)
    return pd.DataFrame({"DMS_id": dms_id, "feature": np.arange(n_latents),
                         "spearman": rhos, "frac_active": frac_active.astype(np.float32),
                         "n_pos": int(len(good))})


def analyze(output_dir: Path, maxact_path: Path, pairings_path: Path) -> None:
    parts = sorted(output_dir.glob("*__per_res_wt.parquet"))
    if not parts:
        log.warning("nothing to analyze"); return
    ff = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
    ff.to_csv(output_dir / "per_residue_wt_spearman.csv", index=False)
    ff["rho"] = ff["spearman"].abs().fillna(0.0)

    maxact = torch.load(maxact_path, map_location="cpu").numpy()
    dead = set(np.where(maxact == 0)[0].tolist())
    paired = set(pd.read_csv(pairings_path)["feature"].unique())

    def is_live(f): return f not in dead
    rows = []
    for dms, g in ff.groupby("DMS_id"):
        feat_rho = g.set_index("feature")["rho"]
        feat_act = g.set_index("feature")["frac_active"]
        live = [f for f in feat_rho.index if is_live(f)]
        active30 = [f for f in feat_rho.index if feat_act[f] >= ACTIVE_FRAC]

        def best(feats):
            s = feat_rho.reindex(feats).dropna()
            return float(s.max()) if len(s) else np.nan
        rows.append({
            "DMS_id": dms,
            "all_live_best":        best([f for f in live]),
            "live_unpaired_best":   best([f for f in live if f not in paired]),
            "live_paired_best":     best([f for f in live if f in paired]),
            "dead_best":            best(sorted(dead)),
            "active30_best":        best([f for f in active30 if is_live(f)]),          # MotifAE pool
            "active30_paired_best": best([f for f in active30 if f in paired]),
            "active30_unpaired_best": best([f for f in active30 if is_live(f) and f not in paired]),
            "n_active30": int(len([f for f in active30 if is_live(f)])),
        })
    out = pd.DataFrame(rows)
    out.to_csv(_REPO_ROOT / "data/proteingym/per_residue_wt_summary.csv", index=False)
    print(f"\nper-residue WT, {len(out)} assays — median best-feature-per-assay |rho|:")
    print(f"  [MotifAE per-residue best-match, ESM2: 0.41; plain SAE: 0.33]")
    for c in ["active30_best", "active30_paired_best", "active30_unpaired_best",
              "all_live_best", "live_paired_best", "live_unpaired_best", "dead_best"]:
        v = out[c].dropna()
        print(f"  {c:<24s} {v.median():.3f}   (n={len(v)}, >0.3={int((v>0.3).sum())}, >0.5={int((v>0.5).sum())})")
    print(f"  median # live features active on >=30% of residues: {int(out['n_active30'].median())}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crosscoder_dir", type=Path, required=True)
    ap.add_argument("--checkpoint", type=str, default="ae_normalized.pt")
    ap.add_argument("--matches", type=Path, required=True)
    ap.add_argument("--dms_dir", type=Path, required=True)
    ap.add_argument("--maxact", type=Path, required=True)
    ap.add_argument("--pairings", type=Path, required=True)
    ap.add_argument("--output_dir", type=Path, required=True)
    ap.add_argument("--max_seq_len", type=int, default=2048)
    ap.add_argument("--device", type=str, default=None)
    ap.add_argument("--skip_existing", action="store_true")
    ap.add_argument("--analyze_only", action="store_true")
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not args.analyze_only:
        matches = pd.read_csv(args.matches)
        if args.max_seq_len:
            matches = matches[matches["target_seq_len"] <= args.max_seq_len]
        assay_ids = sorted(matches["DMS_id"].unique())

        if args.device:
            device = torch.device(args.device)
        elif torch.cuda.is_available():
            device = torch.device("cuda")
        elif torch.backends.mps.is_available():
            device = torch.device("mps")
        else:
            device = torch.device("cpu")
        log.info("Device: %s; %d assays", device, len(assay_ids))
        embedder = ProtT5CrosscoderEmbedder(device=str(device))
        crosscoder = load_sae(args.crosscoder_dir, device=str(device), model_name=args.checkpoint)
        crosscoder.eval()
        n_latents = crosscoder.dict_size

        for i, dms_id in enumerate(assay_ids):
            out_path = args.output_dir / f"{dms_id}__per_res_wt.parquet"
            if args.skip_existing and out_path.exists():
                continue
            log.info("[%d/%d] %s", i + 1, len(assay_ids), dms_id)
            try:
                df = harvest_one(dms_id, args.dms_dir, embedder, crosscoder, n_latents, device)
            except Exception as e:
                log.exception("  failed: %s", e); continue
            if df is not None:
                df.to_parquet(out_path, index=False)
            gc.collect()

    analyze(args.output_dir, args.maxact, args.pairings)


if __name__ == "__main__":
    sys.exit(main())
