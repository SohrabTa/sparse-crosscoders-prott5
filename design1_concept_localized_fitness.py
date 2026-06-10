"""
Design (1) — concept-localized fitness enrichment, scoped to point/site concepts,
with a firing-density-matched random-feature null. (experiment 03, the label-validating
ProteinGym test; see 03-proteingym.md "Open questions").

QUESTION. A feature labelled with a *point/site* concept (disulfide bond, zinc finger,
glycosylation site, motif, modified residue) claims to mark a functional/structural site.
If the label is biologically valid, mutations at the residues where that feature FIRES
should be more consequential (larger |DMS effect|) than at residues where it does not.

STATISTIC (per assay x labelled feature), no best-of-N — the feature is fixed by its label:
  - m[p] = mean DMS_score over substitutions at position p (>= MIN_SUBS subs), as in
    harvest_per_residue_wt.py. effect magnitude dev[p] = |m[p] - median(m)| (sign-robust:
    works whether the assay scores gain or loss of fitness).
  - firing set = positions (among the scored ones) where the feature's WT activation a[p] > 0
    (BatchTopK emits exact zeros off the top-k). Point/site features fire sparsely.
  - localization AUC = P(dev at a firing position > dev at a non-firing position)
    = Mann-Whitney U / (k*(n-k)). AUC > 0.5 => the label sits on functionally important residues.

NULLS (the whole point — the old readout had none):
  A) firing-density-matched random LIVE features: among live features, take the ones whose
     firing-count k_g on THIS protein's scored positions is closest to the labelled feature's
     k; compute each one's AUC. Empirical p = frac(AUC_null >= AUC_obs), plus a z-score. This
     removes the "denser features correlate more" confound (exp 03, 2026-06-10: firing density
     is the real within-class driver, rho up to +0.53).
  B) within-sequence mask permutation: draw k random positions as a pseudo-firing set, AUC,
     repeat N_PERM. Controls the protein's own dev distribution.

NOTE / TODO (not in smoke): conservation matching. Functional residues are intrinsically more
mutation-sensitive; null A partially controls it (real + random features share the protein's
conservation landscape) but a positions-matched-on-entropy null would be stricter. Flag in the doc.

INPUTS (read-only):
  --crosscoder_dir   <ckpt>/  (load_sae -> BatchTopK crosscoder)        [ae_normalized.pt]
  --matches          data/proteingym/concept_matches.csv  (labelled features per (assay, concept))
  --dms_dir          data/DMS_ProteinGym_substitutions/<DMS_id>.csv
  --maxact           <ckpt>/.../max_activations_per_feature.pt  (live mask for null A)

OUTPUT:
  <output_dir>/design1_per_feature.csv  (one row per (assay, concept, feature): auc_obs, nullA_*, nullB_*)

Smoke (local M1, from crosscoder root):
  uv run --project repos/crosscode --with scipy --with pyarrow python \
    repos/sparse-crosscoders-prott5/design1_concept_localized_fitness.py \
      --crosscoder_dir model_checkpoints/crosscoder_l8192_k32_bs512_full_2026-03-12_06-03-41/crashed_epoch_0_step_2519836 \
      --matches data/proteingym/concept_matches.csv \
      --dms_dir data/DMS_ProteinGym_substitutions \
      --maxact model_checkpoints/crosscoder_l8192_k32_bs512_full_2026-03-12_06-03-41/crashed_epoch_0_step_2519836/uniprotkb_modern_score45_67k/max_activations_per_feature.pt \
      --output_dir data/proteingym/design1_localized \
      --assays A0A2Z5U3Z0_9INFA_Doud_2016,C6KNH7_9INFA_Lee_2018 --smoke

Disposition: COMMIT — this is the Design (1) reproducibility script. The headline run is on the
AuxK-fixed crosscoder + its new pairings (experiment 05); the --smoke run on the bugged model is
a pipeline check only, not a reported result.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import torch
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
from harvest_per_residue_wt import (  # noqa: E402  reuse exact conventions
    MIN_POS, MIN_SUBS, encode_wt, parse_mutant_positions,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "InterPLM" / "src"))
from interplm.embedders.prott5 import ProtT5CrosscoderEmbedder  # noqa: E402
from interplm.sae.inference import load_sae  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("design1")

# point/site concepts = those with a 3D locus (sparse, localized firing). Distributed concepts
# (Region_Disordered, Coiled coil, Domain_*) are out of scope for this test (exp 03 locality cut).
POINT_SITE_PREFIXES = ("disulfide", "glycosylation", "zinc finger", "motif", "modified residue")
N_NULL_FEATURES = 200   # density-matched random live features per labelled feature
N_PERM = 2000           # mask-permutation draws


def is_point_site(concept: str) -> bool:
    return str(concept).lower().startswith(POINT_SITE_PREFIXES)


def auc_from_firing(rank_dev: np.ndarray, fire_mask: np.ndarray) -> float:
    """Mann-Whitney AUC = P(dev_firing > dev_nonfiring), via summed ranks. rank_dev: ranks of
    |effect| over all scored positions (avg ties). fire_mask: bool over the same positions."""
    n = len(rank_dev)
    k = int(fire_mask.sum())
    if k < 2 or k > n - 2:
        return np.nan
    R = float(rank_dev[fire_mask].sum())
    U = R - k * (k + 1) / 2.0
    return U / (k * (n - k))


def auc_for_counts(rank_dev_sorted_cumsum, rank_dev: np.ndarray, ks: np.ndarray,
                   rng: np.random.Generator, n_draws: int) -> np.ndarray:
    """Mask-permutation AUCs for a vector of k values (one draw each)."""
    n = len(rank_dev)
    out = np.full(len(ks), np.nan)
    for i, k in enumerate(ks):
        k = int(k)
        if k < 2 or k > n - 2:
            continue
        idx = rng.choice(n, size=k, replace=False)
        U = rank_dev[idx].sum() - k * (k + 1) / 2.0
        out[i] = U / (k * (n - k))
    return out


def load_live_mask(maxact_path: Path, n_latents: int) -> np.ndarray:
    obj = torch.load(maxact_path, map_location="cpu")
    if isinstance(obj, dict):
        obj = obj.get("max_activations", next(iter(obj.values())))
    arr = np.asarray(obj).reshape(-1)[:n_latents]
    return arr > 0.0


def labelled_point_site_features(matches: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in matches.iterrows():
        if not is_point_site(r["matched_concept"]):
            continue
        feats = []
        if pd.notna(r.get("candidate_features")):
            feats = [int(x) for x in str(r["candidate_features"]).split(";") if x != ""]
        elif pd.notna(r.get("top_feature")):
            feats = [int(r["top_feature"])]
        for f in feats:
            rows.append({"DMS_id": r["DMS_id"], "concept": r["matched_concept"], "feature": f})
    return pd.DataFrame(rows).drop_duplicates(["DMS_id", "concept", "feature"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--crosscoder_dir", type=Path, required=True)
    ap.add_argument("--checkpoint", type=str, default="ae_normalized.pt")
    ap.add_argument("--matches", type=Path, required=True)
    ap.add_argument("--dms_dir", type=Path, required=True)
    ap.add_argument("--maxact", type=Path, required=True)
    ap.add_argument("--output_dir", type=Path, required=True)
    ap.add_argument("--assays", type=str, default=None, help="comma-separated DMS_ids to restrict to")
    ap.add_argument("--max_seq_len", type=int, default=2048)
    ap.add_argument("--device", type=str, default=None)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--smoke", action="store_true", help="print extra diagnostics per feature")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    log.info("SEED=%d  N_NULL_FEATURES=%d  N_PERM=%d", args.seed, N_NULL_FEATURES, N_PERM)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    device = args.device or ("cuda" if torch.cuda.is_available()
                             else "mps" if torch.backends.mps.is_available() else "cpu")
    log.info("device=%s", device)

    matches = pd.read_csv(args.matches)
    labelled = labelled_point_site_features(matches)
    if args.assays:
        keep = set(args.assays.split(","))
        labelled = labelled[labelled["DMS_id"].isin(keep)]
    log.info("point/site labelled features: %d rows, %d assays",
             len(labelled), labelled["DMS_id"].nunique())
    if labelled.empty:
        log.error("no point/site labelled features for the requested assays"); return

    embedder = ProtT5CrosscoderEmbedder(device=str(device))
    crosscoder = load_sae(args.crosscoder_dir, device=str(device), model_name=args.checkpoint)
    n_latents = crosscoder.crosscoder.W_dec_LMD.shape[0] if hasattr(crosscoder.crosscoder, "W_dec_LMD") else 8192
    live_mask = load_live_mask(args.maxact, n_latents)
    live_idx = np.where(live_mask)[0]
    log.info("n_latents=%d  live=%d", n_latents, len(live_idx))

    results = []
    for dms_id, grp in labelled.groupby("DMS_id"):
        dms_path = args.dms_dir / f"{dms_id}.csv"
        if not dms_path.exists():
            log.warning("  %s: DMS missing", dms_id); continue
        dms = pd.read_csv(dms_path)
        dms = dms[dms["mutant"].apply(lambda s: len(parse_mutant_positions(s)) == 1)].copy()
        if dms.empty:
            continue
        first = parse_mutant_positions(str(dms["mutant"].iloc[0]))
        wt = list(dms["mutated_sequence"].iloc[0])
        for waa, p, _ in first:
            wt[p - 1] = waa
        wt_str = "".join(wt)
        L = len(wt_str)
        if L > args.max_seq_len:
            log.warning("  %s: len %d > max_seq_len, skip", dms_id, L); continue

        dms["pos"] = dms["mutant"].apply(lambda s: parse_mutant_positions(s)[0][1])
        eff = dms.groupby("pos")["DMS_score"].agg(["mean", "size"])
        good = eff.index[(eff["size"] >= MIN_SUBS) & (eff.index >= 1) & (eff.index <= L)]
        if len(good) < MIN_POS:
            log.warning("  %s: <%d usable positions, skip", dms_id, MIN_POS); continue
        m = eff.loc[good, "mean"].values
        pos0 = good.values - 1
        dev = np.abs(m - np.median(m))
        rank_dev = rankdata(dev)                       # 1..n_pos, avg ties
        n_pos = len(good)

        wt_act = encode_wt(crosscoder, embedder.extract_embeddings([wt_str], batch_size=1), device)
        A = wt_act[pos0, :]                             # (n_pos, 8192) at scored positions
        fire_counts = (A > 0).sum(axis=0)               # k_g for every feature on these positions

        # density-matched live pool, recomputed per labelled feature's k below
        live_k = fire_counts[live_idx]

        log.info("  %s: L=%d  n_pos=%d  labelled=%d", dms_id, L, n_pos, len(grp))
        for _, row in grp.iterrows():
            f = int(row["feature"])
            fire = A[:, f] > 0
            k = int(fire.sum())
            auc_obs = auc_from_firing(rank_dev, fire)
            rec = {"DMS_id": dms_id, "concept": row["concept"], "feature": f,
                   "n_pos": n_pos, "k_fire": k, "auc_obs": auc_obs}
            if not np.isfinite(auc_obs):
                rec.update(nullA_n=0, nullA_mean=np.nan, nullA_p=np.nan, nullA_z=np.nan,
                           nullB_p95=np.nan, nullB_p=np.nan, status="degenerate_k")
                results.append(rec);
                if args.smoke: log.info("    f%d concept=%s k=%d  AUC=nan (degenerate)", f, row["concept"], k)
                continue

            # Null A: nearest-k live features (exclude the labelled feature itself)
            order = np.argsort(np.abs(live_k - k))
            pool = [int(live_idx[j]) for j in order if int(live_idx[j]) != f][:N_NULL_FEATURES]
            aucA = np.array([auc_from_firing(rank_dev, A[:, g] > 0) for g in pool])
            aucA = aucA[np.isfinite(aucA)]
            nullA_mean = float(np.mean(aucA)) if len(aucA) else np.nan
            nullA_p = float((aucA >= auc_obs).mean()) if len(aucA) else np.nan
            nullA_z = float((auc_obs - nullA_mean) / (np.std(aucA) + 1e-9)) if len(aucA) > 1 else np.nan
            realized_k = float(np.median(fire_counts[pool])) if pool else np.nan

            # Null B: mask permutation at fixed k
            aucB = auc_for_counts(None, rank_dev, np.full(N_PERM, k), rng, N_PERM)
            aucB = aucB[np.isfinite(aucB)]
            nullB_p95 = float(np.percentile(aucB, 95)) if len(aucB) else np.nan
            nullB_p = float((aucB >= auc_obs).mean()) if len(aucB) else np.nan

            rec.update(nullA_n=len(aucA), nullA_mean=nullA_mean, nullA_p=nullA_p, nullA_z=nullA_z,
                       nullA_median_k=realized_k, nullB_p95=nullB_p95, nullB_p=nullB_p, status="ok")
            results.append(rec)
            if args.smoke:
                log.info("    f%d %s  k=%d/%d  AUC=%.3f  nullA mean=%.3f p=%.3f z=%.2f (k~%.0f)  nullB p95=%.3f p=%.3f",
                         f, row["concept"], k, n_pos, auc_obs, nullA_mean, nullA_p, nullA_z,
                         realized_k, nullB_p95, nullB_p)

    out = pd.DataFrame(results)
    out_path = args.output_dir / "design1_per_feature.csv"
    out.to_csv(out_path, index=False)
    log.info("wrote %s (%d rows)", out_path, len(out))
    ok = out[out["status"] == "ok"] if "status" in out else out
    if len(ok):
        log.info("\n=== summary (status=ok, %d rows) ===", len(ok))
        log.info("median AUC_obs=%.3f  median nullA_mean=%.3f  frac nullA_p<0.05=%.2f  frac nullB_p<0.05=%.2f",
                 ok["auc_obs"].median(), ok["nullA_mean"].median(),
                 (ok["nullA_p"] < 0.05).mean(), (ok["nullB_p"] < 0.05).mean())


if __name__ == "__main__":
    main()
