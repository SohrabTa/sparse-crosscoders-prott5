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

CONFOUND CONTROLS (the point/site concepts have a sharp one — amino-acid identity):
  C) AA-identity-restricted AUC: restrict both groups to positions whose WT residue is one the
     feature fires on (disulfide -> cysteines), then ask whether the FIRING residues have higher
     |effect| than non-firing residues OF THE SAME TYPE. This defuses "it just detects cysteines"
     — the identity prior that survives even a random-init model (InterPLM). No extra data needed.
     Degenerate (auc_id=nan) when the feature fires on ~all residues of its type — itself the
     finding that the label is identity-level, not function-level.
  D) Conservation-residualized AUC + conservation-stratified mask null: with --conservation
     <CSV: DMS_id,pos,score> (per-position MSA Shannon entropy from harvest_msa_conservation.py),
     rank-residualize |effect| on conservation and recompute AUC; the stratified null matches the
     firing set's conservation distribution. Asks for fitness importance BEYOND conservation.
     Optional — the MSA files are not local (ProteinGym download); the hook is ready for them.

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


def auc_from_ranks(rank_vals: np.ndarray, mask: np.ndarray) -> float:
    """AUC = P(val[mask] > val[~mask]) from precomputed ranks over the full set."""
    n = len(rank_vals)
    k = int(mask.sum())
    if k < 2 or k > n - 2:
        return np.nan
    R = float(rank_vals[mask].sum())
    return (R - k * (k + 1) / 2.0) / (k * (n - k))


def mask_perm_p(rank_vals: np.ndarray, k: int, auc_obs: float,
                rng: np.random.Generator, n_perm: int,
                strata: np.ndarray = None) -> tuple[float, float]:
    """Mask-permutation null: draw k positions (optionally stratified to match a covariate's
    distribution over the firing set) and recompute AUC. Returns (p95, empirical_p)."""
    n = len(rank_vals)
    if k < 2 or k > n - 2:
        return (np.nan, np.nan)
    aucs = []
    if strata is None:
        for _ in range(n_perm):
            msk = np.zeros(n, bool); msk[rng.choice(n, size=k, replace=False)] = True
            aucs.append(auc_from_ranks(rank_vals, msk))
    else:
        # match the firing set's per-stratum counts (set externally via strata == firing strata)
        bin_ids, bin_counts = strata
        for _ in range(n_perm):
            sel = []
            for b, c in bin_counts.items():
                pool = np.where(bin_ids == b)[0]
                sel.extend(rng.choice(pool, size=min(c, len(pool)), replace=False))
            msk = np.zeros(n, bool); msk[np.asarray(sel, int)] = True
            aucs.append(auc_from_ranks(rank_vals, msk))
    aucs = np.array([a for a in aucs if np.isfinite(a)])
    if not len(aucs):
        return (np.nan, np.nan)
    return (float(np.percentile(aucs, 95)), float((aucs >= auc_obs).mean()))


def restricted_auc(dev: np.ndarray, wt_aa: np.ndarray, fire: np.ndarray) -> tuple[float, int, int]:
    """AA-identity-restricted AUC: among positions whose WT residue is one the feature fires on
    (e.g. only cysteines for a disulfide feature), do the FIRING positions have higher |effect|
    than the non-firing ones? Defuses 'it just detects the residue type' — the identity prior that
    survives even in a random-init model. Returns (auc, n_subset, k_in_subset)."""
    aa_fire = set(wt_aa[fire].tolist())
    sub = np.isin(wt_aa, list(aa_fire))
    n_sub = int(sub.sum())
    rd = rankdata(dev[sub])
    fmask = fire[sub]
    return (auc_from_ranks(rd, fmask), n_sub, int(fmask.sum()))


def residualize(dev: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """Rank-residualize |effect| on a covariate (conservation): rank(dev) minus its linear fit on
    rank(cov). AUC on the residual asks for fitness importance BEYOND what conservation explains."""
    rd, rc = rankdata(dev), rankdata(cov)
    a, b = np.polyfit(rc, rd, 1)
    return rd - (a * rc + b)


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
    ap.add_argument("--conservation", type=Path, default=None,
                    help="optional CSV (DMS_id,pos,score) of per-position conservation (e.g. MSA "
                         "Shannon entropy from harvest_msa_conservation.py); higher score = more "
                         "variable. Enables the conservation-residualized AUC + stratified null.")
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

    cons_map = {}
    if args.conservation:
        cons_df = pd.read_csv(args.conservation)
        for did, g in cons_df.groupby("DMS_id"):
            cons_map[did] = dict(zip(g["pos"].astype(int), g["score"].astype(float)))
        log.info("conservation loaded for %d assays", len(cons_map))

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
        wt_aa = np.array(list(wt_str))[pos0]           # WT residue at each scored position

        # optional per-position conservation aligned to the scored positions (1-indexed)
        cons = None
        if dms_id in cons_map:
            cvals = np.array([cons_map[dms_id].get(int(p), np.nan) for p in good.values])
            if np.isfinite(cvals).sum() >= MIN_POS:
                cons = cvals
        rank_resid = None
        cons_strata = None
        if cons is not None:
            ok_c = np.isfinite(cons)
            # residualized |effect| ranks (NaN conservation positions kept at rank-residual 0)
            rr = np.zeros(n_pos)
            rr[ok_c] = rankdata(residualize(dev[ok_c], cons[ok_c]))
            rank_resid = rr
            bins = pd.qcut(rankdata(np.where(ok_c, cons, np.nanmedian(cons))),
                           q=min(5, max(2, ok_c.sum() // 5)), labels=False, duplicates="drop")
            cons_strata = bins

        wt_act = encode_wt(crosscoder, embedder.extract_embeddings([wt_str], batch_size=1), device)
        A = wt_act[pos0, :]                             # (n_pos, 8192) at scored positions
        fire_counts = (A > 0).sum(axis=0)               # k_g for every feature on these positions

        # density-matched live pool, recomputed per labelled feature's k below
        live_k = fire_counts[live_idx]

        log.info("  %s: L=%d  n_pos=%d  labelled=%d  cons=%s", dms_id, L, n_pos, len(grp),
                 "yes" if cons is not None else "no")
        for _, row in grp.iterrows():
            f = int(row["feature"])
            fire = A[:, f] > 0
            k = int(fire.sum())
            auc_obs = auc_from_firing(rank_dev, fire)
            rec = {"DMS_id": dms_id, "concept": row["concept"], "feature": f,
                   "n_pos": n_pos, "k_fire": k, "auc_obs": auc_obs}
            if not np.isfinite(auc_obs):
                rec.update(status="degenerate_k")
                results.append(rec)
                if args.smoke: log.info("    f%d concept=%s k=%d  AUC=nan (degenerate)", f, row["concept"], k)
                continue

            # Null A: density-matched random live features (nearest-k, exclude the labelled feature)
            order = np.argsort(np.abs(live_k - k))
            pool = [int(live_idx[j]) for j in order if int(live_idx[j]) != f][:N_NULL_FEATURES]
            aucA = np.array([auc_from_firing(rank_dev, A[:, g] > 0) for g in pool])
            aucA = aucA[np.isfinite(aucA)]
            nullA_mean = float(np.mean(aucA)) if len(aucA) else np.nan
            nullA_p = float((aucA >= auc_obs).mean()) if len(aucA) else np.nan
            nullA_z = float((auc_obs - nullA_mean) / (np.std(aucA) + 1e-9)) if len(aucA) > 1 else np.nan
            realized_k = float(np.median(fire_counts[pool])) if pool else np.nan

            # Null B: within-sequence mask permutation at fixed k
            nullB_p95, nullB_p = mask_perm_p(rank_dev, k, auc_obs, rng, N_PERM)

            # AA-identity-restricted AUC (compares firing vs non-firing positions of the SAME WT
            # residue type) + its own mask-permutation p within that subset
            auc_id, n_id, k_id = restricted_auc(dev, wt_aa, fire)
            aa_fire = set(wt_aa[fire].tolist())
            sub = np.isin(wt_aa, list(aa_fire))
            nullID_p = np.nan
            if np.isfinite(auc_id):
                _, nullID_p = mask_perm_p(rankdata(dev[sub]), k_id, auc_id, rng, N_PERM)

            # Conservation-residualized AUC + conservation-stratified mask null (only if --conservation)
            auc_cr, nullCons_p = np.nan, np.nan
            if rank_resid is not None:
                auc_cr = auc_from_ranks(rank_resid, fire)
                if np.isfinite(auc_cr) and cons_strata is not None:
                    bc = pd.Series(cons_strata[fire]).value_counts().to_dict()
                    _, nullCons_p = mask_perm_p(rank_resid, k, auc_cr, rng, N_PERM,
                                                strata=(cons_strata, bc))

            rec.update(nullA_n=len(aucA), nullA_mean=nullA_mean, nullA_p=nullA_p, nullA_z=nullA_z,
                       nullA_median_k=realized_k, nullB_p95=nullB_p95, nullB_p=nullB_p,
                       auc_id=auc_id, n_id=n_id, k_id=k_id, nullID_p=nullID_p,
                       auc_consresid=auc_cr, nullCons_p=nullCons_p, status="ok")
            results.append(rec)
            if args.smoke:
                log.info("    f%d %s k=%d/%d AUC=%.3f | nullA m=%.3f p=%.3f z=%.2f | nullB p=%.3f | "
                         "ID-restr AUC=%.3f (n=%d,k=%d) p=%.3f | consResid AUC=%s p=%s",
                         f, row["concept"], k, n_pos, auc_obs, nullA_mean, nullA_p, nullA_z, nullB_p,
                         auc_id, n_id, k_id, nullID_p,
                         f"{auc_cr:.3f}" if np.isfinite(auc_cr) else "na",
                         f"{nullCons_p:.3f}" if np.isfinite(nullCons_p) else "na")

    out = pd.DataFrame(results)
    out_path = args.output_dir / "design1_per_feature.csv"
    out.to_csv(out_path, index=False)
    log.info("wrote %s (%d rows)", out_path, len(out))
    ok = out[out["status"] == "ok"] if "status" in out else out
    if len(ok):
        idok = ok[ok["auc_id"].notna()]
        log.info("\n=== summary (status=ok, %d rows) ===", len(ok))
        log.info("median AUC_obs=%.3f  median nullA_mean=%.3f  frac nullA_p<0.05=%.2f  frac nullB_p<0.05=%.2f",
                 ok["auc_obs"].median(), ok["nullA_mean"].median(),
                 (ok["nullA_p"] < 0.05).mean(), (ok["nullB_p"] < 0.05).mean())
        if len(idok):
            log.info("ID-restricted: %d rows  median AUC_id=%.3f  frac nullID_p<0.05=%.2f",
                     len(idok), idok["auc_id"].median(), (idok["nullID_p"] < 0.05).mean())
        if ok["auc_consresid"].notna().any():
            cr = ok[ok["auc_consresid"].notna()]
            log.info("cons-residualized: %d rows  median AUC=%.3f  frac nullCons_p<0.05=%.2f",
                     len(cr), cr["auc_consresid"].median(), (cr["nullCons_p"] < 0.05).mean())


if __name__ == "__main__":
    main()
