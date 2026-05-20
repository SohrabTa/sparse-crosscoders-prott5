"""
Score crosscoder features against ProteinGym DMS assays.

For a curated set of (DMS_assay, paired_concept) targets:
  1. Load the wild-type sequence and its DMS variants (mutant, mutated_sequence,
     DMS_score) from ProteinGym.
  2. Identify the candidate crosscoder features paired to the concept, using the
     SAME thresholds the InterPLM dashboard applies (the input pairings file
     is already filtered to f1_per_domain > 0.5 ∩ valid/test; we sort by
     f1_per_domain desc to match the dashboard's "Top Feature" + "Other
     Candidates" ordering).
  3. Run ProtT5 once on the wild-type sequence to get baseline crosscoder
     activations at every residue, for the candidate feature subset only.
  4. Run ProtT5 on each variant sequence (batched), encode through the
     crosscoder, and extract activations at the mutated position(s).
  5. Compute per-variant Δactivation = act_mutant[pos] - act_wt[pos] for each
     candidate feature.
  6. Aggregate: Spearman correlation between |Δactivation| and the variant's
     |DMS_score - median(DMS_score)| (magnitude-of-effect hypothesis), and
     between signed Δactivation and DMS_score directly.

Outputs (per assay):
  * <output_dir>/<DMS_id>__<concept>__per_variant.parquet
      long-format: mutant, position(s), DMS_score, feature, act_wt, act_mut,
      delta_act
  * <output_dir>/summary.csv
      one row per (DMS_id, concept, feature): n_variants, spearman_abs,
      spearman_signed, pearson_abs, plus the feature's training f1_per_domain.

Designed to be run on a single GPU (H100 takes ~3 min/assay average for the
recommended 10-assay set).

Usage (from repo root):
    uv run python repos/sparse-crosscoders-prott5/score_proteingym_features.py \\
        --crosscoder_dir model_checkpoints/.../crashed_epoch_0_step_2519836 \\
        --pairings data/crosscoder_eval/uniprotkb_modern_score45_67k/test_counts/heldout_all_top_pairings.csv \\
        --matches data/proteingym_concept_matches.csv \\
        --dms_dir data/DMS_ProteinGym_substitutions \\
        --output_dir data/proteingym_scoring \\
        --batch_size 4
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
from scipy.stats import ConstantInputWarning, pearsonr, spearmanr

# Many candidate features stay silent on most proteins (zero Δactivation
# across all variants), which makes scipy emit ConstantInputWarning per
# correlation call. The resulting NaN is the desired signal; the warning
# itself just floods stderr.
warnings.filterwarnings("ignore", category=ConstantInputWarning)

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "repos" / "InterPLM"))
sys.path.insert(0, str(_REPO_ROOT / "repos" / "crosscode"))

from interplm.embedders.prott5 import ProtT5CrosscoderEmbedder  # noqa: E402
from interplm.sae.inference import load_sae  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("score_pg")


# Curated 10-assay subset useful for narrow / quick experiments. The script
# defaults to scoring ALL (assay, concept) rows from --matches; to restrict to
# just this set, pass --assays pointing at a file containing these rows.
FOCUSED_ASSAYS = [
    # Kinase set (4): same feature 4455 / 1322 across four independent assays
    ("SRC_HUMAN_Ahler_2019", "Domain_Protein kinase"),
    ("SRC_HUMAN_Chakraborty_2023_binding-DAS_25uM", "Domain_Protein kinase"),
    ("SRC_HUMAN_Nguyen_2022", "Domain_Protein kinase"),
    ("MK01_HUMAN_Brenan_2016", "Domain_Protein kinase"),
    # Single-domain functional set (3)
    ("PPARG_HUMAN_Majithia_2016", "Domain_NR LBD"),
    ("NUD15_HUMAN_Suiter_2020", "Domain_Nudix hydrolase"),
    ("AACC1_PSEAI_Dandage_2018", "Domain_N-acetyltransferase"),
    # Disulfide-bond control: same protein, same feature, two selection types
    ("VKOR1_HUMAN_Chiasson_2020_abundance", "Disulfide bond"),
    ("VKOR1_HUMAN_Chiasson_2020_activity", "Disulfide bond"),
    # Region_Disordered baseline (1)
    ("NKX31_HUMAN_Tsuboyama_2023_2L9R", "Region_Disordered"),
]


MUTANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z\*])$")


def parse_mutant_positions(mutant_field: str) -> List[Tuple[str, int, str]]:
    out: List[Tuple[str, int, str]] = []
    for token in mutant_field.split(":"):
        m = MUTANT_RE.match(token.strip())
        if not m:
            return []
        out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


def load_candidate_features(
    pairings_path: Path, concept: str
) -> List[Tuple[int, float]]:
    """Return [(feature, f1_per_domain), ...] for `concept`, sorted as the
    dashboard does (f1_per_domain, recall_per_domain, tp) descending. Includes
    all candidates already present in the file (the file is pre-filtered to
    f1_per_domain > 0.5 on both valid and test splits by InterPLM's pipeline)."""
    df = pd.read_csv(pairings_path)
    sub = df[df["concept"] == concept].copy()
    if sub.empty:
        return []
    # Apply dashboard's tp filter (no-op for this file but kept for safety).
    sub = sub.query("tp_per_domain >= 2 or tp >= 2")
    sub = sub.sort_values(
        ["f1_per_domain", "recall_per_domain", "tp"], ascending=False
    ).drop_duplicates("feature", keep="first")
    return list(
        zip(
            sub["feature"].astype(int).tolist(),
            sub["f1_per_domain"].astype(float).tolist(),
        )
    )


@torch.no_grad()
def encode_sequence(
    embedder: ProtT5CrosscoderEmbedder,
    crosscoder,
    sequence: str,
    feat_idx: torch.Tensor,
    device: torch.device,
) -> np.ndarray:
    """Embed one sequence through ProtT5 (24 hookpoints) and the crosscoder,
    returning per-residue activations for the selected feature subset.

    Returns: np.ndarray of shape (seq_len, n_feats)
    """
    # extract_embeddings returns (seq_len, M=1, P=24, D=1024) for batch_size=1
    emb = embedder.extract_embeddings([sequence], batch_size=1)
    emb = emb.to(device=device, dtype=next(crosscoder.parameters()).dtype)
    # Crosscoder treats each amino-acid token as a batch element.
    latents = crosscoder.encode_feat_subset(emb, feat_idx.tolist())
    return latents.float().cpu().numpy()


@torch.no_grad()
def encode_sequences_batched(
    embedder: ProtT5CrosscoderEmbedder,
    crosscoder,
    sequences: List[str],
    feat_idx: torch.Tensor,
    device: torch.device,
    batch_size: int,
) -> List[np.ndarray]:
    """Embed a list of sequences (each potentially different lengths) and return
    a list of per-residue activations (seq_len_i, n_feats). Uses
    extract_embeddings_with_boundaries to recover per-sequence slices."""
    bundle = embedder.extract_embeddings_with_boundaries(
        sequences, layer=-1, batch_size=batch_size
    )
    emb_all = bundle["embeddings"].to(
        device=device, dtype=next(crosscoder.parameters()).dtype
    )
    boundaries: List[Tuple[int, int]] = bundle["boundaries"]
    latents_all = crosscoder.encode_feat_subset(emb_all, feat_idx.tolist())
    latents_all = latents_all.float().cpu().numpy()
    return [latents_all[s:e] for s, e in boundaries]


def _per_feature_summary(
    per_variant: pd.DataFrame,
    dms_id: str,
    concept: str,
    feat_f1pd: Dict[int, float],
) -> List[dict]:
    """Compute Spearman/Pearson correlations for each feature in `per_variant`."""
    summary_rows = []
    median_score = per_variant["DMS_score"].median()
    for feat_id, grp in per_variant.groupby("feature"):
        if len(grp) < 5:
            continue
        delta = grp["delta_act"].values
        score = grp["DMS_score"].values
        score_abs = np.abs(score - median_score)
        delta_abs = np.abs(delta)
        try:
            rho_abs = spearmanr(delta_abs, score_abs, nan_policy="omit").statistic
        except Exception:
            rho_abs = float("nan")
        try:
            rho_signed = spearmanr(delta, score, nan_policy="omit").statistic
        except Exception:
            rho_signed = float("nan")
        try:
            r_abs = pearsonr(delta_abs, score_abs)[0]
        except Exception:
            r_abs = float("nan")
        summary_rows.append(
            {
                "DMS_id": dms_id,
                "concept": concept,
                "feature": int(feat_id),
                "n_variants": int(len(grp)),
                "f1_per_domain": feat_f1pd[int(feat_id)],
                "spearman_abs": float(rho_abs) if rho_abs is not None else float("nan"),
                "spearman_signed": float(rho_signed)
                if rho_signed is not None
                else float("nan"),
                "pearson_abs": float(r_abs),
                "mean_abs_delta": float(np.mean(np.abs(delta))),
                "max_abs_delta": float(np.max(np.abs(delta))),
            }
        )
    return summary_rows


def score_one_assay_multi_concept(
    dms_id: str,
    concepts: List[str],
    dms_dir: Path,
    pairings_path: Path,
    embedder: ProtT5CrosscoderEmbedder,
    crosscoder,
    device: torch.device,
    output_dir: Path,
    batch_size: int,
    max_variants: Optional[int] = None,
) -> List[dict]:
    """Score all candidate features for all concepts paired to this assay in
    a single pass. ProtT5 + crosscoder run *once* over the WT and over each
    variant; the result is then sliced per concept's feature subset. Saves
    ~(n_concepts × variants) duplicate forward passes vs the per-concept loop.

    Returns: combined summary rows across all concepts.
    """
    # Resolve candidate features per concept + union for the single forward pass
    feats_per_concept: Dict[str, List[int]] = {}
    f1pd_per_concept: Dict[str, Dict[int, float]] = {}
    for concept in concepts:
        candidates = load_candidate_features(pairings_path, concept)
        if not candidates:
            log.warning("  no candidate features for %s; skipping", concept)
            continue
        feats_per_concept[concept] = [f for f, _ in candidates]
        f1pd_per_concept[concept] = {f: f1 for f, f1 in candidates}
    if not feats_per_concept:
        return []

    union_feats = sorted({f for fs in feats_per_concept.values() for f in fs})
    union_feat_pos = {f: i for i, f in enumerate(union_feats)}
    union_idx = torch.tensor(union_feats, dtype=torch.long)
    log.info(
        "  %d concept(s); %d unique candidate features across union",
        len(feats_per_concept),
        len(union_feats),
    )

    dms_path = dms_dir / f"{dms_id}.csv"
    if not dms_path.exists():
        log.error("  DMS file missing: %s", dms_path)
        return []
    dms = pd.read_csv(dms_path)
    if max_variants and len(dms) > max_variants:
        dms = dms.sample(max_variants, random_state=42).reset_index(drop=True)
    log.info("  %d variants to score", len(dms))

    # WT seq: take the first variant's mutated_sequence and revert its mutation(s)
    first_mut = parse_mutant_positions(str(dms["mutant"].iloc[0]))
    if not first_mut:
        log.error("  could not parse mutant '%s'", dms["mutant"].iloc[0])
        return []
    wt_seq = list(dms["mutated_sequence"].iloc[0])
    for wt_aa, pos1, _ in first_mut:
        wt_seq[pos1 - 1] = wt_aa
    wt_seq_str = "".join(wt_seq)

    log.info("  computing WT activations (len=%d)", len(wt_seq_str))
    wt_acts = encode_sequence(embedder, crosscoder, wt_seq_str, union_idx, device)

    # Per-variant activations at the mutated position(s), stored once for the
    # union of features. We accumulate into a flat list and slice per concept.
    flat_rows = []  # (mutant, position, DMS_score, act_wt_vec, act_mut_vec)
    n_batches = (len(dms) + batch_size - 1) // batch_size
    for batch_idx in range(n_batches):
        b_start = batch_idx * batch_size
        b_end = min(b_start + batch_size, len(dms))
        batch_df = dms.iloc[b_start:b_end]
        seqs = batch_df["mutated_sequence"].tolist()
        mut_acts_list = encode_sequences_batched(
            embedder, crosscoder, seqs, union_idx, device, batch_size
        )
        for (_, vrow), mut_acts in zip(batch_df.iterrows(), mut_acts_list):
            parsed = parse_mutant_positions(str(vrow["mutant"]))
            if not parsed:
                continue
            for _, pos1, _ in parsed:
                if pos1 - 1 >= mut_acts.shape[0] or pos1 - 1 >= wt_acts.shape[0]:
                    continue
                flat_rows.append(
                    (
                        str(vrow["mutant"]),
                        pos1,
                        float(vrow["DMS_score"]),
                        wt_acts[pos1 - 1],  # (n_union_feats,)
                        mut_acts[pos1 - 1],  # (n_union_feats,)
                    )
                )
        if batch_idx % 100 == 0:
            log.info("    batch %d/%d", batch_idx + 1, n_batches)

    if not flat_rows:
        log.warning("  no scorable variants for %s", dms_id)
        return []

    # Slice per concept and write per-variant parquet + per-feature summary
    summary_rows = []
    for concept, feat_ids in feats_per_concept.items():
        local_pos = [union_feat_pos[f] for f in feat_ids]
        rows = []
        for mutant, pos1, dms_score, wt_vec, mut_vec in flat_rows:
            wt_sel = wt_vec[local_pos]
            mut_sel = mut_vec[local_pos]
            delta = mut_sel - wt_sel
            for fi, feat_id in enumerate(feat_ids):
                rows.append(
                    {
                        "mutant": mutant,
                        "position": pos1,
                        "DMS_score": dms_score,
                        "feature": feat_id,
                        "act_wt": float(wt_sel[fi]),
                        "act_mut": float(mut_sel[fi]),
                        "delta_act": float(delta[fi]),
                    }
                )
        per_variant = pd.DataFrame(rows)
        safe_concept = concept.replace("/", "_").replace(" ", "_")
        out_path = output_dir / f"{dms_id}__{safe_concept}__per_variant.parquet"
        per_variant.to_parquet(out_path, index=False)
        log.info(
            "  %s -> wrote %d (variant, feature) rows",
            concept,
            len(per_variant),
        )
        summary_rows.extend(
            _per_feature_summary(
                per_variant, dms_id, concept, f1pd_per_concept[concept]
            )
        )

    return summary_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--crosscoder_dir",
        type=Path,
        required=True,
        help="Directory with ae.pt / ae_normalized.pt + config.yaml",
    )
    ap.add_argument(
        "--checkpoint",
        type=str,
        default="ae_normalized.pt",
        help="Filename inside crosscoder_dir (default: ae_normalized.pt)",
    )
    ap.add_argument("--pairings", type=Path, required=True)
    ap.add_argument("--matches", type=Path, required=True)
    ap.add_argument("--dms_dir", type=Path, required=True)
    ap.add_argument("--output_dir", type=Path, required=True)
    ap.add_argument(
        "--assays",
        type=Path,
        default=None,
        help="Optional TSV/CSV with columns DMS_id,concept. If omitted, "
        "every (DMS_id, matched_concept) pair from --matches is scored.",
    )
    ap.add_argument("--batch_size", type=int, default=4)
    ap.add_argument(
        "--max_variants", type=int, default=None, help="Cap variants per assay (debug)"
    )
    ap.add_argument(
        "--max_seq_len",
        type=int,
        default=None,
        help="Skip assays whose target_seq_len exceeds this. Recommended: 2048 "
        "(crosscoder was trained on ≤512; signal beyond that is dubious, "
        "and the few mega-assays dominate runtime without contributing usable data).",
    )
    ap.add_argument(
        "--device",
        type=str,
        default=None,
        help="cuda / cpu / mps; auto-detect if omitted",
    )
    ap.add_argument(
        "--skip_existing",
        action="store_true",
        help="Skip assays whose per_variant.parquet files all already exist",
    )
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    matches = pd.read_csv(args.matches)

    if args.assays:
        assay_df = pd.read_csv(
            args.assays, sep="\t" if args.assays.suffix == ".tsv" else ","
        )
        wanted = set(zip(assay_df["DMS_id"], assay_df["concept"]))
        work = matches[
            matches.apply(
                lambda r: (r["DMS_id"], r["matched_concept"]) in wanted, axis=1
            )
        ].copy()
    else:
        work = matches.copy()

    if args.max_seq_len is not None:
        before = work["DMS_id"].nunique()
        excluded = sorted(work[work["target_seq_len"] > args.max_seq_len]["DMS_id"].unique())
        work = work[work["target_seq_len"] <= args.max_seq_len].copy()
        if excluded:
            log.info(
                "Filtered %d/%d assays by max_seq_len=%d: dropped %s",
                len(excluded), before, args.max_seq_len, excluded,
            )

    # Group by assay so we do one ProtT5 + crosscoder pass per assay, slicing
    # the cached activations for each matched concept.
    assays_grouped: Dict[str, List[str]] = defaultdict(list)
    for _, r in work.iterrows():
        assays_grouped[r["DMS_id"]].append(r["matched_concept"])

    log.info(
        "Scoring %d assays / %d (assay, concept) pairs",
        len(assays_grouped),
        sum(len(v) for v in assays_grouped.values()),
    )

    if args.device:
        device = torch.device(args.device)
    else:
        if torch.cuda.is_available():
            device = torch.device("cuda")
        elif torch.backends.mps.is_available():
            device = torch.device("mps")
        else:
            device = torch.device("cpu")
    log.info("Device: %s", device)

    log.info("Loading ProtT5...")
    embedder = ProtT5CrosscoderEmbedder(device=str(device))
    log.info("Loading crosscoder from %s/%s", args.crosscoder_dir, args.checkpoint)
    crosscoder = load_sae(
        args.crosscoder_dir, device=str(device), model_name=args.checkpoint
    )
    crosscoder.eval()

    all_summary: List[dict] = []
    for i, (dms_id, concepts) in enumerate(assays_grouped.items()):
        log.info(
            "=== [%d/%d] %s | concepts=%s ===",
            i + 1,
            len(assays_grouped),
            dms_id,
            concepts,
        )
        if args.skip_existing:
            expected = [
                args.output_dir
                / f"{dms_id}__{c.replace('/', '_').replace(' ', '_')}__per_variant.parquet"
                for c in concepts
            ]
            if all(p.exists() for p in expected):
                log.info("  all per_variant files exist; skipping")
                continue
        summary_rows = score_one_assay_multi_concept(
            dms_id=dms_id,
            concepts=concepts,
            dms_dir=args.dms_dir,
            pairings_path=args.pairings,
            embedder=embedder,
            crosscoder=crosscoder,
            device=device,
            output_dir=args.output_dir,
            batch_size=args.batch_size,
            max_variants=args.max_variants,
        )
        all_summary.extend(summary_rows)
        # Persist summary incrementally so a crashed job still leaves results.
        if all_summary:
            pd.DataFrame(all_summary).to_csv(
                args.output_dir / "summary.csv", index=False
            )

    if not all_summary:
        log.warning("No summary rows produced")
        return

    summary = pd.DataFrame(all_summary)
    # For each (assay, concept), surface the best feature by spearman_abs
    summary = summary.sort_values(
        ["DMS_id", "concept", "spearman_abs"], ascending=[True, True, False]
    )
    summary_path = args.output_dir / "summary.csv"
    summary.to_csv(summary_path, index=False)
    log.info("Wrote %d summary rows -> %s", len(summary), summary_path)

    # Highlights: best feature per (assay, concept)
    best = summary.drop_duplicates(["DMS_id", "concept"], keep="first")
    log.info("Per-assay best feature by |Δact| ~ |ΔDMS| Spearman:")
    for _, r in best.iterrows():
        log.info(
            "  %-50s %-30s feat=%-5d f1pd=%.2f  ρ_abs=%+.3f  ρ_signed=%+.3f  n=%d",
            r["DMS_id"],
            r["concept"],
            r["feature"],
            r["f1_per_domain"],
            r["spearman_abs"],
            r["spearman_signed"],
            r["n_variants"],
        )


if __name__ == "__main__":
    sys.exit(main())
