"""
Match ProteinGym DMS assays to crosscoder feature-concept pairings.

For each DMS assay in ProteinGym:
  1. Resolve its WT UniProt entry and fetch annotations using the same UniProt
     fields as `download_uniprot_dataset.py`.
  2. Parse those annotations with the SAME functions InterPLM's eval pipeline
     uses (`process_categorical_feature`, `process_binary_feature`,
     `process_interaction_feature`) — to avoid implementation drift.
  3. Restrict the concept vocabulary to the 51 paired concepts from
     `heldout_all_top_pairings.csv` and intersect with what each protein has.
  4. Count how many DMS variants have their mutated position inside the residue
     range of each matched paired concept.

Output: one row per (DMS_assay, matched_concept), sorted by variants-in-concept.
This identifies the concept-matched subset of assays where the interpretability
hypothesis ("feature F's activation delta tracks fitness delta when F is paired
to concept C and the mutation lands in C's residues") can be tested directly.

Usage (from repo root):
    uv run --with pandas --with requests \
      python repos/sparse-crosscoders-prott5/match_proteingym_concepts.py \
        --pairings data/crosscoder_eval/pre-auxfix/real/uniprotkb_modern_score45_67k/test_counts/heldout_all_top_pairings.csv \
        --dms_reference /tmp/DMS_substitutions.csv \
        --dms_dir data/external/DMS_ProteinGym_substitutions \
        --cache_dir data/external/proteingym_uniprot_cache \
        --output data/proteingym_concept_matches.csv
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
import requests

# Make InterPLM importable without installation
_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "repos" / "InterPLM"))

from interplm.analysis.concepts.parsing_utils import (  # noqa: E402
    process_binary_feature,
    process_categorical_feature,
    process_interaction_feature,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("match_pg")


# UniProt fields — verbatim from `download_uniprot_dataset.py` so the TSV columns
# we receive are identical to what InterPLM's eval pipeline consumed.
UNIPROT_FIELDS = [
    "accession", "reviewed", "protein_name", "length", "sequence", "ec",
    "ft_act_site", "ft_binding", "cc_cofactor", "ft_disulfid", "ft_carbohyd",
    "ft_lipid", "ft_mod_res", "ft_signal", "ft_transit", "ft_helix", "ft_turn",
    "ft_strand", "ft_coiled", "cc_domain", "ft_compbias", "ft_domain", "ft_motif",
    "ft_region", "ft_zn_fing", "xref_alphafolddb",
]

# Parser dispatch: InterPLM-column-name -> (parser_kind, short_tag).
# `short_tag` is the leading word in UniProt FT values (e.g. "DOMAIN 35..96; ...").
# These match how `extract_annotations.py` derives `col_name` from the data.
COLUMN_INFO: Dict[str, Tuple[str, str]] = {
    # categorical (have sub-categories via /note= or /ligand=)
    "Active site": ("categorical", "ACT_SITE"),
    "Binding site": ("categorical", "BINDING"),
    "Cofactor": ("categorical", "COFACTOR"),
    "Glycosylation": ("categorical", "CARBOHYD"),
    "Modified residue": ("categorical", "MOD_RES"),
    "Transit peptide": ("categorical", "TRANSIT"),
    "Compositional bias": ("categorical", "COMPBIAS"),
    "Domain [FT]": ("categorical", "DOMAIN"),
    "Region": ("categorical", "REGION"),
    "Zinc finger": ("categorical", "ZN_FING"),
    "Motif": ("categorical", "MOTIF"),
    "Signal peptide": ("categorical", "SIGNAL"),
    # binary
    "Helix": ("binary", "HELIX"),
    "Turn": ("binary", "TURN"),
    "Beta strand": ("binary", "STRAND"),
    "Coiled coil": ("binary", "COILED"),
    "Lipidation": ("binary", "LIPID"),
    # paired interactions
    "Disulfide bond": ("paired", "DISULFID"),
}

# Concept-name prefix -> InterPLM column name. InterPLM strips " [FT]" from
# "Domain [FT]" when assembling concept names, so "Domain_Protein kinase" maps
# back to the "Domain [FT]" UniProt column.
CONCEPT_PREFIX_TO_COLUMN: Dict[str, str] = {
    "Domain": "Domain [FT]",
    "Modified residue": "Modified residue",
    "Region": "Region",
    "Motif": "Motif",
    "Zinc finger": "Zinc finger",
    "Compositional bias": "Compositional bias",
    "Glycosylation": "Glycosylation",
    "Active site": "Active site",
    "Binding site": "Binding site",
    "Cofactor": "Cofactor",
    "Transit peptide": "Transit peptide",
    "Signal peptide": "Signal peptide",
}

BARE_COLUMN_CONCEPTS: Set[str] = {
    name for name, (kind, _) in COLUMN_INFO.items() if kind in ("binary", "paired")
}


MUTANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z\*])$")


def parse_concept_name(concept_name: str) -> Tuple[str, Optional[str]]:
    """Concept name (e.g. 'Domain_Protein kinase') -> (InterPLM column, subcategory_or_None)."""
    if concept_name in BARE_COLUMN_CONCEPTS:
        return concept_name, None
    for prefix, col in CONCEPT_PREFIX_TO_COLUMN.items():
        if concept_name.startswith(f"{prefix}_"):
            return col, concept_name[len(prefix) + 1:]
    raise ValueError(f"Unrecognized concept name: {concept_name!r}")


def concept_name_for(col: str, sub: Optional[str]) -> str:
    """Inverse of parse_concept_name — produce the InterPLM-style concept name."""
    if sub is None:
        return col
    prefix = col.replace(" [FT]", "")
    return f"{prefix}_{sub}"


def load_pairings(path: Path) -> Tuple[Dict[str, List[dict]], Dict[str, List[str]], List[str], List[str]]:
    """Returns:
      pairings: concept_name -> [{feature, f1, precision, recall}, ...]  (sorted by F1 desc)
      categorical_options: InterPLM column -> [subcategory, ...]
      binary_cols, paired_cols: list of InterPLM columns to process
    """
    df = pd.read_csv(path)
    pairings: Dict[str, List[dict]] = {}
    for _, r in df.iterrows():
        pairings.setdefault(r["concept"], []).append({
            "feature": int(r["feature"]),
            "f1": float(r["f1"]),
            "f1_per_domain": float(r["f1_per_domain"]),
            "precision": float(r["precision"]),
            "recall": float(r["recall"]),
            "recall_per_domain": float(r["recall_per_domain"]),
            "tp": float(r["tp"]),
            "tp_per_domain": float(r["tp_per_domain"]),
        })
    # Match the InterPLM dashboard sort: f1_per_domain, recall_per_domain, tp (all desc).
    # This makes the top candidate identical to what the dashboard shows.
    for c in pairings:
        pairings[c].sort(
            key=lambda d: (-d["f1_per_domain"], -d["recall_per_domain"], -d["tp"])
        )

    categorical_options: Dict[str, List[str]] = defaultdict(list)
    binary_cols: List[str] = []
    paired_cols: List[str] = []
    for concept in pairings:
        col, sub = parse_concept_name(concept)
        kind, _ = COLUMN_INFO[col]
        if kind == "categorical":
            if sub not in categorical_options[col]:
                categorical_options[col].append(sub)
        elif kind == "binary":
            if col not in binary_cols:
                binary_cols.append(col)
        elif kind == "paired":
            if col not in paired_cols:
                paired_cols.append(col)
    # InterPLM's process_categorical_feature unconditionally writes to a
    # category_indices["any"] bucket; the training pipeline guarantees "any" is
    # always present in category_options. Mirror that here so behavior matches
    # the eval pipeline exactly (the resulting "<col>_any" concepts are not
    # paired to any feature and get filtered out naturally).
    for col in categorical_options:
        if "any" not in categorical_options[col]:
            categorical_options[col].append("any")
    return pairings, dict(categorical_options), binary_cols, paired_cols


def parse_mutant_positions(mutant_field: str) -> List[Tuple[str, int, str]]:
    """ProteinGym mutant cell -> [(wt, pos_1based, mut), ...]. [] on parse failure."""
    out: List[Tuple[str, int, str]] = []
    for token in mutant_field.split(":"):
        m = MUTANT_RE.match(token.strip())
        if not m:
            return []
        out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


def fetch_uniprot_row(uniprot_id: str, cache_dir: Path, sleep_s: float = 0.2) -> Optional[dict]:
    """Fetch a single UniProt entry by entry name (ProteinGym's `UniProt_ID` column
    holds entry names like 'A4_HUMAN', not accessions). Returns a row dict keyed
    by display column name. Cached as JSON. Returns {} for unmatched / synthetic
    proteins (e.g. ANCSZ), None on transient network failure."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{uniprot_id}.json"
    if cache_path.exists():
        try:
            return json.loads(cache_path.read_text())
        except json.JSONDecodeError:
            cache_path.unlink(missing_ok=True)

    url = (
        "https://rest.uniprot.org/uniprotkb/search"
        f"?query=id:{uniprot_id}&format=tsv&fields={','.join(UNIPROT_FIELDS)}"
    )
    r = None
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=30)
        except requests.RequestException as e:
            log.warning("UniProt fetch network error for %s (attempt %d): %s",
                        uniprot_id, attempt + 1, e)
            time.sleep(2 ** attempt)
            continue
        # Only retry on 5xx; 4xx is a permanent client error (bad query, etc).
        if r.status_code < 400:
            break
        if 400 <= r.status_code < 500:
            log.warning("UniProt %s: HTTP %d (not retrying)", uniprot_id, r.status_code)
            cache_path.write_text(json.dumps({}))
            return {}
        log.warning("UniProt %s: HTTP %d (attempt %d)", uniprot_id, r.status_code, attempt + 1)
        time.sleep(2 ** attempt)
    if r is None or r.status_code >= 400:
        return None

    lines = r.text.strip().splitlines()
    if len(lines) < 2:
        cache_path.write_text(json.dumps({}))
        time.sleep(sleep_s)
        return {}
    header = lines[0].split("\t")
    values = lines[1].split("\t")
    row = dict(zip(header, values))
    cache_path.write_text(json.dumps(row))
    time.sleep(sleep_s)
    return row


def concept_to_residue_set(
    uniprot_row: dict,
    seq_len: int,
    categorical_options: Dict[str, List[str]],
    binary_cols: List[str],
    paired_cols: List[str],
) -> Dict[str, Set[int]]:
    """Run InterPLM's parsers on the UniProt row and return
    concept_name -> set of 1-based residue positions belonging to that concept."""
    out: Dict[str, Set[int]] = {}

    for col, options in categorical_options.items():
        short_tag = COLUMN_INFO[col][1]
        col_data = uniprot_row.get(col, "")
        if not col_data or pd.isna(col_data):
            continue
        # process_categorical_feature returns a list of per-residue vectors (one per
        # subcategory option). Non-zero values are instance indices; treat >0 as "in".
        current_index: Dict[str, int] = {sub: 1 for sub in options}
        vectors, _ = process_categorical_feature(
            col_data, short_tag, options, seq_len, current_index
        )
        for sub, vec in zip(options, vectors):
            positions = {i + 1 for i, v in enumerate(vec) if v}
            if positions:
                out[concept_name_for(col, sub)] = positions

    for col in binary_cols:
        short_tag = COLUMN_INFO[col][1]
        col_data = uniprot_row.get(col, "")
        if not col_data or pd.isna(col_data):
            continue
        vec, _ = process_binary_feature(col_data, short_tag, seq_len, current_index=1)
        positions = {i + 1 for i, v in enumerate(vec) if v}
        if positions:
            out[concept_name_for(col, None)] = positions

    for col in paired_cols:
        short_tag = COLUMN_INFO[col][1]
        col_data = uniprot_row.get(col, "")
        if not col_data or pd.isna(col_data):
            continue
        indices, _pairs = process_interaction_feature(col_data, short_tag, seq_len)
        positions = {i + 1 for i, v in enumerate(indices) if v}
        if positions:
            out[concept_name_for(col, None)] = positions

    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairings", type=Path, required=True)
    ap.add_argument("--dms_reference", type=Path, required=True,
                    help="ProteinGym DMS_substitutions.csv reference file")
    ap.add_argument("--dms_dir", type=Path, required=True,
                    help="Directory containing the per-assay DMS_*.csv files")
    ap.add_argument("--cache_dir", type=Path, required=True,
                    help="Directory to cache UniProt API responses (JSON per accession)")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--limit", type=int, default=None, help="Process only first N assays")
    args = ap.parse_args()

    pairings, categorical_options, binary_cols, paired_cols = load_pairings(args.pairings)
    log.info(
        "Loaded %d paired concepts (%d feature pairings): "
        "%d categorical cols (%d subcategories), %d binary cols, %d paired cols",
        len(pairings), sum(len(v) for v in pairings.values()),
        len(categorical_options), sum(len(s) for s in categorical_options.values()),
        len(binary_cols), len(paired_cols),
    )

    ref = pd.read_csv(args.dms_reference)
    log.info("Loaded %d DMS assays", len(ref))
    if args.limit:
        ref = ref.head(args.limit)

    rows = []
    for i, r in ref.reset_index(drop=True).iterrows():
        dms_id = r["DMS_id"]
        dms_filename = r["DMS_filename"]
        full_id = r["UniProt_ID"]
        target_seq = r["target_seq"]
        selection = r.get("coarse_selection_type", "")
        fn_class = r.get("selection_type", "")

        log.info("[%d/%d] %s (id %s, len=%d, %s)",
                 i + 1, len(ref), dms_id, full_id, len(target_seq), selection)

        uni = fetch_uniprot_row(full_id, args.cache_dir)
        if not uni:
            log.warning("  no UniProt record for entry name %s; skipping", full_id)
            continue
        accession = uni.get("Entry", "")

        try:
            canonical_len = int(uni.get("Length", "0") or 0)
        except ValueError:
            canonical_len = 0
        # Annotation coords are 1-based into the canonical UniProt sequence;
        # DMS mutant positions are also 1-based but may be offset within it.
        try:
            raw = r.get("raw_mut_offset", "")
            offset = int(float(raw)) if raw and str(raw) != "nan" else 0
        except (ValueError, TypeError):
            offset = 0

        seq_len = canonical_len or len(target_seq)
        ranges_by_concept = concept_to_residue_set(
            uni, seq_len, categorical_options, binary_cols, paired_cols
        )
        if not ranges_by_concept:
            log.info("  no annotations match paired concepts")
            continue

        dms_path = args.dms_dir / dms_filename
        if not dms_path.exists():
            log.warning("  DMS file missing: %s", dms_path)
            continue
        dms = pd.read_csv(dms_path, usecols=["mutant"])
        positions_per_row: List[List[int]] = []
        for mutant in dms["mutant"]:
            parsed = parse_mutant_positions(str(mutant))
            positions_per_row.append([pos + offset for _, pos, _ in parsed])
        n_total_variants = len(dms)

        for concept, residue_set in ranges_by_concept.items():
            if concept not in pairings:
                # e.g. "Domain_any" — synthesized by InterPLM's parser but not
                # paired to any feature; skip.
                continue
            n_variants_in = sum(
                1 for ps in positions_per_row if any(p in residue_set for p in ps)
            )
            features = pairings[concept]
            rows.append({
                "DMS_id": dms_id,
                "UniProt_ID": full_id,
                "accession": accession,
                "target_seq_len": len(target_seq),
                "canonical_seq_len": canonical_len,
                "coarse_selection_type": selection,
                "selection_type": fn_class,
                "matched_concept": concept,
                "n_candidate_features": len(features),
                "top_feature": features[0]["feature"],
                "top_feature_f1": round(features[0]["f1"], 4),
                "top_feature_f1_per_domain": round(features[0]["f1_per_domain"], 4),
                "candidate_features": ";".join(str(f["feature"]) for f in features),
                "candidate_f1_per_domain": ";".join(
                    f"{f['f1_per_domain']:.3f}" for f in features
                ),
                "n_concept_residues": len(residue_set),
                "n_total_variants": n_total_variants,
                "n_variants_in_concept": n_variants_in,
                "frac_variants_in_concept": (
                    n_variants_in / n_total_variants if n_total_variants else 0.0
                ),
                "raw_mut_offset": offset,
            })
            log.info(
                "  %s -> %d/%d variants in concept (top feat %d, F1=%.3f)",
                concept, n_variants_in, n_total_variants,
                features[0]["feature"], features[0]["f1"],
            )

    out = pd.DataFrame(rows)
    if out.empty:
        log.warning("No matches found.")
        out.to_csv(args.output, index=False)
        return

    out = out.sort_values(
        ["n_variants_in_concept", "top_feature_f1_per_domain"], ascending=[False, False]
    ).reset_index(drop=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    log.info("Wrote %d (assay, concept) matches to %s", len(out), args.output)

    n_assays = out["DMS_id"].nunique()
    log.info("Summary: %d unique assays matched (%d in reference)", n_assays, len(ref))
    by_class = out.groupby("coarse_selection_type")["DMS_id"].nunique()
    log.info("By coarse_selection_type:\n%s", by_class.to_string())


if __name__ == "__main__":
    sys.exit(main())
