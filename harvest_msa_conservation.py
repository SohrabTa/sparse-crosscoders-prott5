"""
Per-position conservation (Shannon entropy) from ProteinGym MSAs, for Design (1)'s
conservation-matched null (`design1_concept_localized_fitness.py --conservation`).

For each DMS assay, ProteinGym ships the MSA used by its EVE/MSA-Transformer baselines
(`DMS_substitutions.csv` columns MSA_filename / MSA_start / MSA_end). We compute, per MSA
COLUMN that maps onto the target sequence, the Shannon entropy over the 20 AAs (gaps excluded),
optionally sequence-weighted by the ProteinGym `theta` reweighting. Higher entropy = more
variable = less conserved.

OUTPUT CSV (the schema Design 1 expects): DMS_id, pos, score
  pos   = 1-indexed position in the assay's TARGET sequence (aligned to MSA_start)
  score = Shannon entropy (nats) of that column over non-gap AAs; higher = more variable

BLOCKED ON DATA: the MSA files are NOT in the local data/ tree (checked 2026-06-10). Get them from
the ProteinGym release (https://proteingym.org / OATML-Markslab/ProteinGym -> "MSA & weights",
DMS_msa_files.zip, ~a few GB) into --msa_dir. Then this runs locally (no GPU). Until then,
Design 1's identity-restricted control (null C) is the active confound control; this adds null D.

Usage:
  uv run --with biopython --with pandas --with numpy python \
    repos/sparse-crosscoders-prott5/harvest_msa_conservation.py \
      --reference data/external/DMS_substitutions.csv \
      --msa_dir data/external/DMS_msa_files \
      --assays A0A2Z5U3Z0_9INFA_Doud_2016,C6KNH7_9INFA_Lee_2018 \
      --output data/proteingym/msa_conservation.csv

Disposition: COMMIT — feeds a reported null for Design 1 (run once MSAs are fetched).
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("msa_cons")

AAS = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AAS)}


def read_msa(path: Path):
    """Yield (id, seq) from an a2m/a3m/fasta MSA. a3m lowercase = inserts (dropped to match the
    query columns, the standard a3m->alignment convention)."""
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name is not None:
        yield name, "".join(buf)


def seq_weights(aln: np.ndarray, theta: float) -> np.ndarray:
    """ProteinGym/EVE-style reweighting: w_i = 1 / #{j : identity(i,j) > 1-theta}. O(N^2) — only
    affordable on a subsample (see --reweight/--max_weight_seqs). theta NaN/0 -> uniform."""
    n, L = aln.shape
    if not np.isfinite(theta) or theta <= 0:
        return np.ones(n)
    thresh = (1.0 - theta) * L
    w = np.ones(n)
    for i in range(n):
        ident = (aln[i] == aln).sum(axis=1)
        w[i] = 1.0 / np.maximum((ident > thresh).sum(), 1)
    return w


def column_entropies(aln: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Per-column Shannon entropy (nats) over the 20 AAs, gap (-1) excluded, sequence-weighted by w.
    Vectorized: 20 masks over the (N, L) alignment -> (L,) entropies."""
    L = aln.shape[1]
    counts = np.zeros((L, 20))
    for a in range(20):
        counts[:, a] = (w[:, None] * (aln == a)).sum(axis=0)
    tot = counts.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = counts / tot[:, None]
        ent = -np.nansum(np.where(p > 0, p * np.log(p), 0.0), axis=1)
    ent[tot <= 0] = np.nan
    return ent


def process_assay(dms_id: str, ref_row: pd.Series, msa_dir: Path,
                  reweight: bool, max_weight_seqs: int, rng: np.random.Generator) -> pd.DataFrame:
    msa_path = msa_dir / str(ref_row["MSA_filename"])
    if not msa_path.exists():
        log.warning("  %s: MSA missing (%s)", dms_id, msa_path); return None
    seqs = list(read_msa(msa_path))
    if not seqs:
        return None
    # keep only query columns: positions that are upper/gap in the FIRST (query) sequence (a3m convention)
    query = seqs[0][1]
    keep_cols = np.array([j for j, c in enumerate(query) if c == "-" or c.isupper()])
    upper = np.frombuffer(query.encode(), dtype="S1")
    qmap = np.array([AA_IDX.get(c.decode().upper(), -1) for c in upper[keep_cols]])  # query residue per kept col
    is_target_col = np.array([query[j] != "-" for j in keep_cols])                  # query gap -> no target pos

    def encode(s):
        return np.array([AA_IDX.get(s[j].upper(), -1) if s[j] != "-" else -1 for j in keep_cols])
    aln = np.vstack([encode(s) for _, s in seqs])           # (N, Lq), -1 = gap/non-AA
    n_full = aln.shape[0]

    if reweight:
        if n_full > max_weight_seqs:
            sub = rng.choice(n_full, size=max_weight_seqs, replace=False)
            aln = aln[sub]
            log.info("  %s: reweight on %d/%d subsampled seqs (O(N^2))", dms_id, len(sub), n_full)
        w = seq_weights(aln, float(ref_row.get("MSA_theta", np.nan)))
    else:
        w = np.ones(aln.shape[0])                            # unweighted entropy over all seqs (fast)

    ent = column_entropies(aln, w)                           # (Lq,)
    start = int(ref_row.get("MSA_start", 1))
    pos_counter = 0
    rows = []
    for jcol in range(len(keep_cols)):
        if not is_target_col[jcol]:
            continue
        rows.append({"DMS_id": dms_id, "pos": start + pos_counter, "score": float(ent[jcol])})
        pos_counter += 1
    log.info("  %s: %d seqs (reweight=%s), %d target positions", dms_id, n_full, reweight, len(rows))
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--msa_dir", type=Path, required=True)
    ap.add_argument("--assays", type=str, default=None)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--reweight", action="store_true",
                    help="apply theta sequence-reweighting (O(N^2); subsamples to --max_weight_seqs). "
                         "Default off: plain Shannon entropy, a standard conservation measure and a "
                         "fine matching covariate, computed over all sequences in seconds.")
    ap.add_argument("--max_weight_seqs", type=int, default=8000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    log.info("reweight=%s  max_weight_seqs=%d  seed=%d", args.reweight, args.max_weight_seqs, args.seed)
    ref = pd.read_csv(args.reference).set_index("DMS_id")
    ids = args.assays.split(",") if args.assays else list(ref.index)
    out = []
    for dms_id in ids:
        if dms_id not in ref.index:
            log.warning("%s not in reference", dms_id); continue
        df = process_assay(dms_id, ref.loc[dms_id], args.msa_dir,
                           args.reweight, args.max_weight_seqs, rng)
        if df is not None:
            out.append(df)
    if not out:
        log.error("no conservation computed (MSA files present in --msa_dir?)"); return
    res = pd.concat(out, ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.output, index=False)
    log.info("wrote %s (%d rows, %d assays)", args.output, len(res), res["DMS_id"].nunique())


if __name__ == "__main__":
    main()
