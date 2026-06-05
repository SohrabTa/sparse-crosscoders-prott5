"""Check 1 — is the per-residue mutation-effect m[p] predictable from the WT amino acid alone?

The per-residue feature-vs-DMS readout (MotifAE's method) correlates a feature's per-residue WT
activation a[p] against m[p] = mean mutation effect at position p. We found a random-init crosscoder
nearly matches the trained one on it (0.452 vs 0.472), suggesting the readout reflects an
amino-acid-identity / composition prior rather than learned biology. This script tests that directly,
with NO model: how well does the WT amino-acid IDENTITY at p predict m[p]?

Per assay (single mutants; positions with >= MIN_SUBS substitutions; assays with >= MIN_POS such
positions — same filter as harvest_per_residue_wt.py):
  - m[p]   = mean DMS_score over substitutions at p
  - aa[p]  = WT amino acid at p
  - LOO prediction pred[p] = mean of m[q] over q != p with aa[q] == aa[p]  (leave-one-out, so the
    prediction never sees position p itself)
  - Spearman(pred, m) across positions  -> how much of m[p] identity alone explains
  - eta   = sqrt(between-AA variance / total variance) of m[p] grouped by aa  (correlation ratio)

If the median identity-Spearman is ~0.4 (comparable to the per-residue feature numbers 0.45-0.53),
the per-residue readout is largely a non-learned identity prior, and trained ~ random-init follows.

Reads: data/proteingym/concept_matches.csv, data/DMS_ProteinGym_substitutions/*.csv  (no model)
Writes: data/proteingym/per_residue_identity_check.csv ; prints medians.

Usage:  uv run --with pandas --with numpy --with scipy python \
          repos/sparse-crosscoders-prott5/per_residue_identity_check.py
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path("/Users/sohrab.tawana/private/crosscoder")
PG = ROOT / "data/proteingym"
DMS_DIR = ROOT / "data/DMS_ProteinGym_substitutions"
MATCHES = PG / "concept_matches.csv"
OUT = PG / "per_residue_identity_check.csv"
MUT_RE = re.compile(r"^([A-Z])(\d+)([A-Z\*])$")
MIN_SUBS = 10
MIN_POS = 5
MAX_SEQ_LEN = 2048


def main() -> None:
    matches = pd.read_csv(MATCHES)
    matches = matches[matches["target_seq_len"] <= MAX_SEQ_LEN]
    assay_ids = sorted(matches["DMS_id"].unique())

    rows = []
    for dms_id in assay_ids:
        p = DMS_DIR / f"{dms_id}.csv"
        if not p.exists():
            continue
        dms = pd.read_csv(p)
        parsed = dms["mutant"].astype(str).str.match(MUT_RE) & ~dms["mutant"].astype(str).str.contains(":")
        dms = dms[parsed].copy()
        if len(dms) == 0:
            continue
        ext = dms["mutant"].astype(str).str.extract(MUT_RE)
        dms["wt_aa"] = ext[0]
        dms["pos"] = ext[1].astype(int)
        g = dms.groupby("pos").agg(m=("DMS_score", "mean"), n=("DMS_score", "size"),
                                   aa=("wt_aa", "first"))
        g = g[g["n"] >= MIN_SUBS]
        if len(g) < MIN_POS:
            continue
        # leave-one-out AA-mean prediction
        aa_sum = g.groupby("aa")["m"].transform("sum")
        aa_cnt = g.groupby("aa")["m"].transform("size")
        loo = (aa_sum - g["m"]) / (aa_cnt - 1)          # NaN where the AA appears once
        valid = loo.notna()
        if valid.sum() < MIN_POS:
            continue
        rho = spearmanr(loo[valid], g["m"][valid]).statistic
        # correlation ratio eta (between-AA / total)
        grand = g["m"].mean()
        ss_tot = ((g["m"] - grand) ** 2).sum()
        ss_between = g.groupby("aa")["m"].apply(lambda s: len(s) * (s.mean() - grand) ** 2).sum()
        eta = float(np.sqrt(ss_between / ss_tot)) if ss_tot > 0 else np.nan
        rows.append({"DMS_id": dms_id, "n_pos": int(len(g)),
                     "identity_loo_spearman": float(rho) if np.isfinite(rho) else np.nan,
                     "eta_aa": eta})

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    v = out["identity_loo_spearman"].abs().dropna()
    e = out["eta_aa"].dropna()
    print(f"check 1 — identity-only predictor of per-residue m[p], {len(out)} assays")
    print(f"  median |Spearman(LOO AA-mean, m)| = {v.median():.3f}   (>0.3={int((v>0.3).sum())}, >0.5={int((v>0.5).sum())})")
    print(f"  median eta (AA correlation ratio) = {e.median():.3f}")
    print(f"  [compare: per-residue feature readout — trained 0.47, random-init 0.45]")


if __name__ == "__main__":
    main()
