#!/bin/bash
#SBATCH -p lrz-hgx-h100-94x4
#SBATCH --gres=gpu:1
#SBATCH -t 24:00:00
#SBATCH -o logs/proteingym_full_feat_%j.out
#SBATCH -e logs/proteingym_full_feat_%j.err

# Full 8192-feature scan: per-(assay, feature) Spearman across ALL 217 ProteinGym
# substitution assays (213 after the 2048-aa cap drops 4 long proteins). The §5.6
# supplementary: "do ANY of the 8192 features track fitness across the benchmark?"
# Output: ~213 × 8192 ≈ 1.7M (assay, feature, metric) rows. Spearman on
# pool_mean_abs / delta_at_pos only — NO AUC/MCC/NDCG (same confounded signal,
# adds no validity; see exp-03 eta=0.575 prior note).
#
# AUXFIX RE-RUN: defaults are the auxfix full crosscoder + the EVAL-INDEPENDENT
# assay list (built here from DMS_substitutions.csv, NOT from concept_matches.csv,
# so the assay set does not condition on which features InterPLM paired).
#
# MEMORY: no --mem needed. A 1-GPU job here gets ~96 GB HOST RAM by default
# (JobDefaults DefMemPerGPU=95830 MB — host RAM granted per GPU; it equals the H100's
# 94 GB VRAM only because the admins set it so, it is NOT the VRAM). This harvest is
# leak-free (per-assay `del`+gc, accumulators freed each assay) and the SAME code
# completed clean as job 5665789 at MaxRSS 20 GB with no --mem. The earlier OOM (job
# 5661387, MaxRSS 92 GB) was a different script, harvest_pooled_metrics.py, with a
# since-fixed leak (commit 8ce277a: streaming writer + per-assay GC). Host peak here =
# ProtT5/buffers (~15 GB) + per-assay accumulator (3 × max_variants × 8192 × 4 B):
#     cap=50000  -> accum ~5 GB  -> peak ~20 GB (measured, 5665789)
#     cap=200000 -> accum ~20 GB -> peak ~35 GB   (fits the ~96 GB default comfortably)
#     cap=600000 -> accum ~59 GB -> peak ~74 GB   (still fits the default, tighter)
# GPU VRAM is set by batch_size × seq_len^2 (ProtT5 attention), INDEPENDENT of
# max_variants; batch_size 16 already ran fine on the 94 GB H100 in 5665789.
#
# --max_variants 200000 keeps 211/213 assays WHOLE and subsamples only the two true
# mega-assays (SPG1_Olson 537k, HIS7_YEAST 496k), whose Spearman is rock-stable far
# below 200k. Walltime is calibrated from 5665789's MEASURED cost: 744,381 variants
# scanned in 7h29m (the 135-assay/50k job). Scaling to the new 213-assay set, two ways:
#     cap=200000 -> 1.81M scanned -> 14.5 h (by variant*seq_len) .. 18.2 h (by variant)
#     cap=600000 -> 2.44M scanned -> 19.1 h .. 24.5 h
# -t 24:00:00: ~1.3x over the conservative 18.2 h, ~1.65x over the realistic 14.5 h.
# Do NOT over-pad: requested walltime is a scheduling cost (backfill runs shorter jobs
# first), so a bigger -t just queues longer (see insights.md 2026-06). The cheap
# backstop against a too-tight ceiling is RESUMABILITY, not padding: each assay writes
# its own {DMS_id}__full_feat.parquet and --skip_existing skips finished ones, so a
# walltime kill loses only the in-flight assay and a re-submit finishes (and still runs
# the final CSV aggregation over the existing parquets). For cap=600000 use -t 30:00:00.
#
# batch_size 16 is the safe global value: it must hold for the 16 long-protein
# assays (>1024 aa, up to the 2048 cap) on the 94 GB H100. The short-seq majority
# (median 245 aa) could run at 24-32, but a single global batch must be safe for
# the longest, so 16. GPU peak is batch_size × seq_len^2 attention, INDEPENDENT of
# max_variants.

INTERPLM_DIR="/dss/dsshome1/08/ga25ley2/code/InterPLM"
CROSSCODE_DIR="/dss/dsshome1/08/ga25ley2/code/crosscode"
SCC_DIR="/dss/dsshome1/08/ga25ley2/code/sparse-crosscoders-prott5"
DATA_DIR="/dss/dssfs02/lwp-dss-0001/pn67na/pn67na-dss-0000/ga25ley2/data"
CKPT_DIR="/dss/dssfs02/lwp-dss-0001/pn67na/pn67na-dss-0000/ga25ley2/model_checkpoints"

MOUNTS="${INTERPLM_DIR}:/workspace/InterPLM"
MOUNTS="${MOUNTS},${DATA_DIR}:/workspace/data"
MOUNTS="${MOUNTS},${CKPT_DIR}:/workspace/model_checkpoints"
MOUNTS="${MOUNTS},${CROSSCODE_DIR}:/workspace/crosscode"
MOUNTS="${MOUNTS},${SCC_DIR}:/workspace/scc"

# --- re-run override (default = auxfix full run, the hand-in checkpoint) ---
#   export RERUN_FULL_CKPT=/workspace/model_checkpoints/<new_full>/<dir>
CROSSCODER_DIR="${RERUN_FULL_CKPT:-/workspace/model_checkpoints/crosscoder_l8192_k32_bs512_full_auxfix_2026-06-06_07-04-40/jumprelu_global_2519836}"
# EVAL-INDEPENDENT assay list (all 217 from the ProteinGym reference), built below.
# NOT concept_matches.csv: that is derived from InterPLM pairings, which would
# condition the assay set on the eval. See build_proteingym_matches.py.
REFERENCE="/workspace/data/DMS_substitutions.csv"
MATCHES="/workspace/data/proteingym_full_feature_matches.csv"
DMS_DIR="/workspace/data/DMS_ProteinGym_substitutions"
OUTPUT_DIR="/workspace/data/proteingym/full_feature_spearman_auxfix"

export HF_HOME="/workspace/data/hf_home"
export PYTHONPATH="/workspace/InterPLM"

mkdir -p logs

echo "Starting ProteinGym full-feature scan on $(hostname) at $(date)"
START_TIME=$(date +%s)

srun --container-image="nvcr.io/nvidia/pytorch:25.12-py3" \
     --container-mounts="${MOUNTS}" \
     --container-workdir="/workspace/InterPLM" \
     bash -c "uv venv --python 3.12 && \
     source .venv/bin/activate && \
     uv pip install -r requirements.txt && \
     uv pip install -e /workspace/crosscode && \
     uv pip install -e . && \
     uv pip install scipy pyarrow && \
     mkdir -p ${OUTPUT_DIR} && \
     python /workspace/scc/build_proteingym_matches.py \
       --reference ${REFERENCE} \
       --out ${MATCHES} && \
     python /workspace/scc/harvest_full_feature_spearman.py \
       --crosscoder_dir ${CROSSCODER_DIR} \
       --checkpoint ae_normalized.pt \
       --matches ${MATCHES} \
       --dms_dir ${DMS_DIR} \
       --output_dir ${OUTPUT_DIR} \
       --batch_size 16 \
       --max_variants 200000 \
       --max_seq_len 2048 \
       --device cuda \
       --skip_existing"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Full-feature scan finished at $(date)"
echo "Total duration: $((DURATION / 3600))h $((DURATION % 3600 / 60))m $((DURATION % 60))s"
