#!/bin/bash
#SBATCH -p lrz-hgx-h100-94x4
#SBATCH --gres=gpu:1
#SBATCH -t 24:00:00
#SBATCH -o logs/proteingym_full_feat_%j.out
#SBATCH -e logs/proteingym_full_feat_%j.err

# Full 8192-feature ProteinGym scan: per-(assay, feature) Spearman over all 217 assays
# (213 after the 2048-aa cap). Auxfix re-run; eval-independent matches built inline from
# DMS_substitutions.csv. Full rationale (memory, walltime, cap choice) in
# documentation/experiments/03-proteingym.md.
#
# Config notes:
#   --max_variants 200000  keeps 211/213 assays whole (only SPG1_Olson/HIS7 subsampled).
#   -t 24:00:00            ~14.5-18.2 h est; don't over-pad (longer -t queues later).
#   no --mem               ~96 GB default host RAM is plenty (peak ~35 GB, leak-free).
#   resumable              per-assay parquet + --skip_existing; a kill just re-submits.

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
