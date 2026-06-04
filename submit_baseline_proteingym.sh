#!/bin/bash
#SBATCH -p lrz-hgx-h100-94x4
#SBATCH --gres=gpu:1
#SBATCH -t 14:00:00
#SBATCH -o logs/proteingym_baseline_full_feat_%j.out
#SBATCH -e logs/proteingym_baseline_full_feat_%j.err

# Baseline-crosscoder ProteinGym null (SAEBench / InterPLM convention).
#
# Runs the full 8192-feature per-(assay, feature) Spearman scan on the crosscoder trained on a
# RANDOMLY-INITIALIZED ProtT5 (the baseline crosscoder; trained + InterPLM-annotated in
# experiments 01 & 02, ~1 concept / avg F1 0.053 vs 0.278 for the real model). This is the proper
# null for the ProteinGym feature-fitness analysis: if the real crosscoder's best-of-N feature |rho|
# is not above the baseline crosscoder's, the signal is architecture + best-of-N selection rather
# than ProtT5's learned representations. (The within-model random-feature control in
# chance_control.py is a weaker, selection-only null.)
#
# Identical args to submit_full_feat.sh except the checkpoint and the output dir. Compare the
# resulting full_feature_spearman_baseline/ against full_feature_spearman/ with chance_control.py
# (point it at the baseline CSV) or a small best-of-N diff.

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

CROSSCODER_DIR="/workspace/model_checkpoints/crosscoder_l8192_k32_bs512_baseline_2026-05-09_11-50-43/final_epoch_0_step_2519836"
MATCHES="/workspace/data/proteingym/concept_matches.csv"
DMS_DIR="/workspace/data/DMS_ProteinGym_substitutions"
OUTPUT_DIR="/workspace/data/proteingym/full_feature_spearman_baseline"

export HF_HOME="/workspace/data/hf_home"
export PYTHONPATH="/workspace/InterPLM"

mkdir -p logs

echo "Starting BASELINE ProteinGym full-feature scan on $(hostname) at $(date)"
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
     python /workspace/scc/harvest_full_feature_spearman.py \
       --crosscoder_dir ${CROSSCODER_DIR} \
       --checkpoint ae_normalized.pt \
       --matches ${MATCHES} \
       --dms_dir ${DMS_DIR} \
       --output_dir ${OUTPUT_DIR} \
       --batch_size 16 \
       --max_variants 50000 \
       --max_seq_len 2048 \
       --device cuda \
       --skip_existing"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Baseline full-feature scan finished at $(date)"
echo "Total duration: $((DURATION / 3600))h $((DURATION % 3600 / 60))m $((DURATION % 60))s"
