#!/bin/bash
#SBATCH -p lrz-hgx-h100-94x4
#SBATCH --gres=gpu:1
#SBATCH -t 15:00:00
#SBATCH -o logs/proteingym_pooled_%j.out
#SBATCH -e logs/proteingym_pooled_%j.err

# Harvest full-sequence pooled crosscoder metrics for ProteinGym matched assays.
# Computes 8 gated + 8 dense per-(variant, feature) summaries (Δ at position
# plus sequence-pooled mean/max of |Δ| and act_mut), then a per-(assay,
# concept, feature) Spearman summary table — used to test the hypothesis that
# sequence-pooled signal (g_pool_mean) recovers fitness correlation that the
# per-residue Δ misses on the full 135-assay set.

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

CROSSCODER_DIR="/workspace/model_checkpoints/crosscoder_l8192_k32_bs512_full_2026-03-12_06-03-41/crashed_epoch_0_step_2519836"
PAIRINGS="/workspace/InterPLM/results/crosscoder_eval/uniprotkb_modern_score45_67k/test_counts/heldout_all_top_pairings.csv"
MATCHES="/workspace/data/proteingym/concept_matches.csv"
DMS_DIR="/workspace/data/DMS_ProteinGym_substitutions"
OUTPUT_DIR="/workspace/data/proteingym/pooled_metrics"

export HF_HOME="/workspace/data/hf_home"
export PYTHONPATH="/workspace/InterPLM"

mkdir -p logs

echo "Starting ProteinGym pooled-metric harvest on $(hostname) at $(date)"
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
     python /workspace/scc/harvest_pooled_metrics.py \
       --crosscoder_dir ${CROSSCODER_DIR} \
       --checkpoint ae_normalized.pt \
       --pairings ${PAIRINGS} \
       --matches ${MATCHES} \
       --dms_dir ${DMS_DIR} \
       --output_dir ${OUTPUT_DIR} \
       --batch_size 16 \
       --max_seq_len 2048 \
       --device cuda \
       --compute_summary"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Pooled harvest job finished at $(date)"
echo "Total duration: $((DURATION / 3600))h $((DURATION % 3600 / 60))m $((DURATION % 60))s"
