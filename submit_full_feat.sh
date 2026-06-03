#!/bin/bash
#SBATCH -p lrz-hgx-h100-94x4
#SBATCH --gres=gpu:1
#SBATCH -t 14:00:00
#SBATCH -o logs/proteingym_full_feat_%j.out
#SBATCH -e logs/proteingym_full_feat_%j.err

# Full 8192-feature scan: per-(assay, feature) Spearman across all 135 assays.
# Tests whether ProteinGym fitness signal lives in features beyond the InterPLM-
# paired candidates. Output: ~1.1M (assay, feature, metric) rows.
#
# Args parallel submit_pooled.sh wherever possible. The one deliberate divergence
# is --max_variants 50000: this script accumulates a (n_var, 8192) per-metric
# float32 array per assay in RAM, vs pooled.sh's (n_var, ~10). Without a cap
# SPG1_STRSG_Olson_2014 (534k variants) would peak at ~70 GB and OOM the H100
# after 10+ hours of compute — the cap holds peak accumulator ≤ 6.5 GB while
# subsampling only 3 mega-assays whose Spearman is rock-stable well below 50k.
# --skip_existing is set defensively so a re-submission resumes cleanly if
# anything trips the job mid-stream.
#
# Wall projection: ~10.5 h scoring (parallel to pooled.sh's 10h25m) + ~5 min
# aggregation + setup overhead → ~11 h realistic. 14h budget gives ~3h headroom.

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
MATCHES="/workspace/data/proteingym/concept_matches.csv"
DMS_DIR="/workspace/data/DMS_ProteinGym_substitutions"
OUTPUT_DIR="/workspace/data/proteingym/full_feature_spearman"

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
echo "Full-feature scan finished at $(date)"
echo "Total duration: $((DURATION / 3600))h $((DURATION % 3600 / 60))m $((DURATION % 60))s"
