#!/bin/bash
#SBATCH -p lrz-hgx-h100-94x4
#SBATCH --gres=gpu:1
#SBATCH -t 15:00:00
#SBATCH -o logs/proteingym_score_%j.out
#SBATCH -e logs/proteingym_score_%j.err

# Score crosscoder features against the ProteinGym DMS substitutions benchmark.
# By default scores ALL (assay, concept) pairs from --matches (138 assays /
# ~257 pairs / ~1.46M variants). The scoring script groups by assay and runs
# ProtT5 + crosscoder once per assay, slicing the cached activations for each
# matched concept (≈1.9× speedup over a naive per-(assay, concept) loop).
#
# Each per-(variant, feature) row is written to a parquet file as it is
# produced, and summary.csv is rewritten incrementally — a crashed job leaves
# usable partial results. Add --skip_existing to resume.
#
# Inputs (pre-staged on cluster):
#   - /workspace/data/DMS_ProteinGym_substitutions/*.csv  (downloaded from proteingym.org)
#   - /workspace/data/proteingym/concept_matches.csv      (output of match_proteingym_concepts.py;
#                                                          scp from local box, or regenerate via
#                                                          a CPU job — requires UniProt API access)
#   - /workspace/InterPLM/results/crosscoder_eval/uniprotkb_modern_score45_67k/test_counts/heldout_all_top_pairings.csv
#   - /workspace/model_checkpoints/crosscoder_l8192_k32_bs512_full_2026-03-12_06-03-41/crashed_epoch_0_step_2519836
#
# Output:
#   /workspace/data/proteingym/scoring/<DMS_id>__<concept>__per_variant.parquet (one per assay)
#   /workspace/data/proteingym/scoring/summary.csv

# Define Paths (HOST side)
INTERPLM_DIR="/dss/dsshome1/08/ga25ley2/code/InterPLM"
CROSSCODE_DIR="/dss/dsshome1/08/ga25ley2/code/crosscode"
SCC_DIR="/dss/dsshome1/08/ga25ley2/code/sparse-crosscoders-prott5"
DATA_DIR="/dss/dssfs02/lwp-dss-0001/pn67na/pn67na-dss-0000/ga25ley2/data"
CKPT_DIR="/dss/dssfs02/lwp-dss-0001/pn67na/pn67na-dss-0000/ga25ley2/model_checkpoints"

# Mounts: Host:Container
MOUNTS="${INTERPLM_DIR}:/workspace/InterPLM"
MOUNTS="${MOUNTS},${DATA_DIR}:/workspace/data"
MOUNTS="${MOUNTS},${CKPT_DIR}:/workspace/model_checkpoints"
MOUNTS="${MOUNTS},${CROSSCODE_DIR}:/workspace/crosscode"
MOUNTS="${MOUNTS},${SCC_DIR}:/workspace/scc"

# Run-specific paths (CONTAINER side)
CROSSCODER_DIR="/workspace/model_checkpoints/crosscoder_l8192_k32_bs512_full_2026-03-12_06-03-41/crashed_epoch_0_step_2519836"
PAIRINGS="/workspace/InterPLM/results/crosscoder_eval/uniprotkb_modern_score45_67k/test_counts/heldout_all_top_pairings.csv"
MATCHES="/workspace/data/proteingym/concept_matches.csv"
DMS_DIR="/workspace/data/DMS_ProteinGym_substitutions"
OUTPUT_DIR="/workspace/data/proteingym/scoring"

# Env
export HF_HOME="/workspace/data/hf_home"
export PYTHONPATH="/workspace/InterPLM"

mkdir -p logs

echo "Starting ProteinGym crosscoder scoring on $(hostname) at $(date)"
START_TIME=$(date +%s)

# Use Python 3.12 to satisfy crosscode requirements. Same install pattern as
# the other submit scripts: install InterPLM requirements + crosscode (editable)
# + InterPLM (editable), then add scipy + pyarrow for the scoring script.
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
     python /workspace/scc/score_proteingym_features.py \
       --crosscoder_dir ${CROSSCODER_DIR} \
       --checkpoint ae_normalized.pt \
       --pairings ${PAIRINGS} \
       --matches ${MATCHES} \
       --dms_dir ${DMS_DIR} \
       --output_dir ${OUTPUT_DIR} \
       --batch_size 16 \
       --max_seq_len 2048 \
       --device cuda"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo "ProteinGym scoring job finished at $(date)"
echo "Total duration: $((DURATION / 3600))h $((DURATION % 3600 / 60))m $((DURATION % 60))s"
