#!/bin/bash
#SBATCH --job-name=af3_fulldb
#SBATCH --array=0-9
#SBATCH --clusters=htc
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --constraint="gpu_gen:Ampere|gpu_gen:Hopper"
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/af3_fulldb_%A_%a.out
#SBATCH --error=logs/af3_fulldb_%A_%a.err
#
# AF3 on full-database MSAs: the per-chain MSAs from the AFM full_dbs run
# (UniRef90 + BFD/UniRef30 HHblits + MGnify, UniProt for pairing) are put into
# the AF3 JSON by build_af3_fulldb_json.py; AF3 then skips its own (small-BFD)
# MSA search and only runs template search + inference. Seeds are the AF3
# matched run's, so only the MSAs/templates differ from that run.
#
# Submit from $DATA/false_positives_proj/af3_fulldb, chained per array index to
# the AFM full_dbs job so each complex starts once its MSAs exist:
#   sbatch --dependency=aftercorr:9307751 submit_af3_fulldb.sh

set -euo pipefail

BASE="$DATA/false_positives_proj/af3_fulldb"
AFM="$DATA/false_positives_proj/afm_matched/outputs"
SEEDS="$DATA/false_positives_proj/Custom-AlphaFold-3-Workflow/work_afm_matched/json/multimers"
IMG=$DATA/af3/containers/alphafold3_custom.sif
MODEL_DIR=/apps/common/commercial/AlphaFold3/model
DB_DIR=/apps/datasets/alphafold3/public_databases

COMPLEXES=(
  P69B_Ecp36 P69B_Epi1 P69B_EpiC1 P69B_EpiC2B P69B_SSP1 P69B_Six15 P69B_TIL
  Pip1_Epi1 Pip1_EpiC1 Pip1_EpiC2B
)
name="${COMPLEXES[$SLURM_ARRAY_TASK_ID]}"
seed_json="$(echo "$name" | tr 'A-Z' 'a-z' | sed 's/_/__/').json"   # P69B_Epi1 -> p69b__epi1.json

echo "=== $name on $(hostname) at $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
cd "$BASE"
mkdir -p inputs outputs logs

# 1. Build the AF3 JSON from the AFM full_dbs MSAs.
singularity exec \
  --bind "$AFM/$name/msas":/root/afm_msas \
  --bind "$SEEDS":/root/seeds \
  --bind "$BASE":/root/work \
  "$IMG" \
  python3 /root/work/build_af3_fulldb_json.py \
    --name "$name" --msas /root/afm_msas \
    --seed_json "/root/seeds/$seed_json" \
    --out "/root/work/inputs/$name.json"

# 2. Template search (MSAs supplied, so no MSA search) + inference.
singularity exec \
  --nv \
  --env JAX_PLATFORMS=cuda \
  --bind "$DB_DIR":/root/public_databases \
  --bind "$MODEL_DIR":/root/models \
  --bind "$BASE/inputs":/root/af_input \
  --bind "$BASE/outputs":/root/af_output \
  "$IMG" \
  python3 /app/alphafold/run_alphafold.py \
    --json_path="/root/af_input/$name.json" \
    --db_dir=/root/public_databases \
    --model_dir=/root/models \
    --output_dir=/root/af_output \
    --run_data_pipeline=true \
    --run_inference=true

echo "=== done $name at $(date)"
