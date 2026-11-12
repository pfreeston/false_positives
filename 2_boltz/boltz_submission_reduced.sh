#!/bin/bash
#SBATCH --clusters=htc
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --constraint="gpu_gen:Ampere|gpu_gen:Hopper|gpu_gen:Lovelace"
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --job-name=boltz_reduced
#SBATCH --array=0-9
#SBATCH --output=logs/boltz_reduced_%A_%a.out
#
# Boltz-2 on the AFM-matched sequences, but with MSAs from the AF3 data
# pipeline (UniRef90, MGnify, small BFD, UniProt = AlphaFold reduced_dbs set)
# instead of the ColabFold MMseqs2 server. Inputs are made by
# af3_msa_to_boltz.py. Otherwise identical to boltz_submission.sh.
#
# Run from $DATA/false_positives_proj/boltz_reduced (inputs/ next to this script).

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate boltz_env

export PIP_CACHE_DIR=$DATA/pip_cache
export TMPDIR=$DATA/tmp

YAMLS=(inputs/*.yaml)
INPUT=${YAMLS[$SLURM_ARRAY_TASK_ID]}

echo "Running $INPUT"
nvidia-smi --query-gpu=name,memory.total --format=csv

boltz predict "$INPUT" \
    --cache $DATA/boltz_cache \
    --out_dir outputs \
    --diffusion_samples 5
