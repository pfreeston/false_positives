#!/bin/bash
#SBATCH --clusters=htc
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --constraint="gpu_gen:Ampere|gpu_gen:Hopper|gpu_gen:Lovelace"
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --job-name=boltz
#SBATCH --array=0-9
#SBATCH --output=logs/boltz_%A_%a.out

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate boltz_env

export PIP_CACHE_DIR=$DATA/pip_cache
export TMPDIR=$DATA/tmp

YAMLS=(boltz_inputs_afm_matched/*.yaml)
INPUT=${YAMLS[$SLURM_ARRAY_TASK_ID]}

echo "Running $INPUT"
nvidia-smi --query-gpu=name,memory.total --format=csv

boltz predict "$INPUT" \
    --use_msa_server \
    --cache $DATA/boltz_cache \
    --out_dir boltz_results_afm_matched \
    --diffusion_samples 5