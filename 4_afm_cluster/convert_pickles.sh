#!/bin/bash
#SBATCH --job-name=afm_convert
#SBATCH --array=0-9
#SBATCH --clusters=htc
#SBATCH --partition=short
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/convert_%A_%a.out
#
# Rewrite one complex's result_*.pkl.gz without jax objects (same values, plain
# numpy) using 5_alpha_judge/convert_afm_pickles.py, so AlphaJudge can read them on
# dwight, where no env has jax. Run on the cluster after compress_pickles.sh and
# BEFORE transferring, so only the converted pickles are copied; the jax originals
# stay here in <Complex>/original_jax_pkls/ (exclude that dir from the rsync).
#   sbatch convert_pickles.sh <run dir>     e.g. $DATA/false_positives_proj/afm_matched_reduced
# convert_afm_pickles.py must sit next to this script in the run dir.

set -euo pipefail
BASE="$1"
COMPLEXES=(
  P69B_Ecp36 P69B_Epi1 P69B_EpiC1 P69B_EpiC2B P69B_SSP1 P69B_Six15 P69B_TIL
  Pip1_Epi1 Pip1_EpiC1 Pip1_EpiC2B
)
d="$BASE/outputs/${COMPLEXES[$SLURM_ARRAY_TASK_ID]}"
echo "=== $d on $(hostname) at $(date)"
/usr/bin/python3 "$BASE/convert_afm_pickles.py" "$d"
n=$(ls "$d"/result_*.pkl.gz | wc -l)
b=$(ls "$d"/original_jax_pkls/result_*.pkl.gz | wc -l)
echo "$n converted pickles, $b jax originals kept"
[ "$n" -eq 25 ] && [ "$b" -eq 25 ]
