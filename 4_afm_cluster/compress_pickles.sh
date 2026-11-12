#!/bin/bash
#SBATCH --job-name=afm_gzip
#SBATCH --array=0-9
#SBATCH --clusters=htc
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/gzip_%A_%a.out
#
# gzip result_*.pkl (in place -> result_*.pkl.gz, which AlphaJudge reads) for
# one AFM complex, after its prediction task has finished. Chain per index:
#   sbatch --dependency=aftercorr:<afm job> compress_pickles.sh <run dir>
# e.g. compress_pickles.sh $DATA/false_positives_proj/afm_matched

set -euo pipefail
BASE="$1"
COMPLEXES=(
  P69B_Ecp36 P69B_Epi1 P69B_EpiC1 P69B_EpiC2B P69B_SSP1 P69B_Six15 P69B_TIL
  Pip1_Epi1 Pip1_EpiC1 Pip1_EpiC2B
)
d="$BASE/outputs/${COMPLEXES[$SLURM_ARRAY_TASK_ID]}"
echo "=== $d at $(date)"
n=$(ls "$d"/result_*.pkl 2>/dev/null | wc -l)
echo "$n uncompressed result pickles"
[ "$n" -gt 0 ] && pigz -p "$SLURM_CPUS_PER_TASK" "$d"/result_*.pkl
for f in "$d"/result_*.pkl.gz; do gzip -t "$f"; done
echo "$(ls "$d"/result_*.pkl.gz | wc -l) .pkl.gz, all pass gzip -t; $(du -sh "$d" | cut -f1) total"
