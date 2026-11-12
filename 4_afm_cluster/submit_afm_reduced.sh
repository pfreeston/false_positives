#!/bin/bash
#SBATCH --job-name=afm_reduced
#SBATCH --array=0-9
#SBATCH --clusters=htc
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --constraint="gpu_gen:Ampere"
# Ampere only: the AF 2.3.1 module (CUDA 11.7) cannot compile for Hopper (H100,
# CC 9.0; "ptxas does not support CC 9.0") and inference crashes there with
# "XlaRuntimeError: INTERNAL: Could not find the corresponding function".
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/afm_%A_%a.out
#SBATCH --error=logs/afm_%A_%a.err
#
# AlphaFold-Multimer v3 (AF 2.3.1), REDUCED databases (small BFD, as the
# 4_afm screen used), full pipeline (MSA + templates + inference),
# on the ten false-positive complexes using the AFM-matched sequences
# (1_sequence_processing/afm_matched/, same constructs as the Boltz/AF3 matched
# runs: 635 aa mature P69B, 219 aa Pip1, Ecp36 without the FLAG+linker).
#
# Unlike the screen (1 prediction per model, pickles deleted), this keeps
# 5 seeds x 5 models = 25 predictions per complex with all result_*.pkl
# (PAE etc.) and all 25 relaxed, so convergence of the binding pose can be checked.
#
# Submit from $DATA/false_positives_proj/afm_matched_reduced (fastas/ next to this script):
#   sbatch submit_afm_reduced.sh               # all ten
#   sbatch --array=1 submit_afm_reduced.sh     # just P69B_Epi1
#
# Walltime: reduced_dbs only swaps HHblits (BFD+UniRef30) for Jackhmmer on
# small BFD; the UniRef90/MGnify/UniProt Jackhmmer searches, which dominated
# the full_dbs run (~3-3.7 h MSA under shared-storage I/O), still run. Expect
# ~2-3.5 h MSA + ~1.5-2 h for 25 predictions and relax on an H100/A100.
# 12 h (partition max) leaves headroom.

set -euo pipefail

BASE="$DATA/false_positives_proj/afm_matched_reduced"
D=/data/dops-plant-genomics/bras4304/bras4304/alphafold/v2_3_1_databases

FASTAS=(
  P69B_Ecp36 P69B_Epi1 P69B_EpiC1 P69B_EpiC2B P69B_SSP1 P69B_Six15 P69B_TIL
  Pip1_Epi1 Pip1_EpiC1 Pip1_EpiC2B
)
name="${FASTAS[$SLURM_ARRAY_TASK_ID]}"
# Fixed, recorded base seed (the screen never recorded its seeds).
SEED=20261108

echo "=== $name on $(hostname) at $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

module --ignore-cache load AlphaFold/2.3.1-foss-2022a-CUDA-11.7.0
export ALPHAFOLD_DATA_DIR="$D"
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
# EasyBuild's AlphaFold reads these for the MSA tools' thread counts.
export ALPHAFOLD_JACKHMMER_N_CPU=$SLURM_CPUS_PER_TASK
export ALPHAFOLD_HHBLITS_N_CPU=$SLURM_CPUS_PER_TASK

RUN_AF=$(command -v run_alphafold.py || command -v alphafold)
echo "Using $RUN_AF"

cd "$BASE"
mkdir -p outputs logs

"$RUN_AF" \
  --fasta_paths="fastas/${name}.fasta" \
  --output_dir=outputs \
  --data_dir="$D" \
  --model_preset=multimer \
  --db_preset=reduced_dbs \
  --uniref90_database_path="$D/uniref90/uniref90.fasta" \
  --mgnify_database_path="$D/mgnify/mgy_clusters_2022_05.fa" \
  --small_bfd_database_path="$D/small_bfd/bfd-first_non_consensus_sequences.fasta" \
  --pdb_seqres_database_path="$D/pdb_seqres/pdb_seqres.txt" \
  --uniprot_database_path="$D/uniprot/uniprot.fasta" \
  --template_mmcif_dir="$D/pdb_mmcif/mmcif_files" \
  --obsolete_pdbs_path="$D/pdb_mmcif/obsolete.dat" \
  --max_template_date=2026-11-08 \
  --num_multimer_predictions_per_model=5 \
  --random_seed=$SEED \
  --run_relax=true \
  --use_gpu_relax=true \
  --use_precomputed_msas=true

echo "=== done $name at $(date)"
