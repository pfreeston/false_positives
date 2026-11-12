# Cluster setup

All predictions except AlphaFold Server were run on the Oxford ARC cluster (Slurm
cluster `htc`, partition `short`, one GPU per job), November 2026. Working
directory on the cluster: `$DATA/false_positives_proj/`.
Every method used the same mature, AFM-matched sequences
(`1_sequence_processing/afm_matched/`; per-complex FASTAs in `4_afm_cluster/fastas/`).

| Method | Software | MSAs | Templates | Models per complex | Script |
|---|---|---|---|---|---|
| Boltz-2 | Boltz 2.2.1 (conda env `boltz_env`) | ColabFold MMseqs2 server (`--use_msa_server`) | none | 5 diffusion samples | `2_boltz/boltz_submission.sh` |
| Boltz-2 reduced-DB | Boltz 2.2.1 | AF3 data-pipeline MSAs, converted by `2_boltz/af3_msa_to_boltz.py` | none | 5 diffusion samples | `2_boltz/boltz_submission_reduced.sh` |
| AF3 cluster | AlphaFold 3, commit `dcaa19e` (2024-12-03), Singularity image `alphafold3_custom.sif` | AF3 data pipeline (Jackhmmer: UniRef90, MGnify, small BFD, UniProt) | AF3 template search (PDB seqres 2022-09-28) | 1 seed × 5 samples | Snakemake workflow, see below |
| AF3 cluster full-DB | as AF3 cluster, same seeds | AFM `full_dbs` MSAs (UniRef90, BFD + UniRef30 HHblits, MGnify; UniProt paired), built into the AF3 JSON by `3_af3/cluster/fulldb/build_af3_fulldb_json.py` | AF3 template search | 1 seed × 5 samples | `3_af3/cluster/fulldb/submit_af3_fulldb.sh` |
| AF3 server | AlphaFold Server (October 2026) | server | server | 1 seed × 5 samples | — (web) |
| AFM full-DB | AlphaFold 2.3.1, EasyBuild module `AlphaFold/2.3.1-foss-2022a-CUDA-11.7.0`, multimer v3 weights | `--db_preset=full_dbs` | PDB mmCIF, max template date 2026-11-08 | 5 models × 5 seeds, all Amber-relaxed | `4_afm_cluster/submit_afm.sh` |
| AFM reduced-DB | as AFM full-DB | `--db_preset=reduced_dbs` (small BFD) | as above | as above | `4_afm_cluster/submit_afm_reduced.sh` |

AFM base seed: `--random_seed=20261108`. AF3 model seeds are in each input JSON
(e.g. 153946889 for P69B_Epi1) and in the deposited AF3 outputs on Zenodo.

## Databases

**AlphaFold 3** (`/apps/datasets/alphafold3/public_databases`, ARC shared copy):
UniRef90 2022-05, MGnify clusters 2022-05, small BFD (first non-consensus
sequences), UniProt 2021-04, PDB seqres 2022-09-28 and mmCIF files. AF3 model
weights from `/apps/common/commercial/AlphaFold3/model`.

**AlphaFold 2.3.1** (`.../alphafold/v2_3_1_databases`, group copy): UniRef90,
MGnify clusters 2022-05, BFD (`bfd_metaclust_clu_complete_id30_c90_final_seq`),
UniRef30 2021-03, small BFD, UniProt (downloaded 2023-02), PDB seqres (2024-05)
and mmCIF files (2024-07).

## AF3 cluster run (Snakemake)

The AF3 cluster run used the Custom-AlphaFold-3-Workflow
(`github.com/pfreeston/Custom-AlphaFold-3-Workflow`, private), commit `15b807f`
(2026-11-04), in multimer mode with all-vs-all pairing of the plant (P69B, Pip1)
and pathogen (effector) groups. The run-specific config and target list are copied
into `3_af3/cluster/workflow/`; the input FASTAs are
`3_af3/cluster/fastas_afm_matched/` (A = hydrolases, B = effectors). The workflow
environment is `env/af3-custom-workflow.yaml` in that repository (Python 3.11,
snakemake-minimal, snakemake-executor-plugin-slurm).

```
snakemake --snakefile workflow/Snakefile --executor slurm \
  --configfile config/false_positives.yaml --jobs 20 \
  $(cat targets_afm_matched.txt)
```

(The command was reconstructed from the workflow README and the run's config; the
Snakemake log records only that `config/project.yaml` was extended from the command
line.) The workflow writes lowercase output names (`p69b__epi1`);
`5_alpha_judge/standardise_names.py` renames them to the canonical `P69B_Epi1`.

## GPUs

The AF 2.3.1 module is built against CUDA 11.7, which cannot compile for Hopper
(H100): jobs fail with `ptxas does not support CC 9.0` and
`XlaRuntimeError: INTERNAL: Could not find the corresponding function`. The AFM
scripts therefore request `gpu_gen:Ampere` only. AF3 and Boltz run on Ampere or
Hopper.

## After the runs

1. `4_afm_cluster/compress_pickles.sh <run dir>`: gzip the AFM result pickles.
2. `4_afm_cluster/convert_pickles.sh <run dir>`: rewrite the pickles' old-JAX
   arrays as NumPy (`5_alpha_judge/convert_afm_pickles.py`), so AlphaJudge can read
   them without JAX. Originals are kept in `original_jax_pkls/` on the cluster.
3. Transfer the converted outputs (excluding `original_jax_pkls/`), then run
   AlphaJudge 1.4.5 per method:
   `alphajudge -r <run dir> --models_to_analyse all --force_recompute -o <summary.csv>`.
   Summaries are in `5_alpha_judge/`.
