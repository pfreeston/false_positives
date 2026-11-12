# Predicted structures: hydrolase–effector complexes (P69B, Pip1)

Predicted structures and per-model confidence data for ten tomato hydrolase–effector
complexes, each predicted with seven method/database combinations. These are the
structures behind the AlphaJudge interface metrics reported in the manuscript.
The analysis code is in the accompanying GitHub repository (`8_zenodo/` builds
this deposit).

**Complexes (chain A = hydrolase, chain B = effector):** P69B × Ecp36, Epi1, EpiC1,
EpiC2B, SSP1, Six15, TIL; Pip1 × Epi1, EpiC1, EpiC2B.

**Sequences:** mature constructs, no prodomain and no signal peptide. P69B 635 aa
(`LSLHTT…`), Pip1 219 aa, Ecp36 54 aa without the FLAG tag, EpiC2B 104 aa. Exact
inputs are in `sequences.zip`. The AF3 server runs were submitted separately; their
exact sequences are in each `*_job_request.json` (P69B is 634 aa, `SLHTT…`, in eight
of them).

## Files

| File | Method | Per complex |
|---|---|---|
| `boltz2.zip` | Boltz-2 v2.2.1, MSAs from the ColabFold MMseqs2 server | 5 models: CIF, `confidence_*.json`, `pae_*.npz`, `plddt_*.npz` |
| `boltz2_reduced_db.zip` | Boltz-2 v2.2.1, MSAs from the AF3 data pipeline (reduced BFD), paired by top UniProt hit per taxon | as above |
| `af3_cluster.zip` | AlphaFold 3 (local, commit `dcaa19e`), own data pipeline (UniRef90, MGnify, small BFD, UniProt) | 1 seed × 5 samples: `*_model.cif`, `*_confidences.json`, `*_summary_confidences.json`; `*_ranking_scores.csv`; terms of use |
| `af3_cluster_full_db.zip` | AlphaFold 3, same seed as `af3_cluster`, MSAs taken from the AFM full-database run (unpaired: UniRef90 + BFD/UniRef30 + MGnify; paired: UniProt); templates searched by AF3 | as above |
| `af3_server.zip` | AlphaFold Server (October 2026) | 1 seed × 5 samples as above, plus `*_job_request.json` and the server terms of use |
| `afm_full_db.zip` | AlphaFold-Multimer v3 (AlphaFold 2.3.1), `full_dbs` (BFD + UniRef30, UniRef90, MGnify, UniProt) | 25 predictions (5 models × 5 seeds, base seed 20261108, template cut-off 2026-11-08), all Amber-relaxed: `relaxed_*.pdb`, `confidence_*.npz`, `ranking_debug.json` |
| `afm_reduced_db.zip` | AlphaFold-Multimer v3, `reduced_dbs` (small BFD) | as above |
| `sequences.zip` | input FASTA per complex | |
| `alphajudge.zip` | AlphaJudge 1.4.5 interface metrics, one CSV per method, every model | |
| `MANIFEST.tsv` | path, size and SHA-256 of every file inside the zips | |

Each zip contains `<method>/<Complex>/…`.

### AFM confidence files

The AFM result pickles (about 115 GB in all) are not deposited. Each
`confidence_<model>.npz` holds the values AlphaJudge reads from them:
`plddt` (per residue), `pae` (N × N, Å), `ptm`, `iptm` and `ranking_confidence`
(0.8·ipTM + 0.2·pTM). The manuscript figures use each complex's top five of the 25
predictions by `ranking_confidence`, so that every method contributes five models.

## Terms of use

AlphaFold 3 and AlphaFold Server outputs (`af3_*.zip`) are subject to Google
DeepMind's output terms of use, which are included in each complex folder and must
stay with the files. Among other things, they restrict use to non-commercial
purposes. Boltz-2 and AlphaFold-Multimer outputs carry no such terms.

## Citation

Please cite the accompanying manuscript (reference to be added) and this deposit.
