# False positives in hydrolase–effector complex prediction

Code and data for the analysis of predicted complexes between two tomato
hydrolases (the subtilase P69B and the papain-like cysteine protease Pip1) and
pathogen effectors. Each of ten complexes (four positive controls with
lab-confirmed interactions, six negatives that could not be confirmed) was
predicted with seven method/MSA-database combinations and scored with AlphaJudge.

**Predicted structures are not in this repository.** They are deposited on Zenodo
(DOI to be added); `8_zenodo/` builds that deposit.

## Layout

| Folder | Contents |
|---|---|
| `1_sequence_processing/` | sequences: precursor, mature (SignalP 6 trimmed) and the AFM-matched constructs used for every run |
| `2_boltz/` | Boltz-2 inputs and Slurm scripts (ColabFold MSAs; AF3-pipeline MSAs via `af3_msa_to_boltz.py`) |
| `3_af3/cluster/` | AF3 input FASTAs, Snakemake run config (`workflow/`), full-database MSA variant (`fulldb/`) |
| `4_afm_cluster/` | AlphaFold-Multimer full- and reduced-database Slurm scripts, inputs, pickle post-processing |
| `5_alpha_judge/` | AlphaJudge summaries (one CSV per method, every model), name standardisation, consensus tiers |
| `7_manuscript/` | manuscript tables, figure scripts and the CSVs behind them; see its README |
| `8_zenodo/` | builds the structure deposit; README and draft metadata |

`CLUSTER_SETUP.md` records software versions, databases, seeds and the exact
cluster commands for every method.

## Reproducing the tables and figures

```
conda env create -f environment.yml   # or any Python with pandas + matplotlib
cd 7_manuscript
python3 collect_manuscript_csvs.py && python3 compare_iptm_ptm.py && \
python3 make_iptm_table.py && python3 make_metric_figures.py && \
python3 make_metric_grid.py && python3 make_derived_dotplots.py && \
python3 make_metric_bars.py && python3 make_pdf.py
```

These need only the committed AlphaJudge summaries; figures are written to
`7_manuscript/` and are not tracked. `make_pdf.py` also needs Chromium.
Re-running AlphaJudge itself needs the structures from Zenodo.

## Licence and terms

AlphaFold 3 and AlphaFold Server outputs are subject to Google DeepMind's output
terms of use (non-commercial); see `8_zenodo/README.md`.
