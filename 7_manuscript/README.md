# Manuscript tables and figures

Everything here is generated. Run the scripts in this order with the system `python3`
(needs pandas and matplotlib), then rebuild the PDF:

```
python3 collect_manuscript_csvs.py   # AlphaJudge summaries -> per-model / per-complex CSVs
python3 compare_iptm_ptm.py          # ipTM+pTM comparison incl. the AFM screen (Figure 1)
python3 make_iptm_table.py           # Table 1 (.tex, .md, standalone .tex)
python3 make_metric_figures.py       # Figures 2-3 + metric_stats.json
python3 make_metric_grid.py          # Figures 4-6 (all metrics, one figure per category)
python3 make_derived_dotplots.py     # Figure-1-style dot plot per derived score (9 figures)
python3 make_metric_bars.py          # per-method bar charts (21 figures)
python3 make_pdf.py                  # manuscript_tables_and_figures.{html,pdf}
```

`methods.py` is the single source of truth for the method set, colours, labels and the
number of models per method; `controls.py` is the same for the control labels. Do not
re-declare either in individual scripts.

## Methods

All runs use the mature constructs (no prodomain, no signal peptide; AFM-matched
sequences in `4_afm_cluster/fastas`). Seven AlphaJudge methods, five models each per
complex:

| Key | Label | Run | AlphaJudge summary |
|---|---|---|---|
| `Boltz2` | Boltz-2 | Boltz-2, ColabFold MSA server | `boltz2_matched_summary.csv` |
| `Boltz2_reduced` | Boltz-2 reduced-DB | Boltz-2 given the AF3 data-pipeline MSAs | `boltz2_reduced_summary.csv` |
| `AF3_cluster` | AF3 cluster | local AF3, own MSA search (reduced BFD) | `af3_matched_summary.csv` |
| `AF3_fulldb` | AF3 cluster full-DB | local AF3 given the AFM full-database MSAs | `af3_fulldb_summary.csv` |
| `AF3_server` | AF3 server | AlphaFold Server | `af3_server_summary.csv` |
| `AFM` | AFM full-DB | AlphaFold-Multimer v3, `full_dbs`, 5 models × 5 seeds | `afm_matched_summary.csv` |
| `AFM_reduced` | AFM reduced-DB | as above, `reduced_dbs` | `afm_reduced_summary.csv` |

Summaries are in `5_alpha_judge/`. **AFM produces 25 predictions per complex; only the
top five by AFM's own `ranking_confidence` (0.8·ipTM + 0.2·pTM) are used**, so every
method contributes five models. `collect_manuscript_csvs.py` asserts 10 complexes × 5
models for every method.

The structures themselves are not in git; they are deposited on Zenodo (`8_zenodo/`).

## Generated files

| File | Rows | Contents |
|---|---|---|
| `alphajudge_mature_by_model.csv` | 350 | every scored model, 7 methods, AlphaJudge columns + `method`, `complex`, `source_csv` |
| `alphajudge_mature_by_complex.csv` | 10 | top model per complex per method, with consensus tier |
| `iptm_ptm_by_model.csv` | 400 | 10 complexes × 8 methods (AFM screen + 7 AlphaJudge methods) × 5 models |
| `iptm_ptm_by_complex.csv` | 80 | best / mean / sd / spread per complex per method |
| `iptm_ptm_table.{tex,md}`, `_standalone.tex` | 10 | Table 1 |
| `metric_separation.csv` | 29 | separation value and direction per metric |
| `metric_stats.json` | — | numbers quoted in the figure captions (separation, ρ block, PCA) |
| `fig_*.{pdf,svg,png}` | — | figures (not tracked); PDF and SVG are vector with editable text, PNGs 300 dpi |
| `manuscript_tables_and_figures.pdf` | — | not tracked; A4, 38 pages: contents, Table 1, then one figure per page (36 figures) |

## Table 1 and Figure 1: ipTM + pTM

**Table 1 and Figure 1 show all three AFM runs**: the original screen, and the AFM
full-DB and reduced-DB cluster runs (top five of 25 by AFM's ranking), next to the five
other methods. Screen values come from the screen spreadsheet (`False_positives.xlsx`,
sheet "AFM Interaction Scores"), which matches the screen runs in `4_afm/` for nine of
ten complexes (see Ecp36 below). The spreadsheet is kept local and out of git;
`compare_iptm_ptm.py` writes the 50 per-model screen scores to `afm_screen_scores.csv`
(committed) whenever the spreadsheet is present and reads that CSV when it is not, so the
table and figure rebuild from a clean clone. Complexes are listed positive controls
first, each group ordered by the AFM screen value; every other figure uses this order.

**P69B_Ecp36 screen value.** The sheet gives 0.750 (model_4); the only Ecp36 AFM screen
run on disk is FLAG-tagged and scores 0.537. The run behind the sheet's row is not in the
repository. Decision (2026-11-06): the spreadsheet is authoritative, and the mismatch is
printed on every run and footnoted in Table 1.

Metric comparability: `iptm_ptm == 0.8·iptm + 0.2·ptm` exactly for AF3 and AFM. For
Boltz-2 some models deviate because its `iptm` is `chain_pair` scoped;
`iptm_ptm_recomputed` is carried alongside.

The .tex table needs `\usepackage{booktabs}` and `\usepackage[table]{xcolor}`. Not
compiled here (no LaTeX installed).

## Figures 2-3: control separation and cluster representatives

Separation = (lowest positive − highest negative), inverted where lower is better, over
the metric's range across the ten complexes. Each complex's value is the mean over the
seven methods of each method's top-ranked model.

| Category | Separating | Best | Median |
|---|---|---|---|
| AlphaFold internal | 7/9 | ipTM +0.263 | +0.152 |
| Derived scores | 8/9 | CCC +0.374 | +0.217 |
| Physical & geometric | 1/11 | shape complementarity +0.023 | −0.576 |

16 of 29 metrics separate the controls; the strongest are CCC (+0.374), ipSAE (+0.361),
ipTM (+0.263), meta score (+0.259) and confidence score (+0.250). Physical and geometric
metrics barely separate at all. **P69B_TIL is the highest-scoring negative in 16 of 28
metrics.** S–S bonds has no variance and is labelled as such.

These numbers changed sharply when P69B_SSP1 became a positive control (2026-11-10). With
SSP1 as a negative only confidence score, ipTM+pTM and ipTM separated (+0.07 to +0.12),
because SSP1 outscored the weakest positive on almost every metric.

Across the 350 models ipSAE correlates |ρ| ≥ 0.94 (Spearman; median 0.96) with LIS, iLIS,
cLIS, both contact-probability metrics, ipTM, meta score, pDockQ2, ipTM+pTM and interface
PAE. Over 26 numeric metrics, 4 principal components carry 80% of the variance and 6
carry 90% (PC1 52%). The Figure 3 panels are one representative per cluster.

These numbers are regenerated into `metric_stats.json` and the captions on every run;
the values above are from the 2026-11-10 build.

## Figures 4-6 and the bar charts

`make_metric_grid.py`: one figure per category (9, 11 and 9 panels), complexes on x, all
five models open with the top-ranked filled. **Boltz-2 and Boltz-2 reduced-DB supply no
contact-probability metrics** (`interface_contact_prob_*`, `interface_expected_contacts`
are NaN), so those panels show the AF3 and AFM methods only.

`make_derived_dotplots.py`: one figure per derived score (pDockQ/mpDockQ, pDockQ2, ipSAE,
LIS, cLIS, iLIS, CCC, interface score, meta score), named `fig_derived_<metric>`, in the
same layout as Figure 1 (`dotplot.py` holds the shared plotting code). The filled marker
is each method's top-ranked model by ipTM + pTM, which is not always its highest value
for that score. The AFM screen is absent because AlphaJudge could not score it. PDF
Figures 7-15.

`make_metric_bars.py`: one figure per method per category, 21 in all, named
`fig_bars_<slug>_<category>` (slugs: boltz2, boltz2_reduced, af3_cluster, af3_fulldb,
af3_server, afm, afm_reduced). Five bars per complex ordered by the method's own
ranking; the error bar is mean ± 1 SD of the five models. Bars start at zero, which
compresses metrics such as pLDDT; the dot plots show those at full resolution.

## PDF

`make_pdf.py` lays everything out in HTML with print CSS and renders it with headless
Chromium. Captions are built from `metric_stats.json` and the by-model CSV, so they stay
in step with the data. Figures embed the vector SVGs, so the whole PDF is
vector (no raster images) and its text is real, embedded TrueType text: it opens for
editing in CorelDRAW or Illustrator. For single figures use `fig_*.svg` (text stays
text, `svg.fonttype = none`) or `fig_*.pdf` (TrueType, `pdf.fonttype = 42`).
The wide `fig_metric_panel` gets its own landscape page.

## Not used

- **Energetics (PRODIGY)** was dropped from the manuscript (2026-11-10): it did not prove
  useful. Its scripts and outputs are kept locally and are not in the repository.
- **Prodomain runs** were deleted (2026-11-10).
- **Old AFM reruns** (`4_afm/from_cluster`, `4_afm/rerun_from_features`) are superseded
  by the AFM full-DB / reduced-DB runs and are not used.

## Caveats

- **Controls.** Positives: P69B_Epi1, Pip1_EpiC2B, Pip1_EpiC1 and P69B_SSP1. SSP1 was
  moved to the positives on 2026-11-10: it binds P69B, but not at the predicted
  interface, and the interaction does not affect P69B function. Negatives: the other
  six, meaning *could not be shown in the lab*, which is absence of evidence rather
  than demonstrated non-interaction.
- **Separation values are descriptive.** Four positives and six negatives cannot
  support an AUC or a fitted threshold.
- **AF3 server P69B length** is 634 aa (`SLHTT…`) in eight runs and 635 aa in the others;
  both are mature.
