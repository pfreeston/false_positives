#!/usr/bin/env python3
"""Typeset every table and figure as one print-ready A4 PDF.

No LaTeX is installed here, so the document is laid out in HTML with print CSS and
rendered by headless Chromium. Table data comes from the same loader as the .tex and
.md versions (make_iptm_table), so all three stay in step. Caption numbers (control
separation, correlation block, PCA) are read from metric_stats.json, written by
make_metric_figures.py, so captions follow the data.

The PRODIGY energetics table was dropped on 2026-11-10; it did not prove useful.

Layout: title + contents, then one table or figure per page, as journals normally ask.
Tall figures are fitted by height; the one wide figure (fig_metric_panel) gets a
landscape page so its labels stay legible.

Outputs: manuscript_tables_and_figures.html / .pdf
"""
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_iptm_table as IP        # noqa: E402
from controls import NEGATIVE, POSITIVE, control  # noqa: E402
from methods import AJ_METHODS, N_MODELS, PRETTY, SLUG  # noqa: E402
from make_derived_dotplots import DERIVED  # noqa: E402

MINUS = "−"
OUT = "manuscript_tables_and_figures"

CATEGORY_NAME = {"af_internal": "AlphaFold-internal", "physical": "physical and geometric",
                 "derived": "derived"}
CONTACT_PROB = ["interface_contact_prob_max", "interface_contact_prob_top10_mean",
                "interface_expected_contacts"]

CONTROLS_SENTENCE = (f"Positive controls (+), shaded blue, are interactions demonstrated "
                     f"in the lab; negative controls ({MINUS}), shaded amber, are those "
                     f"that could not be demonstrated — absence of evidence, not "
                     f"proven non-interaction. P69B_SSP1 is a positive control: SSP1 "
                     f"binds P69B, but not at the predicted interface, and the interaction "
                     f"does not affect P69B function.")
METHODS_SENTENCE = (
    "Methods: Boltz-2 (ColabFold MSA server) and Boltz-2 reduced-DB (given the AF3 "
    "pipeline MSAs); AF3 cluster (AF3's own MSA search, reduced BFD) and AF3 cluster "
    "full-DB (given the AFM full-database MSAs, BFD&nbsp;+&nbsp;UniRef30); AF3 server; "
    "AFM full-DB and AFM reduced-DB (AlphaFold-Multimer v3, 5 models &times; 5 seeds, "
    "of which the top five by AFM's own ranking are shown).")
AFM_NOTE = ("These are the 2026-11-09 cluster reruns, not the original screen; the five "
            "bars are the top five of 25 predictions (5 models &times; 5 seeds) by AFM's "
            "own ranking.")


def _num(n):
    return ["zero", "one", "two", "three", "four", "five", "six", "seven",
            "eight", "nine", "ten"][n]


def _list(xs):
    xs = list(xs)
    return xs[0] if len(xs) == 1 else ", ".join(xs[:-1]) + " and " + xs[-1]


def _sep(v):
    return f"{v:+.3f}".replace("-", MINUS)


def figure_list():
    st = json.loads((HERE / "metric_stats.json").read_text())
    a = pd.read_csv(HERE / "alphajudge_mature_by_model.csv")
    no_cp = [PRETTY[m] for m in AJ_METHODS
             if a.loc[a.method == m, CONTACT_PROB].isna().all().all()]
    with_cp = [PRETTY[m] for m in AJ_METHODS if PRETTY[m] not in no_cp]
    cats = st["categories"]
    cat_txt = "; ".join(
        f"{name.lower()} {c['separating']}/{c['n']}"
        for name, c in cats.items())
    sep_txt = (_list(f"{r['metric']} ({_sep(r['separation'])})" for r in st["separating"])
               if st["separating"] else "none")
    n_models = st["n_models"]

    figs = [
        ("fig_iptm_ptm.svg", "portrait",
         "Best ipTM&nbsp;+&nbsp;pTM for every model, by complex and method.",
         f"Filled markers are each method's top-ranked model, open markers the other four. "
         f"AFM appears three times: the original screen (one prediction per model, from the "
         f"screen record), and the AFM full-DB and reduced-DB cluster runs (top five of 25 "
         f"predictions by AFM's own ranking). The grey band spans the positive controls' "
         f"best values in the AFM screen. Positive controls are listed first, and each "
         f"group is ordered by the AFM screen value. {METHODS_SENTENCE} "
         f"{CONTROLS_SENTENCE}"),
        ("fig_metric_separation.svg", "portrait",
         f"Control separation for {st['n_metrics']} AlphaJudge metrics, grouped by category.",
         f"Separation is the lowest positive control minus the highest negative control "
         f"(inverted where lower values indicate a real interaction), divided by the "
         f"metric's range across the ten complexes; positive means the controls are "
         f"correctly ordered. Each complex's value is the mean over the seven methods of "
         f"each method's top-ranked model. Metrics that separate, by category: {cat_txt}. "
         f"Separating metrics: {sep_txt}. {st['blocking_negative']} is the highest-scoring "
         f"negative in {st['blocking_count']} of {st['blocking_of']} metrics. With "
         f"{_num(len(POSITIVE))} positive and {_num(len(NEGATIVE))} negative controls these values are descriptive, not "
         f"validated discriminative power. {METHODS_SENTENCE}"),
        ("fig_metric_panel.svg", "landscape",
         "The five metrics that survive clustering, one panel each.",
         f"Across the {n_models} models, ipSAE correlates |&rho;|&nbsp;&ge;&nbsp;"
         f"{st['ipsae_block_min_abs_rho']:.2f} (Spearman) with {_list(st['ipsae_block'])}; "
         f"over {st['pca_n_metrics']} numeric metrics, {st['pcs_for_80']} principal "
         f"components carry 80% of the variance and {st['pcs_for_90']} carry 90% "
         f"(PC1 alone {st['pc1']}%). Panel titles are coloured by source category. "
         f"Open markers are individual models, filled markers the highest value per method. "
         f"{METHODS_SENTENCE}"),
    ]
    for cat, name in CATEGORY_NAME.items():
        figs.append((f"fig_metrics_{cat}.svg", "portrait",
                     f"All {name} metrics, every method.",
                     f"Complexes on x; each metric in its own panel in native units. Open "
                     f"circles are the five models of each method, filled the top-ranked. "
                     + (f"{_list(no_cp)} supply no contact-probability metrics, so those "
                        f"panels show only {_list(with_cp)}. " if cat == "af_internal" else "")
                     + f"{METHODS_SENTENCE} {CONTROLS_SENTENCE}"))
    for col, slug, name in DERIVED:
        figs.append((f"fig_derived_{slug}.svg", "portrait",
                     f"{name} for every model, by complex and method.",
                     f"Same layout as Figure&nbsp;1. Open markers are individual models; the "
                     f"filled marker is each method's top-ranked model (by its own ipTM + pTM "
                     f"ranking), which is not always its highest {name} value. AFM full-DB and "
                     f"reduced-DB show the top five of 25 predictions; the AFM screen kept no "
                     f"result files, so AlphaJudge could not score it and it is absent. "
                     f"{METHODS_SENTENCE} {CONTROLS_SENTENCE}"))
    for meth in AJ_METHODS:
        for cat, cname in CATEGORY_NAME.items():
            body = (f"Five bars per complex, one per model, ordered left to right by this "
                    f"method's own ranking (rank&nbsp;1 solid, ranks&nbsp;2&ndash;5 lighter). "
                    f"The error bar is the mean of the five models &plusmn;&nbsp;1&nbsp;SD. "
                    f"Bars start at zero, which compresses metrics that never approach it; "
                    f"the dot-plot versions show those at full resolution. "
                    f"{CONTROLS_SENTENCE}")
            if meth.startswith("AFM"):
                body += " " + AFM_NOTE
            figs.append((f"fig_bars_{SLUG[meth]}_{cat}.svg", "portrait",
                         f"{PRETTY[meth]}: all {cname} metrics.", body))
    return figs


def num(v, dp=3, dagger=False):
    if pd.isna(v):
        return "&mdash;"
    s = f"{v:.{dp}f}".replace("-", MINUS)
    return s + "<sup>&dagger;</sup>" if dagger else s


def row_class(kind):
    return {"positive": ' class="pos"', "negative": ' class="neg"'}.get(kind, "")


def label(cx, kind, note=""):
    mark = " (+)" if kind == "positive" else f" ({MINUS})" if kind == "negative" else ""
    name = f"{cx}{mark}"
    if kind == "positive":
        name = f"<strong>{name}</strong>"
    return name + note


def table_iptm(n):
    wide, ctrl = IP.load()
    head = "".join(f"<th>{c}</th>" for c, _ in IP.COLUMNS)
    rows = []
    for cx in wide.index:
        k = ctrl[cx]
        note = '<sup>a</sup>' if cx == "P69B_Ecp36" else ""
        cells = "".join(f"<td>{num(wide.loc[cx, c], 3)}</td>" for _, c in IP.COLUMNS)
        rows.append(f'<tr{row_class(k)}><th scope="row">{label(cx, k, note)}</th>{cells}</tr>')
    return f"""
<section>
<p class="cap"><strong>Table {n}.</strong> Best ipTM&nbsp;+&nbsp;pTM per complex and method.
Each value is the highest-scoring of the five models produced by that method. AFM screen
values are the original screen; AFM full-DB and AFM reduced-DB are the AlphaFold-Multimer v3
cluster runs with full and reduced databases (top five of 25 predictions by AFM's own
ranking); all other values are the mature-construct runs. AF3 Cluster
uses AF3's own MSA search (reduced BFD) and AF3 Cluster full-DB is given the AFM
full-database MSAs (BFD&nbsp;+&nbsp;UniRef30); Boltz-2 uses the ColabFold MSA server and
Boltz-2 reduced-DB is given the AF3 pipeline MSAs. {CONTROLS_SENTENCE}</p>
<table>
  <thead>
    <tr><td class="corner"></td><th colspan="{len(IP.COLUMNS)}" class="spanner">Best ipTM + pTM</th></tr>
    <tr><th scope="col" class="stub">Complex</th>{head}</tr>
  </thead>
  <tbody>{''.join(rows)}</tbody>
</table>
<p class="note"><sup>a</sup> The screen run behind P69B_Ecp36's AFM screen value is not in the
repository; the only Ecp36 run on disk is FLAG-tagged and scores 0.537.</p>
</section>"""


CSS = """
@page { size: A4; margin: 18mm 16mm; }
html { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
body { font-family: "Source Serif 4","Charter","Georgia","Times New Roman",serif;
       font-size: 9.6pt; line-height: 1.42; color: #111; margin: 0; }
h1 { font-size: 13pt; font-weight: 600; margin: 0 0 2mm; letter-spacing: .01em; }
.sub { font-size: 8.4pt; color: #555; margin: 0 0 6mm; }
section { page-break-inside: avoid; margin-bottom: 8mm; }
section + section { page-break-before: always; }
.cap { font-size: 8.4pt; line-height: 1.38; margin: 0 0 2.2mm; text-align: justify; }
.note { font-size: 7.8pt; color: #444; margin: 1.8mm 0 0; line-height: 1.35; }
table { border-collapse: collapse; width: 100%;
        font-variant-numeric: tabular-nums; font-feature-settings: "tnum" 1; }
th, td { padding: 1.5mm 1.8mm; text-align: right; white-space: nowrap; }
thead th[scope="col"]:not(.stub) { white-space: normal; vertical-align: bottom;
        line-height: 1.2; }
th.stub { text-align: left; white-space: nowrap; }
th[scope="row"] { text-align: left; white-space: nowrap; font-weight: 400; }
thead tr:first-child .spanner { text-align: center; font-weight: 600;
        border-bottom: .4pt solid #333; padding-bottom: .8mm; }
thead tr:first-child td { border: 0; }
thead tr:last-child th { border-bottom: .5pt solid #333; font-weight: 600; }
table { border-top: 1.1pt solid #111; }
tbody tr:last-child td, tbody tr:last-child th { border-bottom: 1.1pt solid #111; }
tr.pos td, tr.pos th { background: #DCE9F5; }
tr.neg td, tr.neg th { background: #FBE6D4; }
sup { font-size: .72em; }
/* contents */
.toc { font-size: 8.2pt; border-collapse: collapse; width: 100%; margin-top: 1mm; }
.toc td { padding: .7mm 2mm .7mm 0; text-align: left; white-space: normal; }
.toc td.k { white-space: nowrap; font-weight: 600; width: 20mm; }
/* figures: fitted by height so tall panels are never clipped */
/* 222mm is the largest that keeps every caption+figure on a single page;
   236mm overflows three of the tall ones onto a second page. */
section.fig img { display: block; max-width: 100%; max-height: 222mm;
                  width: auto; height: auto; margin: 1.5mm auto 0; }
/* the one wide figure gets a real landscape page (a CSS-rotated block is paginated
   by its unrotated box and spilled onto the next page) */
@page wide { size: A4 landscape; margin: 14mm 16mm; }
section.rot { page: wide; }
section.rot img { display: block; max-width: 100%; max-height: 150mm;
                  width: auto; height: auto; margin: 1.5mm auto 0; }
"""


def build():
    figs = figure_list()
    missing = [f for f, *_ in figs if not (HERE / f).exists()]
    if missing:
        raise SystemExit(f"missing figures: {missing}")

    toc = [("Table 1", "Best ipTM + pTM per complex and method")]
    body = [table_iptm(1)]
    for i, (name, orient, lead, rest) in enumerate(figs, start=1):
        toc.append((f"Figure {i}", lead.rstrip(".")))
        cap = f'<p class="cap"><strong>Figure {i}.</strong> {lead} {rest}</p>'
        if orient == "landscape":
            body.append(f'<section class="rot">{cap}'
                        f'<img src="{name}" alt="Figure {i}"></section>')
        else:
            body.append(f'<section class="fig">{cap}'
                        f'<img src="{name}" alt="Figure {i}"></section>')
    rows = "".join(f'<tr><td class="k">{k}</td><td>{v}</td></tr>' for k, v in toc)
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>False positives &mdash; tables and figures</title>
<style>{CSS}</style></head>
<body>
<h1>Predicted hydrolase&ndash;effector complexes: prediction confidence and interface metrics</h1>
<p class="sub">Mature constructs (no prodomain, no signal peptide), seven prediction
methods including full- and reduced-database MSA variants. Generated from 7_manuscript;
see README.md for provenance and caveats.</p>
<table class="toc">{rows}</table>
{''.join(body)}
</body></html>"""
    return html, len(toc)


def render(src, pdf):
    exe = next((p for p in ("chromium-browser", "chromium", "google-chrome",
                            "google-chrome-stable") if shutil.which(p)), None)
    if exe is None:
        raise SystemExit("no Chromium/Chrome found to render the PDF")
    last = ""
    for flag in ("--headless=new", "--headless"):
        r = subprocess.run([exe, flag, "--disable-gpu", "--no-sandbox",
                            "--no-pdf-header-footer", f"--print-to-pdf={pdf}",
                            "--virtual-time-budget=20000", src.as_uri()],
                           capture_output=True, text=True, timeout=300)
        if pdf.exists() and pdf.stat().st_size > 2000:
            return exe, flag
        last = r.stderr.strip()[-600:]
    raise SystemExit(f"{exe} produced no usable PDF.\n{last}")


def main():
    html, n_items = build()
    src = HERE / f"{OUT}.html"
    pdf = HERE / f"{OUT}.pdf"
    src.write_text(html, encoding="utf-8")
    if pdf.exists():
        pdf.unlink()
    exe, flag = render(src, pdf)
    data = pdf.read_bytes()
    pages = data.count(b"/Type /Page") - data.count(b"/Type /Pages")
    print(f"rendered with {exe} {flag}")
    print(f"{pdf.name}: {len(data)/1024/1024:.1f} MB | {pages} pages | {n_items} items")


if __name__ == "__main__":
    main()
