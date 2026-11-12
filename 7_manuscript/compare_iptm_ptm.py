#!/usr/bin/env python3
"""Compare iptm+ptm across methods, for each complex.

Methods are methods.IPTM_METHODS: all three AFM runs (the original screen and the
full-DB and reduced-DB cluster runs, top five of 25 by AFM's ranking) plus Boltz-2,
Boltz-2 reduced-DB, AF3 cluster, AF3 cluster full-DB and AF3 server.

AFM comes from the SCREEN SPREADSHEET (False_positives.xlsx, sheet "AFM Interaction
Scores"), which is the record of the original screen: five per-model iptm+ptm values
and the best model per complex. The screen kept no result pickles, so iptm+ptm is all
AFM offers -- no raw iptm, no PAE. 4_afm/from_cluster and 4_afm/rerun_from_features
are NOT used.

The spreadsheet is cross-checked against 4_afm/*/ranking_debug.json on every run.
Nine complexes agree to 1e-3. P69B_Ecp36 does not: the sheet records
0.4747/0.5736/0.4813/0.7502/0.4809 (best 0.750, model_4) while the only Ecp36 run on
disk, 4_afm/Solyc08g079870.3.1_flagKAH3648627.1, gives 0.466-0.537 for a FLAG-TAGGED
construct (chain B 64 aa, DYKDDDDKGG+KRGSGG...). The sheet's Ecp36 run is not in the
repository. The sheet is used, and the mismatch is reported by `cross_check`.

Comparability of the metric:
  AF3_cluster, AF3_server  iptm_ptm == 0.8*iptm + 0.2*ptm exactly
  Boltz2                   19/50 models deviate (mean 0.021, max 0.099); its `iptm` is
                           chain_pair scoped, so `iptm_ptm_recomputed` is carried
                           alongside the reported value. The deviation is concentrated
                           in the weak P69B interfaces (Ecp36, EpiC2B, Six15, EpiC1).
  AFM, AFM_reduced         iptm_ptm == 0.8*iptm + 0.2*ptm (AFM's ranking_confidence)
  AFM_screen               spreadsheet value only; no iptm/ptm components available.

Outputs (this folder):
  iptm_ptm_by_model.csv    every model, 8 methods x 10 complexes
  iptm_ptm_by_complex.csv  best / mean / sd / spread per complex per method
  fig_iptm_ptm.pdf/.png    dot plot, all models, controls annotated
"""
import json
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from controls import POS, NEG, control  # noqa: E402
from methods import IPTM_METHODS  # noqa: E402
from dotplot import dot_plot  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
XLSX = HERE / "False_positives.xlsx"
SHEET = "AFM Interaction Scores"
# The spreadsheet is kept out of git (user decision 2026-11-10). Its 50 per-model screen
# scores are exported to this CSV, which IS committed; without the .xlsx the script
# reads the CSV instead, so Table 1 and Figure 1 rebuild from a clean clone.
SCREEN_CSV = HERE / "afm_screen_scores.csv"

METHODS = IPTM_METHODS

# screen run directory -> complex, used only to cross-check the spreadsheet
SCREEN_DIRS = {
    "matureP69B_matureEpi1": "P69B_Epi1",
    "matureP69B_maturePiEpiC1": "P69B_EpiC1",
    "matureP69B_maturePiEpiC2B": "P69B_EpiC2B",
    "maturePip1_matureEpi1": "Pip1_Epi1",
    "maturePip1_maturePiEpiC1": "Pip1_EpiC1",
    "maturePip1_maturePiEpiC2B": "Pip1_EpiC2B",
    "Solyc08g079870.3.1_APP91304.1": "P69B_Six15",
    "Solyc08g079870.3.1_WP046932418.1": "P69B_SSP1",
    "Solyc08g079870.3.1_XP018243121.1": "P69B_TIL",
    "Solyc08g079870.3.1_flagKAH3648627.1": "P69B_Ecp36",
}
NOTE = {"P69B_Ecp36": "spreadsheet run not in 4_afm; only on-disk Ecp36 run is FLAG-tagged"}


# ---------------------------------------------------------------- xlsx reader
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def _col(ref):
    n = 0
    for ch in re.match(r"([A-Z]+)", ref).group(1):
        n = n * 26 + ord(ch) - 64
    return n - 1


def read_sheet(path, want):
    """Minimal .xlsx reader -- avoids an openpyxl dependency."""
    z = zipfile.ZipFile(path)
    shared = ["".join(t.text or "" for t in si.iter(NS + "t"))
              for si in ET.fromstring(z.read("xl/sharedStrings.xml"))] \
        if "xl/sharedStrings.xml" in z.namelist() else []
    names = [s.get("name") for s in ET.fromstring(z.read("xl/workbook.xml")).iter(NS + "sheet")]
    if want not in names:
        raise SystemExit(f"sheet {want!r} not in {names}")
    rows = []
    for row in ET.fromstring(z.read(f"xl/worksheets/sheet{names.index(want) + 1}.xml")).iter(NS + "row"):
        cells = {}
        for c in row.iter(NS + "c"):
            v, t = c.find(NS + "v"), c.get("t")
            if t == "inlineStr":
                node = c.find(NS + "is")
                val = "".join(x.text or "" for x in node.iter(NS + "t")) if node is not None else None
            elif v is None:
                val = None
            elif t == "s":
                val = shared[int(v.text)]
            else:
                try:
                    val = float(v.text)
                except (TypeError, ValueError):
                    val = v.text
            cells[_col(c.get("r"))] = val
        rows.append([cells.get(j) for j in range(max(cells) + 1 if cells else 0)])
    return rows


def afm_screen():
    rows = read_sheet(XLSX, SHEET)
    header = [re.sub(r"\s+", " ", str(h or "")).strip() for h in rows[0]]
    model_cols = {}
    for i, h in enumerate(header):
        m = re.match(r"AFM model_(\d+)", h)
        if m:
            model_cols[f"model_{m.group(1)}"] = i
    if len(model_cols) != 5:
        raise SystemExit(f"expected 5 AFM model columns, found {sorted(model_cols)}")
    out = []
    for r in rows[1:]:
        if not r or not r[0]:
            continue
        cx = str(r[0]).strip()
        for model, i in model_cols.items():
            v = r[i] if i < len(r) else None
            if v is None:
                continue
            out.append(dict(complex=cx, method="AFM_screen", model=model,
                            iptm_ptm=round(float(v), 6), iptm_ptm_recomputed=np.nan,
                            iptm=np.nan, ptm=np.nan))
    return pd.DataFrame(out)


def afm_screen_scores():
    """Spreadsheet if present (and refresh the CSV from it), else the committed CSV."""
    if XLSX.exists():
        afm = afm_screen()
        afm[["complex", "model", "iptm_ptm"]].to_csv(SCREEN_CSV, index=False)
        print(f"read {XLSX.name}; wrote {SCREEN_CSV.name}")
        return afm
    if not SCREEN_CSV.exists():
        raise SystemExit(f"neither {XLSX.name} nor {SCREEN_CSV.name} found")
    d = pd.read_csv(SCREEN_CSV)
    d["method"] = "AFM_screen"
    for c in ("iptm_ptm_recomputed", "iptm", "ptm"):
        d[c] = np.nan
    print(f"{XLSX.name} absent; read {SCREEN_CSV.name}")
    return d[["complex", "method", "model", "iptm_ptm", "iptm_ptm_recomputed", "iptm", "ptm"]]


def cross_check(afm):
    """Compare the spreadsheet against ranking_debug.json; report, do not overwrite."""
    sheet_best = afm.groupby("complex").iptm_ptm.max()
    print("cross-check vs 4_afm/*/ranking_debug.json")
    for d, cx in SCREEN_DIRS.items():
        f = ROOT / "4_afm" / d / "ranking_debug.json"
        if not f.exists():
            print(f"  {cx:14s} no ranking_debug.json on disk")
            continue
        o = json.load(open(f))
        disk = max(o[next(k for k in o if k != "order")].values())
        s = sheet_best.get(cx, np.nan)
        flag = "OK" if abs(disk - s) < 1e-3 else f"MISMATCH (disk {disk:.3f})"
        print(f"  {cx:14s} sheet {s:.3f}  {flag}")


def others():
    a = pd.read_csv(HERE / "alphajudge_mature_by_model.csv")
    a = a[a.method.isin(METHODS)].copy()
    a["iptm_ptm_recomputed"] = (0.8 * a.iptm + 0.2 * a.ptm).round(4)
    return a.rename(columns={"model_used": "model"})[
        ["complex", "method", "model", "iptm_ptm", "iptm_ptm_recomputed", "iptm", "ptm"]]


def main():
    afm = afm_screen_scores()
    cross_check(afm)
    d = pd.concat([afm, others()], ignore_index=True)
    missing = set(d[d.method != "AFM_screen"].complex) - set(afm.complex)
    if missing:
        raise SystemExit(f"complexes absent from the spreadsheet: {sorted(missing)}")
    d["control"] = d.complex.map(control)
    d["construct_note"] = d.apply(
        lambda r: NOTE.get(r.complex, "") if r.method == "AFM_screen" else "", axis=1)
    d["method"] = pd.Categorical(d.method, METHODS, ordered=True)
    d = d.sort_values(["complex", "method", "iptm_ptm"], ascending=[True, True, False])
    d.to_csv(HERE / "iptm_ptm_by_model.csv", index=False)

    grp = d.groupby(["complex", "method"], observed=True).iptm_ptm
    g = grp.agg(best="max", mean="mean", sd="std", n="size").round(3).reset_index()
    g["spread"] = (grp.max() - grp.min()).round(3).values
    g["control"] = g.complex.map(control)
    wide = g.pivot(index="complex", columns="method", values="best")
    # positive controls first, then negatives; each group by AFM screen value
    order = sorted(wide.index, key=lambda c: (control(c) != "positive",
                                              -wide.loc[c, "AFM_screen"]))
    g.set_index("complex").loc[order].reset_index().to_csv(
        HERE / "iptm_ptm_by_complex.csv", index=False)

    # ---- figure ----
    dot_plot(d, "iptm_ptm", METHODS, order, "iptm + ptm", HERE / "fig_iptm_ptm",
             band=(min(wide.loc[c, "AFM_screen"] for c in POS),
                   max(wide.loc[c, "AFM_screen"] for c in POS)),
             band_label="positive-control range", xlim=(0.2, 1.0))

    print()
    print(wide.loc[order].round(3).to_string())
    print(f"\n{len(d)} model rows | {d.complex.nunique()} complexes | {len(METHODS)} methods")


if __name__ == "__main__":
    main()
