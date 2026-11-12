#!/usr/bin/env python3
"""Two figures over the AlphaJudge metrics.

Fig A  fig_metric_separation  how well each metric separates the positive from the
                              negative controls, grouped by category.
Fig B  fig_metric_panel       the five cluster representatives, one panel each.

Methods are methods.AJ_METHODS (seven, five models each over 10 complexes): Boltz-2,
Boltz-2 reduced-DB, AF3 cluster, AF3 cluster full-DB, AF3 server, AFM full-DB and AFM
reduced-DB (AFM = top 5 of 25 by its own ranking). The AFM screen has no AlphaJudge
metrics and is only in fig_iptm_ptm.

Separation is (lowest positive - highest negative), or its inverse where lower values
mean a real interaction, divided by the metric's range across the 10 complexes.
Positive = the controls are cleanly ordered. With two positives and two negatives this
is DESCRIPTIVE ONLY -- it is not validated discriminative power.

The numbers quoted in the figure captions (separation counts, blocking negative,
ipSAE correlation block, PCA) are computed here and written to metric_stats.json,
which make_pdf.py reads, so captions cannot drift from the data.

Complex order is taken from iptm_ptm_by_complex.csv so every figure in this folder
lists the complexes in the same order.

Outputs: fig_metric_separation.{pdf,png}, fig_metric_panel.{pdf,png},
         metric_separation.csv, metric_stats.json
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"  # editable text in CorelDRAW/Illustrator:
matplotlib.rcParams["pdf.fonttype"] = 42      # SVG text as text, PDF fonts as TrueType
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from controls import POS, NEG, control  # noqa: E402
from methods import AJ_METHODS, COLOUR, PRETTY  # noqa: E402

HERE = Path(__file__).resolve().parent

METHODS = AJ_METHODS

CATEGORIES = {
    "AlphaFold internal": ["iptm", "ptm", "iptm_ptm", "confidence_score",
                           "average_interface_pae", "interface_average_plddt",
                           "interface_contact_prob_max",
                           "interface_contact_prob_top10_mean",
                           "interface_expected_contacts"],
    "Physical & geometric": ["interface_area", "interface_contact_pairs",
                             "interface_num_intf_residues", "interface_hb",
                             "interface_sb", "interface_ss", "interface_sc",
                             "interface_polar", "interface_hydrophobic",
                             "interface_charged", "interface_solv_en"],
    "Derived scores": ["pDockQ/mpDockQ", "interface_pDockQ2", "interface_ipSAE",
                       "interface_LIS", "interface_cLIS", "interface_iLIS",
                       "interface_ccc", "interface_score", "interface_meta_score"],
}
CAT_COLOUR = {"AlphaFold internal": "#4C72B0", "Physical & geometric": "#C44E52",
              "Derived scores": "#55A868"}

# one representative per correlation cluster (see README)
REPRESENTATIVES = [("confidence_score", "Confidence score", "AlphaFold internal"),
                   ("interface_ipSAE", "ipSAE", "Derived scores"),
                   ("interface_ccc", "CCC", "Derived scores"),
                   ("pDockQ/mpDockQ", "pDockQ / mpDockQ", "Derived scores"),
                   ("interface_sc", "Shape complementarity", "Physical & geometric")]


NAMES = {
    "iptm": "ipTM", "ptm": "pTM", "iptm_ptm": "ipTM + pTM",
    "confidence_score": "Confidence score", "average_interface_pae": "Interface PAE",
    "interface_average_plddt": "Interface pLDDT",
    "interface_contact_prob_max": "Contact prob. (max)",
    "interface_contact_prob_top10_mean": "Contact prob. (top-10 mean)",
    "interface_expected_contacts": "Expected contacts",
    "interface_area": "Buried area", "interface_contact_pairs": "Contact pairs",
    "interface_num_intf_residues": "Interface residues", "interface_hb": "H-bonds",
    "interface_sb": "Salt bridges", "interface_ss": "S–S bonds",
    "interface_sc": "Shape complementarity", "interface_polar": "Polar fraction",
    "interface_hydrophobic": "Hydrophobic fraction",
    "interface_charged": "Charged fraction", "interface_solv_en": "Solvation energy",
    "pDockQ/mpDockQ": "pDockQ / mpDockQ", "interface_pDockQ2": "pDockQ2",
    "interface_ipSAE": "ipSAE", "interface_LIS": "LIS", "interface_cLIS": "cLIS",
    "interface_iLIS": "iLIS", "interface_ccc": "CCC",
    "interface_score": "Interface score", "interface_meta_score": "Meta score",
}


def pretty(m):
    return NAMES.get(m, m.replace("interface_", "").replace("_", " "))


def load():
    a = pd.read_csv(HERE / "alphajudge_mature_by_model.csv")
    a = a[a.method.isin(METHODS)].copy()
    a["method"] = pd.Categorical(a.method, METHODS, ordered=True)
    order = pd.read_csv(HERE / "iptm_ptm_by_complex.csv").complex.drop_duplicates().tolist()
    return a, order


def separation(best_per_complex, metric):
    v = best_per_complex[metric].dropna()
    rng = v.max() - v.min()
    if not rng or np.isnan(rng):
        return 0.0, "flat"
    p, n = v[POS], v[NEG]
    hi, lo = p.min() - n.max(), n.min() - p.max()
    return (hi / rng, "higher") if hi >= lo else (lo / rng, "lower")


# ------------------------------------------------------------------ figure A
def fig_separation(sep):
    rows, y, ticks, labels = [], 0.0, [], []
    for cat, metrics in CATEGORIES.items():
        sub = sep[sep.category == cat].sort_values("separation", ascending=False)
        y -= 1.0
        rows.append(("header", cat, y))
        for r in sub.itertuples():
            y -= 1.0
            rows.append(("metric", r, y))
            ticks.append(y)
            suffix = (" ↓" if r.direction == "lower"
                      else "  (no variance)" if r.direction == "flat" else "")
            labels.append(pretty(r.metric) + suffix)
    fig, ax = plt.subplots(figsize=(7.1, 8.4))
    ax.axvline(0, color="0.25", lw=1.1, zorder=3)
    for kind, obj, yy in rows:
        if kind == "header":
            ax.text(-0.76, yy, obj, fontsize=9.5, fontweight="bold",
                    color=CAT_COLOUR[obj], va="center")
        else:
            c = CAT_COLOUR[obj.category]
            ax.hlines(yy, 0, obj.separation, color=c, lw=2.2, alpha=.85, zorder=2)
            ax.plot(obj.separation, yy, "o", ms=6.5, color=c, zorder=4)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.78, 0.78)
    ax.set_ylim(min(t for t in ticks) - 1.0, 0.2)
    ax.set_xlabel("control separation   (negative ← fails | separates → positive)"
                  "\n↓ lower values indicate a real interaction", fontsize=9.5)
    ax.grid(axis="x", color="0.92", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ("pdf", "png", "svg"):
        fig.savefig(HERE / f"fig_metric_separation.{ext}", dpi=300, bbox_inches="tight")


# ------------------------------------------------------------------ figure B
def fig_panel(a, order):
    n = len(REPRESENTATIVES)
    fig, axes = plt.subplots(1, n, figsize=(12.0, 6.4), sharey=True)
    off = dict(zip(METHODS, np.linspace(0.33, -0.33, len(METHODS))))
    for ax, (metric, label, cat) in zip(axes, REPRESENTATIVES):
        for i, cx in enumerate(order):
            y0 = len(order) - 1 - i
            if control(cx) != "unlabelled":
                ax.axhspan(y0 - 0.45, y0 + 0.45, color="0.93", zorder=0)
            for m in METHODS:
                sub = a[(a.complex == cx) & (a.method == m)]
                if sub.empty or sub[metric].isna().all():
                    continue
                yy = y0 + off[m]
                ax.plot(sub[metric], [yy] * len(sub), "o", ms=3.0, mfc="none",
                        mec=COLOUR[m], mew=0.9, zorder=2)
                ax.plot([sub[metric].max()], [yy], "o", ms=4.8, color=COLOUR[m], zorder=3)
        ax.set_title(label, fontsize=9.5, color=CAT_COLOUR[cat], fontweight="bold", pad=7)
        ax.grid(axis="x", color="0.92", lw=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(axis="both", labelsize=8, length=0)
    labels = [f"{c}  (+)" if control(c) == "positive"
              else f"{c}  (−)" if control(c) == "negative" else c for c in order]
    axes[0].set_yticks(range(len(order) - 1, -1, -1))
    axes[0].set_yticklabels(labels, fontsize=8.5)
    for t, c in zip(axes[0].get_yticklabels(), order):
        if control(c) != "unlabelled":
            t.set_fontweight("bold")
    axes[0].set_ylim(-0.6, len(order) - 0.4)
    handles = [plt.Line2D([], [], marker="o", ls="", color=COLOUR[m], ms=6, label=PRETTY[m])
               for m in METHODS]
    handles.append(plt.Line2D([], [], marker="o", ls="", mfc="none", mec="0.35", ms=3.0,
                              label="individual models"))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.legend(handles=handles, fontsize=8, frameon=False, ncol=4,
               loc="lower center", bbox_to_anchor=(0.5, 0.0))
    for ext in ("pdf", "png", "svg"):
        fig.savefig(HERE / f"fig_metric_panel.{ext}", dpi=300, bbox_inches="tight")


def blocking(best_per_complex, sep):
    """Per separating-or-failing metric, the negative that scores most like a positive."""
    counts = {}
    for r in sep.itertuples():
        if r.direction == "flat":
            continue
        v = best_per_complex[r.metric].dropna()[NEG]
        cx = v.idxmax() if r.direction == "higher" else v.idxmin()
        counts[cx] = counts.get(cx, 0) + 1
    return counts


def correlation_and_pca(a):
    metrics = [m for ms in CATEGORIES.values() for m in ms]
    num = a[metrics].apply(pd.to_numeric, errors="coerce")
    rho = num.corr(method="spearman")["interface_ipSAE"].drop("interface_ipSAE")
    block = rho[rho.abs() >= 0.9].sort_values(key=abs, ascending=False)
    complete = num.dropna(axis=1, how="any")
    complete = complete.loc[:, complete.std() > 0]
    z = (complete - complete.mean()) / complete.std()
    ev = np.sort(np.linalg.eigvalsh(np.corrcoef(z.values, rowvar=False)))[::-1]
    cum = np.cumsum(ev) / ev.sum()
    return dict(
        n_models=int(len(a)),
        ipsae_block=[pretty(m) for m in block.index],
        ipsae_block_min_abs_rho=round(float(block.abs().min()), 2),
        ipsae_block_median_abs_rho=round(float(block.abs().median()), 2),
        pca_n_metrics=int(complete.shape[1]),
        pca_dropped=[pretty(m) for m in metrics if m not in complete.columns],
        pc1=round(float(cum[0]) * 100),
        pcs_for_80=int(np.searchsorted(cum, 0.80) + 1),
        pcs_for_90=int(np.searchsorted(cum, 0.90) + 1),
    )


def main():
    a, order = load()
    best = (a.sort_values("iptm_ptm", ascending=False)
             .groupby(["complex", "method"], as_index=False, observed=True).first())
    g = best.groupby("complex").mean(numeric_only=True)

    rows = []
    for cat, metrics in CATEGORIES.items():
        for m in metrics:
            if m not in g.columns:
                raise SystemExit(f"metric not in the AlphaJudge file: {m}")
            s, d = separation(g, m)
            rows.append(dict(category=cat, metric=m, direction=d, separation=round(s, 3)))
    sep = pd.DataFrame(rows)
    sep.sort_values(["category", "separation"], ascending=[True, False]).to_csv(
        HERE / "metric_separation.csv", index=False)

    fig_separation(sep)
    fig_panel(a, order)

    cats = {}
    for cat, sub in sep.groupby("category", sort=False):
        top = sub.loc[sub.separation.idxmax()]
        cats[cat] = dict(n=int(len(sub)), separating=int((sub.separation > 0).sum()),
                         best_metric=pretty(top.metric), best=float(top.separation),
                         median=round(float(sub.separation.median()), 3))
    blk = blocking(g, sep)
    worst = max(blk, key=blk.get)
    stats = dict(
        n_metrics=int(len(sep)),
        n_separating=int((sep.separation > 0).sum()),
        separating=[dict(metric=pretty(r.metric), separation=float(r.separation))
                    for r in sep[sep.separation > 0].sort_values("separation", ascending=False).itertuples()],
        flat=[pretty(m) for m in sep[sep.direction == "flat"].metric],
        categories=cats,
        blocking_negative=worst, blocking_count=int(blk[worst]),
        blocking_of=int((sep.direction != "flat").sum()),
        methods=METHODS,
        **correlation_and_pca(a),
    )
    (HERE / "metric_stats.json").write_text(json.dumps(stats, indent=1))

    print(sep.groupby("category").separation
             .agg(n="size", separating=lambda x: int((x > 0).sum()),
                  best="max", median="median").round(3).to_string())
    print("\ntop 6 overall:")
    print(sep.nlargest(6, "separation")[["metric", "category", "separation"]].to_string(index=False))
    print(f"\nblocking negative: {worst} in {blk[worst]}/{stats['blocking_of']} metrics  {blk}")
    print(f"ipSAE block (|rho|>=0.9): {stats['ipsae_block']}  min {stats['ipsae_block_min_abs_rho']}")
    print(f"PCA over {stats['pca_n_metrics']} metrics: PC1 {stats['pc1']}%, "
          f"80% at {stats['pcs_for_80']} PCs, 90% at {stats['pcs_for_90']}; dropped {stats['pca_dropped']}")
    print(f"\nmodels plotted: {len(a)} | complexes {a.complex.nunique()} | methods {METHODS}")


if __name__ == "__main__":
    main()
