#!/usr/bin/env python3
"""One figure per metric category: complexes on x, every metric as its own panel row.

Each panel is one metric in its native units. Within a panel, each complex has one
cluster per method (methods.AJ_METHODS, seven), and every one of that method's five
models is drawn -- open
circles for all five, filled for the method's top-ranked model (highest ipTM+pTM,
the same definition used by the other figures here). Control complexes are shaded.

AFM here is the 2026-11-09 cluster rerun (full-DB and reduced-DB, top 5 of 25 models
by AFM's own ranking), not the original screen, which kept no result pickles.

Outputs: fig_metrics_af_internal, fig_metrics_physical, fig_metrics_derived (.pdf/.png)
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"  # editable text in CorelDRAW/Illustrator:
matplotlib.rcParams["pdf.fonttype"] = 42      # SVG text as text, PDF fonts as TrueType
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from controls import POS, NEG, control  # noqa: E402
from methods import AJ_METHODS, COLOUR, PRETTY  # noqa: E402

HERE = Path(__file__).resolve().parent

CATEGORIES = {
    "af_internal": ("AlphaFold internal", "#4C72B0", [
        "iptm", "ptm", "iptm_ptm", "confidence_score", "average_interface_pae",
        "interface_average_plddt", "interface_contact_prob_max",
        "interface_contact_prob_top10_mean", "interface_expected_contacts"]),
    "physical": ("Physical & geometric", "#C44E52", [
        "interface_area", "interface_contact_pairs", "interface_num_intf_residues",
        "interface_hb", "interface_sb", "interface_ss", "interface_sc",
        "interface_polar", "interface_hydrophobic", "interface_charged",
        "interface_solv_en"]),
    "derived": ("Derived scores", "#55A868", [
        "pDockQ/mpDockQ", "interface_pDockQ2", "interface_ipSAE", "interface_LIS",
        "interface_cLIS", "interface_iLIS", "interface_ccc", "interface_score",
        "interface_meta_score"]),
}

NAMES = {
    "iptm": "ipTM", "ptm": "pTM", "iptm_ptm": "ipTM + pTM",
    "confidence_score": "Confidence score (backend-specific)",
    "average_interface_pae": "Interface PAE ↓",
    "interface_average_plddt": "Interface pLDDT",
    "interface_contact_prob_max": "Contact prob. (max)",
    "interface_contact_prob_top10_mean": "Contact prob. (top-10 mean)",
    "interface_expected_contacts": "Expected contacts",
    "interface_area": "Buried area (Å²)",
    "interface_contact_pairs": "Contact pairs",
    "interface_num_intf_residues": "Interface residues",
    "interface_hb": "H-bonds", "interface_sb": "Salt bridges",
    "interface_ss": "S–S bonds", "interface_sc": "Shape complementarity",
    "interface_polar": "Polar fraction", "interface_hydrophobic": "Hydrophobic fraction",
    "interface_charged": "Charged fraction",
    "interface_solv_en": "Solvation energy ↓",
    "pDockQ/mpDockQ": "pDockQ / mpDockQ", "interface_pDockQ2": "pDockQ2",
    "interface_ipSAE": "ipSAE", "interface_LIS": "LIS", "interface_cLIS": "cLIS",
    "interface_iLIS": "iLIS", "interface_ccc": "CCC (confident contact count)",
    "interface_score": "Interface score", "interface_meta_score": "Meta score",
}


def load():
    a = pd.read_csv(HERE / "alphajudge_mature_by_model.csv")
    methods = [m for m in AJ_METHODS if m in set(a.method)]
    a = a.copy()
    a["method"] = pd.Categorical(a.method, methods, ordered=True)
    order = [c for c in pd.read_csv(HERE / "iptm_ptm_by_complex.csv")
             .complex.drop_duplicates() if c in set(a.complex)]
    # the method's own top-ranked model, as used by every other figure here
    best = (a.sort_values("iptm_ptm", ascending=False)
             .groupby(["complex", "method"], observed=True).head(1).index)
    a["is_best"] = a.index.isin(best)
    return a, order, methods


def draw(key, a, order, methods):
    title, colour, metrics = CATEGORIES[key]
    metrics = [m for m in metrics if m in a.columns]
    n = len(metrics)
    fig, axes = plt.subplots(n, 1, figsize=(10.5, 1.02 * n + 1.9), sharex=True)
    axes = np.atleast_1d(axes)
    off = dict(zip(methods, np.linspace(-0.36, 0.36, len(methods))))
    for ax, metric in zip(axes, metrics):
        flat = a[metric].notna().sum() and a[metric].nunique() <= 1
        for i, cx in enumerate(order):
            if control(cx) != "unlabelled":
                ax.axvspan(i - 0.45, i + 0.45, color="0.93", zorder=0)
            for m in methods:
                sub = a[(a.complex == cx) & (a.method == m)]
                if sub.empty or sub[metric].isna().all():
                    continue
                x = i + off[m]
                ax.plot([x] * len(sub), sub[metric], "o", ms=3.0, mfc="none",
                        mec=COLOUR[m], mew=0.8, zorder=2)
                b = sub[sub.is_best]
                if not b.empty:
                    ax.plot([x], [b[metric].iloc[0]], "o", ms=4.8,
                            color=COLOUR[m], zorder=3)
        label = NAMES.get(metric, metric)
        if flat:
            label += "  (no variance)"
        ax.set_title(label, fontsize=8.5, loc="left", pad=3, fontweight="semibold")
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.tick_params(axis="y", labelsize=7, length=2)
        ax.grid(axis="y", color="0.93", lw=0.5)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.set_xlim(-0.6, len(order) - 0.4)
    labels = [f"{c}  (+)" if control(c) == "positive"
              else f"{c}  (−)" if control(c) == "negative" else c for c in order]
    axes[-1].set_xticks(range(len(order)))
    axes[-1].set_xticklabels(labels, fontsize=8, rotation=38, ha="right")
    for t, c in zip(axes[-1].get_xticklabels(), order):
        if control(c) != "unlabelled":
            t.set_fontweight("bold")
    handles = [plt.Line2D([], [], marker="o", ls="", color=COLOUR[m], ms=6,
                          label=PRETTY[m]) for m in methods]
    handles += [plt.Line2D([], [], marker="o", ls="", mfc="none", mec="0.35", ms=3.4,
                           label="each of the 5 models"),
                plt.Line2D([], [], marker="o", ls="", color="0.35", ms=5.6,
                           label="top-ranked model")]
    fig.legend(handles=handles, fontsize=8, frameon=False, ncol=5,
               loc="upper center", bbox_to_anchor=(0.5, 1.0 + 0.3 / fig.get_figheight()))
    fig.suptitle(title, fontsize=11, fontweight="bold", color=colour,
                 x=0.012, ha="left", y=1.0 + 0.55 / fig.get_figheight())
    fig.tight_layout()
    fig.subplots_adjust(top=1 - 0.3 / fig.get_figheight())
    for ext in ("pdf", "png", "svg"):
        fig.savefig(HERE / f"fig_metrics_{key}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return n


def main():
    a, order, methods = load()
    print(f"methods: {methods} | complexes: {len(order)} | models: {len(a)}")
    for key in CATEGORIES:
        n = draw(key, a, order, methods)
        print(f"  fig_metrics_{key}: {n} metric panels")
    miss = a.groupby("method", observed=True).apply(
        lambda d: [c for c in d.columns
                   if c in NAMES and d[c].isna().all()], include_groups=False)
    for m, cols in miss.items():
        if cols:
            print(f"  note: {m} has no values for {cols}")


if __name__ == "__main__":
    main()
