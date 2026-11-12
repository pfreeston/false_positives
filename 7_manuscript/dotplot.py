#!/usr/bin/env python3
"""Shared layout for the per-complex dot plots (Figure 1 and the derived-score figures).

One row per complex (controls marked and bolded), one offset line per method inside
each row. Every model is an open marker; the method's chosen model (its top-ranked
model) is filled. An optional grey band marks a reference range on x.
"""
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"  # editable text in CorelDRAW/Illustrator:
matplotlib.rcParams["pdf.fonttype"] = 42      # SVG text as text, PDF fonts as TrueType
import matplotlib.pyplot as plt
import numpy as np

from controls import control
from methods import COLOUR, PRETTY


def dot_plot(d, value, methods, order, xlabel, out_stem, *, best_col=None,
             band=None, band_label=None, xlim=None, title=None):
    """d: long table with complex, method and `value`; best_col marks the filled model
    (default: the maximum per complex and method). Writes out_stem.{pdf,png}."""
    fig, ax = plt.subplots(figsize=(7.1, 0.51 * len(order) * len(methods) / 8 + 3.5))
    if band is not None:
        ax.axvspan(*band, color="0.88", zorder=0)
    off = dict(zip(methods, np.linspace(0.36, -0.36, len(methods))))
    for i, cx in enumerate(order):
        y0 = len(order) - 1 - i
        for m in methods:
            sub = d[(d.complex == cx) & (d.method == m)].dropna(subset=[value])
            if sub.empty:
                continue
            y = y0 + off[m]
            ax.plot(sub[value], [y] * len(sub), "o", ms=3.2, mfc="none",
                    mec=COLOUR[m], mew=1.0, zorder=2)
            best = sub[sub[best_col]][value] if best_col else [sub[value].max()]
            ax.plot(list(best), [y] * len(best), "o", ms=5.0, color=COLOUR[m], zorder=3)
    labels = [f"{c}  (+)" if control(c) == "positive"
              else f"{c}  (−)" if control(c) == "negative" else c for c in order]
    ax.set_yticks(range(len(order) - 1, -1, -1))
    ax.set_yticklabels(labels, fontsize=9)
    for t, c in zip(ax.get_yticklabels(), order):
        if control(c) != "unlabelled":
            t.set_fontweight("bold")
    ax.set_xlabel(xlabel, fontsize=10)
    if xlim is not None:
        ax.set_xlim(*xlim)
    ax.set_ylim(-0.6, len(order) - 0.4)
    for k in range(1, len(order)):
        ax.axhline(k - 0.5, color="0.9", lw=0.5, zorder=0)
    ax.grid(axis="x", color="0.9", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    handles = [plt.Line2D([], [], marker="o", ls="", color=COLOUR[m], ms=6, label=PRETTY[m])
               for m in methods]
    handles.append(plt.Line2D([], [], marker="o", ls="", mfc="none", mec="0.35", ms=4.2,
                              label="individual models"))
    if band is not None:
        handles.append(plt.Rectangle((0, 0), 1, 1, color="0.88", label=band_label))
    leg = ax.legend(handles=handles, fontsize=8, frameon=False, ncol=4,
                    loc="lower center", bbox_to_anchor=(0.42, 1.0))
    if title:
        fig.suptitle(title, fontsize=11, fontweight="bold", x=0.02, ha="left")
    fig.tight_layout()
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{out_stem}.{ext}", dpi=300, bbox_inches="tight",
                    bbox_extra_artists=[leg])
    plt.close(fig)
