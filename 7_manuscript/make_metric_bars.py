#!/usr/bin/env python3
"""Bar versions of the per-metric grids, one figure per method per category.

Complexes on x, one panel row per metric in native units. Within each complex
there are five bars -- one per model, ordered left to right by that method's own
ranking (rank 1 = highest ipTM+pTM). Rank 1 is solid, ranks 2-5 lighter.

A single model has no spread of its own, so the error bar is drawn at the
complex level: a dash at the mean of that method's five models with whiskers at
+/- 1 SD.

Every figure keeps the same ten x positions in the same order, so figures can be
laid side by side and read across methods.

Methods are methods.AJ_METHODS (seven), so there are 21 figures. AFM full-DB and
reduced-DB are the 2026-11-09 cluster reruns (top 5 of 25 by AFM's own ranking),
not the original screen, which kept no result pickles.

Bars start at zero in every panel. For metrics whose values never approach zero
(interface pLDDT, shape complementarity) that compresses the differences; the
dot-plot versions (fig_metrics_*) show those at full resolution.

Outputs: fig_bars_<method>_<category>.{pdf,png}
         with PER_METHOD = False instead: fig_bars_<category>.{pdf,png}
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"  # editable text in CorelDRAW/Illustrator:
matplotlib.rcParams["pdf.fonttype"] = 42      # SVG text as text, PDF fonts as TrueType
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from controls import control  # noqa: E402
from make_metric_grid import CATEGORIES, NAMES  # noqa: E402
from methods import AJ_METHODS, COLOUR, N_MODELS, PRETTY, SLUG  # noqa: E402

HERE = Path(__file__).resolve().parent

PER_METHOD = True

SHADE = {"positive": "#DCE9F5", "negative": "#FBE6D4"}
AFM_NOTE = "top 5 of 25 predictions (5 models × 5 seeds) by AFM's own ranking"


def load():
    a = pd.read_csv(HERE / "alphajudge_mature_by_model.csv")
    a = a.copy()
    a["rank"] = (a.groupby(["complex", "method"]).iptm_ptm
                  .rank(ascending=False, method="first").astype(int))
    methods = [m for m in AJ_METHODS if m in set(a.method)]
    # every figure keeps all ten complexes so they can be read across
    order = pd.read_csv(HERE / "iptm_ptm_by_complex.csv").complex.drop_duplicates().tolist()
    return a, order, methods


def panel(ax, a, order, metric, methods, group_off, model_off, width):
    vals = a[metric]
    for i, cx in enumerate(order):
        kind = control(cx)
        if kind in SHADE:
            ax.axvspan(i - 0.47, i + 0.47, color=SHADE[kind], alpha=0.55, zorder=0)
        for j, m in enumerate(methods):
            sub = a[(a.complex == cx) & (a.method == m)].sort_values("rank")
            if sub.empty or sub[metric].isna().all():
                continue
            for _, row in sub.iterrows():
                if pd.isna(row[metric]):
                    continue
                x = i + group_off[j] + model_off[int(row["rank"]) - 1]
                ax.bar(x, row[metric], width=width, color=COLOUR[m],
                       alpha=1.0 if row["rank"] == 1 else 0.45, linewidth=0, zorder=2)
            v = sub[metric].dropna()
            if len(v) > 1:
                ax.errorbar(i + group_off[j], v.mean(), yerr=v.std(), fmt="_",
                            color="0.15", ms=11 if len(methods) == 1 else 9,
                            lw=1.0, capsize=2.4, capthick=1.0, zorder=5)
    finite = vals.dropna()
    lo, hi = (finite.min(), finite.max()) if len(finite) else (0.0, 0.0)
    pad = 0.06 * (hi - lo if hi > lo else abs(hi) or 1)
    ax.set_ylim(min(0, lo - pad), max(0, hi + pad))
    ax.axhline(0, color="0.5", lw=0.6, zorder=1)
    flat = vals.notna().sum() and vals.nunique() <= 1
    ax.set_title(NAMES.get(metric, metric) + ("  (no variance)" if flat else ""),
                 fontsize=8.5, loc="left", pad=3, fontweight="semibold")
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.tick_params(axis="y", labelsize=7, length=2)
    ax.grid(axis="y", color="0.92", lw=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlim(-0.6, len(order) - 0.4)


def draw(key, a, order, methods, stem, heading, note=""):
    title, colour, metrics = CATEGORIES[key]
    metrics = [m for m in metrics if m in a.columns]
    n = len(metrics)
    single = len(methods) == 1
    group_off = [0.0] if single else np.linspace(-0.295, 0.295, len(methods))
    model_off = np.linspace(-0.30, 0.30, N_MODELS) if single else np.linspace(-0.108, 0.108, N_MODELS)
    width = 0.128 if single else 0.048

    fig, axes = plt.subplots(n, 1, figsize=(8.6 if single else 13.5, 1.08 * n + 1.7),
                             sharex=True)
    axes = np.atleast_1d(axes)
    for ax, metric in zip(axes, metrics):
        panel(ax, a, order, metric, methods, group_off, model_off, width)

    labels = [c + {"positive": "  (+)", "negative": "  (−)"}.get(control(c), "")
              for c in order]
    axes[-1].set_xticks(range(len(order)))
    axes[-1].set_xticklabels(labels, fontsize=8.5, rotation=38, ha="right")
    for t, c in zip(axes[-1].get_xticklabels(), order):
        t.set_fontweight("bold" if control(c) == "positive" else "normal")

    handles = []
    if not single:
        handles += [plt.Rectangle((0, 0), 1, 1, color=COLOUR[m], label=PRETTY.get(m, m))
                    for m in methods]
    else:
        handles.append(plt.Rectangle((0, 0), 1, 1, color=COLOUR[methods[0]],
                                     label=PRETTY.get(methods[0], methods[0])))
    handles += [
        plt.Rectangle((0, 0), 1, 1, color="0.35", label="rank 1 (top-ranked model)"),
        plt.Rectangle((0, 0), 1, 1, color="0.35", alpha=0.45, label="ranks 2–5"),
        plt.Line2D([], [], marker="_", ls="", color="0.15", ms=9,
                   label="mean ± SD over the 5 models"),
        plt.Rectangle((0, 0), 1, 1, color=SHADE["positive"], label="positive control"),
        plt.Rectangle((0, 0), 1, 1, color=SHADE["negative"], label="negative control"),
    ]
    fig.legend(handles=handles, fontsize=8, frameon=False, ncol=len(handles),
               loc="upper center", bbox_to_anchor=(0.5, 1.0))
    fig.suptitle(heading, fontsize=11, fontweight="bold", color=colour,
                 x=0.008, ha="left", y=1.038)
    if note:
        fig.text(0.008, 1.012, note, fontsize=7.6, color="0.3", ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.962))
    for ext in ("pdf", "png", "svg"):
        fig.savefig(HERE / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return n


def main():
    a, order, methods = load()
    print(f"methods: {methods} | complexes: {len(order)}")
    if not PER_METHOD:
        for key in CATEGORIES:
            n = draw(key, a, order, methods, f"fig_bars_{key}", CATEGORIES[key][0])
            print(f"  fig_bars_{key}: {n} panels")
        return
    for m in methods:
        sub = a[a.method == m]
        covered = sorted(set(sub.complex))
        note = AFM_NOTE if m.startswith("AFM") else ""
        for key in CATEGORIES:
            title = CATEGORIES[key][0]
            stem = f"fig_bars_{SLUG[m]}_{key}"
            draw(key, sub, order, [m], stem, f"{PRETTY[m]} — {title}", note)
        print(f"  {PRETTY[m]}: 3 figures | {len(covered)}/{len(order)} complexes")


if __name__ == "__main__":
    main()
