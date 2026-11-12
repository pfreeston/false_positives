#!/usr/bin/env python3
"""One Figure-1-style dot plot per AlphaJudge derived score.

Same layout as fig_iptm_ptm (dotplot.dot_plot): one row per complex in the shared
order (positive controls first), one offset line per AlphaJudge method, every model
open and the method's top-ranked model filled. The AFM screen is absent: it kept no
result pickles, so AlphaJudge could not score it.

Metrics are make_metric_grid.CATEGORIES["derived"]. Writes
fig_derived_<slug>.{pdf,png}; slugs are listed in DERIVED below.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from dotplot import dot_plot  # noqa: E402
from make_metric_grid import CATEGORIES, NAMES, load  # noqa: E402
from methods import IPTM_METHODS  # noqa: E402

# column -> file slug, in the order of make_metric_grid's derived category
SLUGS = {"pDockQ/mpDockQ": "pdockq", "interface_pDockQ2": "pdockq2",
         "interface_ipSAE": "ipsae", "interface_LIS": "lis", "interface_cLIS": "clis",
         "interface_iLIS": "ilis", "interface_ccc": "ccc", "interface_score": "interface_score",
         "interface_meta_score": "meta_score"}
DERIVED = [(m, SLUGS[m], NAMES[m]) for m in CATEGORIES["derived"][2]]


def main():
    a, order, _ = load()
    # Figure 1's method order, without the screen (no AlphaJudge data for it)
    methods = [m for m in IPTM_METHODS if m != "AFM_screen"]
    a = a.copy()
    a["method"] = a.method.astype(str)
    for col, slug, name in DERIVED:
        dot_plot(a, col, methods, order, name, HERE / f"fig_derived_{slug}",
                 best_col="is_best")
        print(f"fig_derived_{slug}: {name} | {a[col].notna().sum()} models")


if __name__ == "__main__":
    main()
