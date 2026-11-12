#!/usr/bin/env python3
"""Assemble the manuscript CSVs from the AlphaJudge summaries (mature constructs).

"Mature" = no prodomain and no signal peptide, i.e. the AFM-matched construct set:
  P69B  mature catalytic domain (634/635 aa, LSLHTT.../SLHTT...)  -- NOT the 723 aa prodomain form
  Pip1  mature (219 aa, VTEVPNSMDW...)                            -- NOT the 316 aa prodomain form
  EpiC2B without its 21 aa signal peptide (104 aa)
  Ecp36 without AFM's N-terminal FLAG+linker (54 aa)

Sources: one AlphaJudge summary per method, listed in methods.py (SUMMARY), plus
  5_alpha_judge/consensus_summary_matched.csv  per-complex consensus + tiers
AFM methods ran 25 predictions per complex; only the top methods.N_MODELS by AFM's own
ranking (iptm + ptm) are kept, so every method contributes five models.

Outputs (in this folder):
  alphajudge_mature_by_model.csv    every scored model, all methods, method column added
  alphajudge_mature_by_complex.csv  top model per complex per method, with tiers
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AJ = ROOT / "5_alpha_judge"

import sys
sys.path.insert(0, str(HERE))
from methods import AJ_METHODS, N_MODELS, SUMMARY, TOP_N  # noqa: E402

MATURE_SUMMARIES = {m: AJ / SUMMARY[m] for m in AJ_METHODS}


def alphajudge_by_model():
    frames = []
    for method, path in MATURE_SUMMARIES.items():
        d = pd.read_csv(path)
        d.insert(0, "method", method)
        d.insert(1, "complex", d.pop("jobs"))
        d["source_csv"] = path.relative_to(ROOT).as_posix()
        if method in TOP_N:
            d = (d.sort_values("confidence_score", ascending=False)
                  .groupby("complex", sort=False).head(N_MODELS))
        n = d.groupby("complex").size()
        if d.complex.nunique() != 10 or (n != N_MODELS).any():
            raise SystemExit(f"{method}: expected 10 complexes x {N_MODELS} models, got\n{n}")
        frames.append(d)
    out = pd.concat(frames, ignore_index=True)
    return out.sort_values(["complex", "method", "interface_ipSAE"],
                           ascending=[True, True, False]).reset_index(drop=True)


def main():
    aj = alphajudge_by_model()
    aj.to_csv(HERE / "alphajudge_mature_by_model.csv", index=False)

    cons = pd.read_csv(AJ / "consensus_summary_matched.csv")
    cons.to_csv(HERE / "alphajudge_mature_by_complex.csv", index=False)

    outs = [("alphajudge_mature_by_model", aj), ("alphajudge_mature_by_complex", cons)]

    for name, d in outs:
        ncx = d.complex.nunique() if "complex" in d.columns else "-"
        print(f"{name:34s} {len(d):4d} rows  {len(d.columns):3d} cols  {ncx} complexes")
    print("\nmodels per method (AlphaJudge):")
    print(aj.groupby("method").agg(models=("complex", "size"), complexes=("complex", "nunique")).to_string())


if __name__ == "__main__":
    main()
