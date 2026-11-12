#!/usr/bin/env python3
"""Cross-method consensus ranking of interface confidence.

For each complex and method, takes the top-ranked model (highest
confidence_score, i.e. the model each predictor ranks first) from the
per-source AlphaJudge summaries and reports interface metrics. Complexes are
ranked by the mean interface_meta_score across the methods available.

Tiers (on top-ranked models):
  High      every method has ipSAE >= 0.5 and meta >= 0.6
  Mixed     at least one method has ipSAE >= 0.3 (method/construct-dependent)
  Low       no method reaches ipSAE 0.3

Only the AFM-matched (mature-construct) runs are compared. The prodomain runs
(P69B 723 aa / Pip1 316 aa) and their summaries were deleted on 2026-11-10, so the
old default mode is gone; --matched is still accepted and does nothing extra.

Sources: AFM (afm_matched_summary.csv,
the 2026-11-09 cluster rerun: multimer v3, full_dbs, 25 predictions), AF3 server
(af3_server_summary.csv; P69B 634 aa, and the EpiC2B runs keep the 21 aa
signal peptide), AF3 cluster
(af3_matched_summary.csv) and Boltz2 (boltz2_matched_summary.csv), all using
the AFM sequences (P69B 635 aa mature, Pip1 mature; Ecp36 without the FLAG
tag the old AFM screen carried). Writes consensus_summary_matched.csv.
"""
import argparse
import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
MATCHED_SOURCES = {  # file stem -> (label, construct)
    "afm_matched": ("AFM", "mature"),
    "af3_server": ("AF3_server", "mature"),
    "af3_matched": ("AF3_cluster", "mature"),
    "boltz2_matched": ("Boltz2", "mature"),
}
SOURCES = MATCHED_SOURCES
METRICS = ["interface_meta_score", "interface_ipSAE", "iptm", "average_interface_pae", "interface_pDockQ2"]


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def mean(xs):
    xs = [x for x in xs if not math.isnan(x)]
    return sum(xs) / len(xs) if xs else math.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matched", action="store_true", help="accepted for back-compat; always on")
    ap.parse_args()
    out_name = "consensus_summary_matched.csv"
    data = {}
    for stem, (label, _) in SOURCES.items():
        for r in csv.DictReader((HERE / f"{stem}_summary.csv").open()):
            data.setdefault(r["jobs"], {}).setdefault(label, []).append(r)

    out = []
    for job, by_src in data.items():
        row = {"complex": job}
        metas, ipsaes = [], []
        for stem, (label, construct) in SOURCES.items():
            rs = by_src.get(label)
            if not rs:
                continue
            top = max(rs, key=lambda r: f(r["confidence_score"]))
            row[f"{label}_model_used"] = top["model_used"]
            row[f"{label}_structure_file"] = top["structure_file"]
            for m in METRICS:
                row[f"{label}_{m}"] = round(f(top[m]), 3)
            row[f"{label}_ipSAE_models_ge_0.5"] = f"{sum(f(r['interface_ipSAE']) >= 0.5 for r in rs)}/{len(rs)}"
            metas.append(f(top["interface_meta_score"]))
            ipsaes.append(f(top["interface_ipSAE"]))
        row["n_methods"] = len(metas)
        row["mean_meta"] = round(mean(metas), 3)
        row["min_meta"] = round(min(metas), 3)
        row["mean_ipSAE"] = round(mean(ipsaes), 3)
        row["min_ipSAE"] = round(min(ipsaes), 3)
        row["max_ipSAE"] = round(max(ipsaes), 3)
        if all(i >= 0.5 for i in ipsaes) and min(metas) >= 0.6:
            row["tier"] = "High"
        elif max(ipsaes) >= 0.3:
            row["tier"] = "Mixed"
        else:
            row["tier"] = "Low"
        out.append(row)

    tier_order = {"High": 0, "Mixed": 1, "Low": 2}
    out.sort(key=lambda r: (tier_order[r["tier"]], -r["mean_meta"]))
    for i, r in enumerate(out, 1):
        r["rank"] = i

    lead = ["rank", "complex", "tier", "n_methods", "mean_meta", "min_meta", "mean_ipSAE", "min_ipSAE", "max_ipSAE"]
    rest = [f"{label}_{m}" for label, _ in SOURCES.values()
            for m in ["model_used", "structure_file"] + METRICS + ["ipSAE_models_ge_0.5"]]
    with (HERE / out_name).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=lead + rest)
        w.writeheader()
        w.writerows(out)

    labels = [l for l, _ in SOURCES.values()]
    print(f"{'#':>2} {'complex':12s} {'tier':6s} {'meanMeta':>8s} " + " ".join(f"{l:>15s}" for l in labels))
    print(f"{'':>2} {'':12s} {'':6s} {'':>8s} " + " ".join(f"{'meta/ipSAE':>15s}" for _ in labels))
    for r in out:
        cells = []
        for l in labels:
            if f"{l}_interface_meta_score" in r:
                cells.append(f"{r[f'{l}_interface_meta_score']:.2f}/{r[f'{l}_interface_ipSAE']:.2f}")
            else:
                cells.append("-")
        print(f"{r['rank']:>2} {r['complex']:12s} {r['tier']:6s} {r['mean_meta']:>8.2f} " + " ".join(f"{c:>15s}" for c in cells))


if __name__ == "__main__":
    main()
