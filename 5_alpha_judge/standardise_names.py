#!/usr/bin/env python3
"""Standardise prediction output names so AlphaJudge job names match the
canonical complex names in False_positives.xlsx (e.g. P69B_SSP1, Pip1_EpiC2B).

AlphaJudge uses the *run directory name* as the job name, so every run
directory is renamed to <Complex>. Layouts after running:

  2_boltz/boltz_results/boltz_results_<Complex>/predictions/<Complex>/   (unchanged, already canonical)
  3_af3/cluster/outputs/<Complex>/<Complex>/<Complex>_ranking_scores.csv ...
  3_af3/server/<Complex>/<Complex>_ranking_scores.csv
                        /seed-<seed>_sample-<n>/<Complex>_seed-<seed>_sample-<n>_{model.cif,confidences.json,summary_confidences.json}
  4_afm/from_cluster/<Complex>/ranking_debug.json.gz ...

AF3 server downloads (fold_<date>_<job>_model_N.cif, full_data_N.json,
summary_confidences_N.json) are not recognised by AlphaJudge, so they are
restructured into the local-AF3 layout. Server model_N is the server's rank N
(single seed), stored as sample-N; a <Complex>_ranking_scores.csv is written
from each model's ranking_score (rows in server order, so ties keep it).

Every move is appended to rename_log.tsv (old<TAB>new) so it can be reversed.
Usage: python standardise_names.py [--apply]   (default is a dry run)
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LOG = Path(__file__).resolve().parent / "rename_log.tsv"

# Canonical names from the spreadsheet "Complex" column (+ Pip1 controls run on AF3 cluster).
HYDROLASES = {"p69b": "P69B", "pip1": "Pip1"}
EFFECTORS = {
    "ecp36": "Ecp36", "epi1": "Epi1", "epic1": "EpiC1", "epic2b": "EpiC2B",
    "six15": "Six15", "ssp1": "SSP1", "til": "TIL",
    "fotil": "TIL", "xpssp1": "SSP1",  # AF3 server job names
}


def canonical(raw: str) -> str:
    parts = [p for p in re.split(r"_+", raw.lower()) if p]
    if len(parts) != 2 or parts[0] not in HYDROLASES or parts[1] not in EFFECTORS:
        sys.exit(f"Cannot map '{raw}' to a canonical complex name")
    return f"{HYDROLASES[parts[0]]}_{EFFECTORS[parts[1]]}"


class Mover:
    def __init__(self, apply: bool):
        self.apply = apply
        self.moves: list[tuple[Path, Path]] = []

    def mv(self, src: Path, dst: Path):
        if src == dst:
            return
        if dst.exists():
            sys.exit(f"Refusing to overwrite existing {dst}")
        self.moves.append((src, dst))
        print(f"  {src.relative_to(ROOT)}  ->  {dst.relative_to(ROOT)}")
        if self.apply:
            dst.parent.mkdir(parents=True, exist_ok=True)
            src.rename(dst)
            with LOG.open("a") as f:
                f.write(f"{src}\t{dst}\n")

    def write(self, path: Path, text: str):
        print(f"  [new] {path.relative_to(ROOT)}")
        if self.apply:
            if path.exists():
                sys.exit(f"Refusing to overwrite existing {path}")
            path.write_text(text)
            with LOG.open("a") as f:
                f.write(f"<created>\t{path}\n")


def af3_cluster(m: Mover, subdir: str = "3_af3/cluster/outputs"):
    base = ROOT / subdir
    for outer in sorted(p for p in base.iterdir() if p.is_dir()):
        old = outer.name
        new = canonical(old)
        if old == new:
            continue
        print(f"[af3 cluster] {old} -> {new}")
        inner = outer / old
        # rename files carrying the old prefix (deepest first), then dirs
        for f in sorted(inner.rglob(f"{old}_*"), key=lambda p: -len(p.parts)):
            m.mv(f, f.with_name(new + f.name[len(old):]))
        m.mv(inner, outer / new)
        m.mv(outer, base / new)


def af3_server(m: Mover):
    base = ROOT / "3_af3/server"
    for run in sorted(p for p in base.iterdir() if p.is_dir()):
        req = next(run.glob("fold_*_job_request.json"), None)
        if req is None:
            continue  # already converted
        prefix = req.name[: -len("_job_request.json")]  # fold_<date>_<job>
        job = json.loads(req.read_text())
        job = job[0] if isinstance(job, list) else job
        seeds = job.get("modelSeeds") or []
        if len(seeds) != 1:
            sys.exit(f"{run}: expected exactly one seed, got {seeds}")
        seed = seeds[0]
        new = canonical(run.name)
        print(f"[af3 server] {run.name} -> {new} (seed {seed})")

        rows = []
        for cif in sorted(run.glob(f"{prefix}_model_*.cif"), key=lambda p: int(p.stem.rsplit("_", 1)[1])):
            n = int(cif.stem.rsplit("_", 1)[1])
            sample = f"seed-{seed}_sample-{n}"
            summ = run / f"{prefix}_summary_confidences_{n}.json"
            full = run / f"{prefix}_full_data_{n}.json"
            rows.append((seed, n, json.loads(summ.read_text())["ranking_score"]))
            sd = run / sample
            m.mv(cif, sd / f"{new}_{sample}_model.cif")
            m.mv(full, sd / f"{new}_{sample}_confidences.json")
            m.mv(summ, sd / f"{new}_{sample}_summary_confidences.json")

        csv_text = "seed,sample,ranking_score\n" + "".join(f"{s},{n},{r}\n" for s, n, r in rows)
        m.write(run / f"{new}_ranking_scores.csv", csv_text)
        m.mv(req, run / f"{new}_job_request.json")
        for f in sorted(list(run.glob("msas/*")) + list(run.glob("templates/*"))):
            if f.name.startswith(prefix):
                m.mv(f, f.with_name(new + f.name[len(prefix):]))
        m.mv(run, base / new)


def afm_cluster(m: Mover):
    base = ROOT / "4_afm/from_cluster"
    if not base.is_dir():   # moved to the top-level from_cluster/ after renaming
        return
    for run in sorted(p for p in base.iterdir() if p.is_dir()):
        new = canonical(run.name)
        if run.name != new:
            print(f"[afm cluster] {run.name} -> {new}")
            m.mv(run, base / new)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="perform the moves (default: dry run)")
    args = ap.parse_args()
    m = Mover(args.apply)
    # 3_af3/cluster/outputs (prodomain runs) was deleted on 2026-11-10
    af3_cluster(m, "3_af3/cluster/outputs_afm_matched")
    af3_server(m)
    afm_cluster(m)
    print(f"\n{len(m.moves)} moves {'applied' if args.apply else 'planned (dry run)'}")


if __name__ == "__main__":
    main()
