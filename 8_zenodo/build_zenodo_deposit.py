#!/usr/bin/env python3
"""Stage the predicted structures and confidence data for a Zenodo deposit.

The structures are deliberately NOT in git (user decision 2026-11-10); this script
collects them from the run directories into deposit/ (gitignored) as one zip per
method, plus the input sequences, the AlphaJudge summaries, the README and a
checksum manifest. Nothing is uploaded -- upload deposit/ from the Zenodo web UI.

Per method and complex (all ten complexes, mature AFM-matched sequences):
  Boltz-2, Boltz-2 reduced-DB    model CIFs, confidence JSONs, PAE + pLDDT npz
  AF3 cluster, AF3 cluster full-DB, AF3 server
                                 every seed/sample CIF + confidences + summary
                                 confidences, ranking scores, AF3 terms of use
  AFM full-DB, AFM reduced-DB    all 25 relaxed PDBs, ranking_debug.json, and per
                                 model confidence_<model>.npz (plddt, pae, ptm,
                                 iptm, ranking_confidence) extracted from the
                                 result pickles, which are too large to deposit

Usage: python3 8_zenodo/build_zenodo_deposit.py [--workers 8]
"""
import argparse
import gzip
import hashlib
import pickle
import shutil
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / "deposit"
STAGE = HERE / "staging"

COMPLEXES = ["P69B_Ecp36", "P69B_Epi1", "P69B_EpiC1", "P69B_EpiC2B", "P69B_SSP1",
             "P69B_Six15", "P69B_TIL", "Pip1_Epi1", "Pip1_EpiC1", "Pip1_EpiC2B"]

# zip stem -> (kind, run root, AlphaJudge summary)
METHODS = {
    "boltz2": ("boltz", "2_boltz/boltz_results_afm_matched", "boltz2_matched_summary.csv"),
    "boltz2_reduced_db": ("boltz", "2_boltz/boltz_results_reduced", "boltz2_reduced_summary.csv"),
    "af3_cluster": ("af3_nested", "3_af3/cluster/outputs_afm_matched", "af3_matched_summary.csv"),
    "af3_cluster_full_db": ("af3", "3_af3/cluster/outputs_fulldb", "af3_fulldb_summary.csv"),
    "af3_server": ("af3", "3_af3/server", "af3_server_summary.csv"),
    "afm_full_db": ("afm", "4_afm/afm_matched", "afm_matched_summary.csv"),
    "afm_reduced_db": ("afm", "4_afm/afm_matched_reduced", "afm_reduced_summary.csv"),
}
EXPECTED = {"boltz": 5, "af3": 5, "af3_nested": 5, "afm": 25}   # models per complex


def copy(src: Path, dst_dir: Path):
    dst_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst_dir / src.name)


def stage_boltz(root: Path, cx: str, dst: Path) -> int:
    pred = root / f"boltz_results_{cx}" / "predictions" / cx
    n = 0
    for pat in ("*.cif", "confidence_*.json", "pae_*.npz", "plddt_*.npz"):
        for f in sorted(pred.glob(pat)):
            copy(f, dst)
            n += f.suffix == ".cif"
    return n


def stage_af3(root: Path, cx: str, dst: Path, nested: bool) -> int:
    base = root / cx / cx if nested else root / cx
    n = 0
    for sample in sorted(base.glob("seed-*_sample-*")):
        for f in sorted(sample.glob("*")):
            if f.name.endswith(("_model.cif", "_confidences.json")):
                copy(f, dst / sample.name)
                n += f.name.endswith("_model.cif")
    for pat in ("*_ranking_scores.csv", "TERMS_OF_USE.md", "terms_of_use.md", "*_job_request.json"):
        for f in base.glob(pat):
            copy(f, dst)
    return n


def afm_confidence(args):
    pkl, out = args
    with gzip.open(pkl, "rb") as fh:
        d = pickle.load(fh)
    np.savez_compressed(out, plddt=np.asarray(d["plddt"], np.float32),
                        pae=np.asarray(d["predicted_aligned_error"], np.float32),
                        ptm=np.float32(d["ptm"]), iptm=np.float32(d["iptm"]),
                        ranking_confidence=np.float32(d["ranking_confidence"]))
    return out


def stage_afm(root: Path, cx: str, dst: Path, jobs: list) -> int:
    src = root / cx
    dst.mkdir(parents=True, exist_ok=True)
    for name in ("ranking_debug.json", "relax_metrics.json", "timings.json"):
        if (src / name).exists():
            copy(src / name, dst)
    pdbs = sorted(src.glob("relaxed_*.pdb"))
    for f in pdbs:
        copy(f, dst)
        model = f.stem[len("relaxed_"):]
        pkl = src / f"result_{model}.pkl.gz"
        if not pkl.exists():
            raise SystemExit(f"missing {pkl}")
        jobs.append((pkl, dst / f"confidence_{model}.npz"))
    return len(pdbs)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def zip_dir(src: Path, out: Path):
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for f in sorted(src.rglob("*")):
            if f.is_file():
                z.write(f, f.relative_to(src.parent))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    for d in (OUT, STAGE):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)

    afm_jobs = []
    for stem, (kind, root, _) in METHODS.items():
        root = ROOT / root
        for cx in COMPLEXES:
            dst = STAGE / stem / cx
            if kind == "boltz":
                n = stage_boltz(root, cx, dst)
            elif kind == "afm":
                n = stage_afm(root, cx, dst, afm_jobs)
            else:
                n = stage_af3(root, cx, dst, nested=kind == "af3_nested")
            if n != EXPECTED[kind]:
                raise SystemExit(f"{stem} {cx}: {n} models, expected {EXPECTED[kind]}")
        print(f"staged {stem}")

    print(f"extracting {len(afm_jobs)} AFM confidence files ...")
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for i, _ in enumerate(ex.map(afm_confidence, afm_jobs, chunksize=4), 1):
            if i % 50 == 0:
                print(f"  {i}/{len(afm_jobs)}")

    # inputs and AlphaJudge summaries
    for f in sorted((ROOT / "4_afm_cluster/fastas").glob("*.fasta")):
        copy(f, STAGE / "sequences")
    (STAGE / "alphajudge").mkdir()
    for stem, (_, _, summary) in METHODS.items():
        shutil.copy2(ROOT / "5_alpha_judge" / summary, STAGE / "alphajudge" / f"{stem}_alphajudge.csv")

    manifest = ["zip\tpath\tbytes\tsha256"]
    for sub in sorted(p for p in STAGE.iterdir() if p.is_dir()):
        out = OUT / f"{sub.name}.zip"
        zip_dir(sub, out)
        for f in sorted(sub.rglob("*")):
            if f.is_file():
                manifest.append(f"{out.name}\t{f.relative_to(STAGE)}\t{f.stat().st_size}\t{sha256(f)}")
        print(f"{out.name}: {out.stat().st_size / 1e6:.0f} MB")
    (OUT / "MANIFEST.tsv").write_text("\n".join(manifest) + "\n")
    shutil.copy2(HERE / "README.md", OUT / "README.md")
    shutil.rmtree(STAGE)

    total = sum(f.stat().st_size for f in OUT.iterdir())
    print(f"\n{len(list(OUT.iterdir()))} files in {OUT.relative_to(ROOT)}, {total / 1e9:.2f} GB")
    print("zip sha256:")
    for f in sorted(OUT.glob("*.zip")):
        print(f"  {sha256(f)}  {f.name}")


if __name__ == "__main__":
    main()
