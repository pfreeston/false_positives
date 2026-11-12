#!/usr/bin/env python3
"""Build Boltz-2 inputs whose MSAs come from the AF3 data pipeline (reduced DBs).

The AF3 data pipeline searches UniRef90, MGnify, small BFD and UniProt with
Jackhmmer, i.e. the same database set as AlphaFold-Multimer's reduced_dbs.
Boltz-2 normally gets its MSAs from the ColabFold MMseqs2 server instead; this
script replaces those with the AF3 MSAs so Boltz-2 runs on reduced-DB MSAs.

Per complex and chain it writes a Boltz CSV (columns key,sequence), the same
layout boltz writes for server MSAs:
  key 0       the query (paired with the other chain's query)
  key 1..n    paired rows: top UniProt hit per species (OX= taxonomy id) in
              each chain, for species present in both chains, ordered by the
              chain A hit rank (AlphaFold-Multimer's pairing rule)
  key -1      every other row: the unpairedMsa rows, then the remaining
              pairedMsa (UniProt) rows
and a YAML pointing at the two CSVs.

Usage (on the cluster):
  python3 af3_msa_to_boltz.py --pipeline_dir <work_afm_matched/json/pipeline> \
      --yaml_dir <boltz_inputs_afm_matched> --out_dir <boltz_reduced/inputs>
"""
import argparse
import json
import re
from pathlib import Path

COMPLEXES = [
    "P69B_Ecp36", "P69B_Epi1", "P69B_EpiC1", "P69B_EpiC2B", "P69B_SSP1",
    "P69B_Six15", "P69B_TIL", "Pip1_Epi1", "Pip1_EpiC1", "Pip1_EpiC2B",
]
OX = re.compile(r"\bOX=(\d+)")


def read_a3m(text):
    """Return [(header, sequence)] from an a3m string."""
    rows, header = [], None
    for line in (text or "").splitlines():
        if line.startswith(">"):
            header = line[1:]
            rows.append([header, ""])
        elif header is not None:
            rows[-1][1] += line.strip()
    return [tuple(r) for r in rows]


def load_chain(pipeline_dir, name):
    f = pipeline_dir / name.lower() / name.lower() / f"{name.lower()}_data.json"
    prot = json.load(open(f))["sequences"][0]["protein"]
    unpaired = read_a3m(prot.get("unpairedMsa"))
    paired = read_a3m(prot.get("pairedMsa"))
    query = prot["sequence"]
    # Drop the query row both MSAs start with.
    unpaired = [r for i, r in enumerate(unpaired) if not (i == 0 and r[1] == query)]
    paired = [r for i, r in enumerate(paired) if not (i == 0 and r[1] == query)]
    return query, unpaired, paired


def top_hit_per_species(paired):
    """{taxid: row index} for the first (best) UniProt hit of each species."""
    best = {}
    for i, (h, _) in enumerate(paired):
        m = OX.search(h)
        if m and m.group(1) not in best:
            best[m.group(1)] = i
    return best


def write_csv(path, rows):
    with open(path, "w") as fh:
        fh.write("key,sequence\n")
        for key, seq in rows:
            fh.write(f"{key},{seq}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pipeline_dir", required=True, type=Path)
    ap.add_argument("--yaml_dir", required=True, type=Path,
                    help="existing boltz_inputs_afm_matched (sequences are checked against it)")
    ap.add_argument("--out_dir", required=True, type=Path)
    args = ap.parse_args()
    msa_dir = args.out_dir / "msa"
    msa_dir.mkdir(parents=True, exist_ok=True)

    for cx in COMPLEXES:
        a_name, b_name = cx.split("_")
        qa, ua, pa = load_chain(args.pipeline_dir, a_name)
        qb, ub, pb = load_chain(args.pipeline_dir, b_name)

        ref = (args.yaml_dir / f"{cx}.yaml").read_text()
        for q, n in ((qa, a_name), (qb, b_name)):
            if f"sequence: {q}\n" not in ref + "\n":
                raise SystemExit(f"{cx}: AF3 pipeline sequence for {n} differs from {args.yaml_dir}/{cx}.yaml")

        ba, bb = top_hit_per_species(pa), top_hit_per_species(pb)
        shared = sorted(set(ba) & set(bb), key=lambda t: ba[t])
        used_a = {ba[t] for t in shared}
        used_b = {bb[t] for t in shared}

        rows_a = [(0, qa)] + [(k, pa[ba[t]][1]) for k, t in enumerate(shared, 1)]
        rows_b = [(0, qb)] + [(k, pb[bb[t]][1]) for k, t in enumerate(shared, 1)]
        rows_a += [(-1, s) for _, s in ua] + [(-1, s) for i, (_, s) in enumerate(pa) if i not in used_a]
        rows_b += [(-1, s) for _, s in ub] + [(-1, s) for i, (_, s) in enumerate(pb) if i not in used_b]

        ca, cb = msa_dir / f"{cx}_A.csv", msa_dir / f"{cx}_B.csv"
        write_csv(ca, rows_a)
        write_csv(cb, rows_b)
        (args.out_dir / f"{cx}.yaml").write_text(
            "version: 1\nsequences:\n"
            f"  - protein:\n      id: A\n      sequence: {qa}\n      msa: {ca.resolve()}\n"
            f"  - protein:\n      id: B\n      sequence: {qb}\n      msa: {cb.resolve()}\n"
        )
        print(f"{cx}: paired {len(shared)} species; A {len(rows_a)} rows, B {len(rows_b)} rows")


if __name__ == "__main__":
    main()
