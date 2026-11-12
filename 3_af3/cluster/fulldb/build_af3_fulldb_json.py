#!/usr/bin/env python3
"""Build an AF3 input JSON whose MSAs are AlphaFold-Multimer full_dbs MSAs.

AF3's public data pipeline only searches small BFD (Jackhmmer). The AF3 paper
also used Uniclust30+BFD (HHblits) at training and for the PoseBusters
evaluation (Supplement 2.2, Table 1). To run AF3 on full-database MSAs, this
script takes the per-chain MSAs that AF 2.3.1 run_alphafold.py computed with
--db_preset=full_dbs (afm_matched/outputs/<Complex>/msas/{A,B}) and puts them in
the AF3 JSON:

  unpairedMsa = UniRef90 (<=10,000) + BFD/UniRef30 HHblits (no cap) +
                MGnify (<=5,000; AF 2.3.1 itself stored at most 500)
                stacked in that order and deduplicated, as AF3's pipeline does
  pairedMsa   = UniProt (<=50,000), not deduplicated, as AF3's pipeline does
  templates   = unset, so AF3 runs its own template search from unpairedMsa

The caps are AF3's (Supplement Table 1). Model seeds are copied from the AF3
matched run's input JSON so only the MSAs (and templates derived from them)
differ from that run.

Runs inside the AF3 container (needs the alphafold3 package):
  python3 build_af3_fulldb_json.py --name P69B_Epi1 --msas <.../msas> \
      --seed_json <.../p69b__epi1.json> --out inputs/P69B_Epi1.json
"""
import argparse
import json
from pathlib import Path

from alphafold3.data import msa, parsers
from alphafold3.constants import mmcif_names

CAPS = {"uniref90": 10_000, "mgnify": 5_000, "uniprot": 50_000}


def sto_to_a3m(path: Path, cap: int) -> str:
    with open(path) as fh:
        return parsers.convert_stockholm_to_a3m(fh, max_sequences=cap)


def chain_msas(msa_dir: Path, query: str):
    a3ms = [
        sto_to_a3m(msa_dir / "uniref90_hits.sto", CAPS["uniref90"]),
        (msa_dir / "bfd_uniref_hits.a3m").read_text(),
        sto_to_a3m(msa_dir / "mgnify_hits.sto", CAPS["mgnify"]),
    ]
    unpaired = msa.Msa.from_multiple_a3ms(
        a3ms, chain_poly_type=mmcif_names.PROTEIN_CHAIN, deduplicate=True)
    paired = msa.Msa.from_multiple_a3ms(
        [sto_to_a3m(msa_dir / "uniprot_hits.sto", CAPS["uniprot"])],
        chain_poly_type=mmcif_names.PROTEIN_CHAIN, deduplicate=False)
    for name, m in (("unpaired", unpaired), ("paired", paired)):
        if m.query_sequence != query:
            raise SystemExit(f"{msa_dir}: {name} MSA query does not match the input sequence")
    return unpaired, paired


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True)
    ap.add_argument("--msas", required=True, type=Path, help="AFM msas/ dir with A/ and B/")
    ap.add_argument("--seed_json", required=True, type=Path, help="AF3 matched-run input JSON (seeds + sequences)")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    ref = json.load(open(args.seed_json))
    chain_map = json.load(open(args.msas / "chain_id_map.json"))
    seqs = []
    for letter, ref_seq in zip("AB", ref["sequences"]):
        query = ref_seq["protein"]["sequence"]
        if chain_map[letter]["sequence"] != query:
            raise SystemExit(f"{args.name} chain {letter}: AFM sequence differs from {args.seed_json}")
        unpaired, paired = chain_msas(args.msas / letter, query)
        print(f"{args.name} {letter}: unpaired {unpaired.depth} rows, paired {paired.depth} rows")
        seqs.append({"protein": {
            "id": letter,
            "sequence": query,
            "unpairedMsa": unpaired.to_a3m(),
            "pairedMsa": paired.to_a3m(),
        }})

    out = {
        "name": args.name,
        "modelSeeds": ref["modelSeeds"],
        "sequences": seqs,
        "dialect": "alphafold3",
        "version": ref.get("version", 2),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(args.out, "w"))


if __name__ == "__main__":
    main()
