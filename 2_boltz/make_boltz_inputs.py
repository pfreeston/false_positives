import argparse
import os

complexes = {
    "Pip1_EpiC1":  ("Pip1",  "EpiC1"),
    "Pip1_EpiC2B": ("Pip1",  "EpiC2B"),
    "Pip1_Epi1":   ("Pip1",  "Epi1"),
    "P69B_Epi1":   ("P69B",  "Epi1"),
    "P69B_EpiC1":  ("P69B",  "EpiC1"),
    "P69B_EpiC2B": ("P69B",  "EpiC2B"),
    "P69B_Ecp36":  ("P69B",  "Ecp36"),
    "P69B_Six15":  ("P69B",  "Six15"),
    "P69B_TIL":    ("P69B",  "TIL"),
    "P69B_SSP1":   ("P69B",  "SSP1"),
}

def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)

parser = argparse.ArgumentParser(description="Write Boltz-2 input YAMLs for protein pairs.")
parser.add_argument("seq_dir", help="Directory containing <name>.fasta files")
parser.add_argument("-o", "--out_dir", default="boltz_inputs",
                    help="Output directory for YAML files (default: boltz_inputs)")
args = parser.parse_args()

os.makedirs(args.out_dir, exist_ok=True)

for complex_name, (p1, p2) in complexes.items():
    seq1 = read_fasta(os.path.join(args.seq_dir, f"{p1}.fasta"))
    seq2 = read_fasta(os.path.join(args.seq_dir, f"{p2}.fasta"))
    yaml_content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {seq1}
  - protein:
      id: B
      sequence: {seq2}
"""
    out_path = os.path.join(args.out_dir, f"{complex_name}.yaml")
    with open(out_path, "w") as f:
        f.write(yaml_content)
    print(f"Written: {out_path}  ({len(seq1)} + {len(seq2)} aa)")