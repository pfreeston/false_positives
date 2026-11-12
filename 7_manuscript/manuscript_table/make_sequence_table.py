"""Build the manuscript sequence-information tables.

Sequences are read from 1_sequence_processing/afm_matched/ (the constructs used
for the AFM-matched runs). Accessions, species and residue ranges were checked
against NCBI Protein / UniProtKB on 2026-11-05: every modelled sequence is an
exact substring of the listed database entry at the stated residue range.

Writes:
  proteins.csv   one row per protein (name, role, species, accession, range, sequence)
  complexes.csv  one row per modelled complex
  tables.md      tables formatted as Markdown
  tables.tex     the same tables in LaTeX (booktabs, tabularx, longtable, seqsplit)
  tables_standalone.tex  minimal document that \\input{tables.tex}, for a compile check
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
FASTA_DIR = HERE.parent / "afm_matched"

# name: (role, description, species, primary accession, other IDs, db residue range)
PROTEINS = {
    "P69B": ("Host hydrolase", "Subtilisin-like serine protease",
             "Solanum lycopersicum", "UniProtKB O04678",
             "Solyc08g079870; EMBL Y10149", (111, 745)),
    "Pip1": ("Host hydrolase", "Papain-like cysteine protease",
             "Solanum lycopersicum", "UniProtKB Q156I2",
             "EMBL DQ631808", (127, 345)),
    "Epi1": ("Effector", "Kazal-like protease inhibitor",
             "Phytophthora infestans T30-4", "UniProtKB D0MVC9",
             "PITG_22681", (17, 149)),
    "EpiC1": ("Effector", "Cystatin-like protease inhibitor",
              "Phytophthora infestans T30-4", "GenBank EEY55256.1",
              "PITG_09169", (22, 126)),
    "EpiC2B": ("Effector", "Cystatin-like protease inhibitor",
               "Phytophthora infestans T30-4", "GenBank EEY55258.1",
               "PITG_09173", (22, 125)),
    "Ecp36": ("Effector", "Extracellular protein Ecp36-1",
              "Fulvia fulva (syn. Cladosporium fulvum)", "GenBank KAH3648627.1",
              "CLAFUR5_14667", (20, 73)),
    "Six15": ("Effector", "Secreted in xylem 15 (SIX15)",
              "Fusarium oxysporum f. sp. lycopersici (race 3)", "GenBank APP91304.1",
              "", (19, 79)),
    "TIL": ("Effector", "Trypsin inhibitor-like cysteine-rich protein",
            "Fusarium oxysporum f. sp. lycopersici 4287", "RefSeq XP_018243121.1",
            "FOXG_19485", (23, 100)),
    "SSP1": ("Effector", "Hypothetical protein",
             "Xanthomonas spp. (predominantly X. euvesicatoria)", "RefSeq WP_046932418.1",
             "", (14, 126)),
}

HYDROLASES = ["P69B", "Pip1"]
EFFECTORS = ["Epi1", "EpiC1", "EpiC2B", "Ecp36", "Six15", "TIL", "SSP1"]


ITALIC = {
    "Solanum lycopersicum": "*Solanum lycopersicum*",
    "Phytophthora infestans T30-4": "*Phytophthora infestans* T30-4",
    "Fulvia fulva (syn. Cladosporium fulvum)": "*Fulvia fulva* (syn. *Cladosporium fulvum*)",
    "Fusarium oxysporum f. sp. lycopersici (race 3)": "*Fusarium oxysporum* f. sp. *lycopersici* (race 3)",
    "Fusarium oxysporum f. sp. lycopersici 4287": "*Fusarium oxysporum* f. sp. *lycopersici* 4287",
    "Xanthomonas spp. (predominantly X. euvesicatoria)": "*Xanthomonas* spp. (predominantly *X. euvesicatoria*)",
}


# Abbreviated names for the complexes table (full names are in the proteins table).
SHORT = {
    "Phytophthora infestans T30-4": "*P. infestans*",
    "Fulvia fulva (syn. Cladosporium fulvum)": "*F. fulva*",
    "Fusarium oxysporum f. sp. lycopersici (race 3)": "*F. oxysporum* f. sp. *lycopersici*",
    "Fusarium oxysporum f. sp. lycopersici 4287": "*F. oxysporum* f. sp. *lycopersici*",
    "Xanthomonas spp. (predominantly X. euvesicatoria)": "*Xanthomonas* spp.",
}


def read_fasta(path):
    return "".join(l.strip() for l in path.read_text().splitlines() if not l.startswith(">"))


def main():
    rows = []
    for name, (role, desc, sp, acc, other, (start, end)) in PROTEINS.items():
        seq = read_fasta(FASTA_DIR / f"{name}.fasta")
        assert len(seq) == end - start + 1, (name, len(seq), start, end)
        rows.append(dict(Protein=name, Role=role, Description=desc, Species=sp,
                         Accession=acc, Other_IDs=other,
                         Modelled_residues=f"{start}–{end}", Length_aa=len(seq),
                         Sequence=seq))

    with open(HERE / "proteins.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    info = {r["Protein"]: r for r in rows}
    complexes = []
    for h in HYDROLASES:
        for e in EFFECTORS:
            complexes.append(dict(
                Complex=f"{h}_{e}",
                Chain_A=h, Chain_A_accession=info[h]["Accession"],
                Chain_A_species=info[h]["Species"], Chain_A_length_aa=info[h]["Length_aa"],
                Chain_B=e, Chain_B_accession=info[e]["Accession"],
                Chain_B_species=info[e]["Species"], Chain_B_length_aa=info[e]["Length_aa"],
            ))
    with open(HERE / "complexes.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=complexes[0].keys())
        w.writeheader()
        w.writerows(complexes)

    md = ["# Sequence information", "",
          "**Table S1. Proteins modelled in this study.** Sequences are the mature "
          "constructs (signal peptide and, for the proteases, prodomain removed) used "
          "for structure prediction. Residue numbering refers to the listed database entry.",
          "",
          "| Protein | Role | Description | Species | Accession | Other IDs | Residues | Length (aa) |",
          "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['Protein']} | {r['Role']} | {r['Description']} | {ITALIC[r['Species']]} | "
                  f"{r['Accession']} | {r['Other_IDs']} | {r['Modelled_residues']} | {r['Length_aa']} |")
    md += ["", "**Table S2. Modelled complexes.** Chain A, host hydrolase; chain B, effector.", "",
           "| Complex | Chain A (accession) | Chain B (accession) | Chain B species | Total length (aa) |",
           "|---|---|---|---|---|"]
    for c in complexes:
        md.append(f"| {c['Complex']} | {c['Chain_A']} ({c['Chain_A_accession']}) | "
                  f"{c['Chain_B']} ({c['Chain_B_accession']}) | {ITALIC[c['Chain_B_species']]} | "
                  f"{c['Chain_A_length_aa'] + c['Chain_B_length_aa']} |")
    md += ["", "**Table S3. Modelled sequences.**", ""]
    for r in rows:
        md += [f"**{r['Protein']}** ({r['Accession']}, residues {r['Modelled_residues']}, "
               f"{r['Length_aa']} aa)", "", "```", r["Sequence"], "```", ""]
    (HERE / "tables.md").write_text("\n".join(md))

    write_latex(rows, complexes)


def tex(s):
    """Markdown-italic species string / plain text -> LaTeX."""
    s = s.replace("_", r"\_").replace("\u2013", "--")
    return re.sub(r"\*([^*]+)\*", r"\\textit{\1}", s)


def write_latex(rows, complexes):
    L = [r"% Generated by make_sequence_table.py -- do not edit by hand.",
         r"% Requires in the preamble: \usepackage{booktabs,tabularx,longtable,seqsplit}",
         "",
         r"\begin{table}[htbp]",
         r"\centering\small",
         r"\caption{Proteins modelled in this study. Sequences are the mature constructs "
         r"(signal peptide and, for the proteases, prodomain removed) used for structure "
         r"prediction. Residue numbering refers to the listed database entry.}",
         r"\label{tab:proteins}",
         r"\begin{tabularx}{\textwidth}{@{}l>{\raggedright\arraybackslash}X"
         r">{\raggedright\arraybackslash}X>{\raggedright\arraybackslash}p{3.3cm}rr@{}}",
         r"\toprule",
         r"Protein & Description & Species & Accession & Residues & Length (aa) \\",
         r"\midrule"]
    for group, label in (("Host hydrolase", "Host hydrolases"), ("Effector", "Effectors")):
        L.append(rf"\multicolumn{{6}}{{@{{}}l}}{{\textit{{{label}}}}} \\")
        for r in rows:
            if r["Role"] != group:
                continue
            acc = tex(r["Accession"])
            if r["Other_IDs"]:
                acc += r"\newline{\footnotesize " + tex(r["Other_IDs"]) + "}"
            L.append(" & ".join([r["Protein"], tex(r["Description"]), tex(ITALIC[r["Species"]]),
                                 acc, tex(r["Modelled_residues"]), str(r["Length_aa"])]) + r" \\")
        if group == "Host hydrolase":
            L.append(r"\addlinespace")
    L += [r"\bottomrule", r"\end{tabularx}", r"\end{table}", "",
          r"\begin{table}[htbp]",
          r"\centering\small",
          r"\caption{Modelled complexes. Chain A, host hydrolase (\textit{Solanum lycopersicum}); "
          r"chain B, effector. Length is the combined length of both chains. Full species "
          r"names and accessions are given in Table~\ref{tab:proteins}.}",
          r"\label{tab:complexes}",
          r"\begin{tabularx}{\textwidth}{@{}lll>{\raggedright\arraybackslash}Xr@{}}",
          r"\toprule",
          r"Complex & Chain A & Chain B (accession) & Chain B source organism & Length (aa) \\",
          r"\midrule"]
    for i, c in enumerate(complexes):
        if i and c["Chain_A"] != complexes[i - 1]["Chain_A"]:
            L.append(r"\addlinespace")
        L.append(" & ".join([tex(c["Complex"]), c["Chain_A"],
                             f"{c['Chain_B']} ({tex(c['Chain_B_accession'].split()[-1])})",
                             tex(SHORT[c["Chain_B_species"]]),
                             str(c["Chain_A_length_aa"] + c["Chain_B_length_aa"])]) + r" \\")
    L += [r"\bottomrule", r"\end{tabularx}", r"\end{table}", "",
          r"{\small",
          r"\begin{longtable}{@{}l>{\ttfamily\footnotesize}p{0.78\textwidth}@{}}",
          r"\caption{Amino-acid sequences of the modelled constructs.}\label{tab:sequences}\\",
          r"\toprule Protein & \normalfont Sequence \\ \midrule \endfirsthead",
          r"\toprule Protein & \normalfont Sequence \\ \midrule \endhead",
          r"\bottomrule \endlastfoot"]
    for r in rows:
        L.append(rf"\begin{{tabular}}[t]{{@{{}}l@{{}}}}{r['Protein']}\\"
                 rf"\scriptsize {tex(r['Modelled_residues'])} ({r['Length_aa']} aa)\end{{tabular}}"
                 rf" & \seqsplit{{{r['Sequence']}}} \\ \addlinespace")
    L += [r"\end{longtable}", r"}", ""]
    (HERE / "tables.tex").write_text("\n".join(L))

    (HERE / "tables_standalone.tex").write_text("\n".join([
        r"\documentclass[11pt]{article}",
        r"\usepackage[a4paper,margin=2cm]{geometry}",
        r"\usepackage{booktabs,tabularx,longtable,seqsplit}",
        r"\begin{document}",
        r"\input{tables.tex}",
        r"\end{document}", ""]))


if __name__ == "__main__":
    main()
