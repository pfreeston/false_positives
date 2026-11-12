#!/usr/bin/env python3
"""Single source of truth for the prediction methods shown in the manuscript.

Added 2026-11-10 when the database-variant runs were added as methods of their own.
Every run uses the mature, AFM-matched sequences (see README). Import from here;
do not re-declare method lists, colours or labels in individual scripts.

key            AlphaJudge summary (5_alpha_judge/)  MSA databases
Boltz2         boltz2_matched_summary.csv           ColabFold MMseqs2 server
Boltz2_reduced boltz2_reduced_summary.csv           AF3 pipeline MSAs (small BFD)
AF3_cluster    af3_matched_summary.csv              AF3 pipeline (small BFD = reduced)
AF3_fulldb     af3_fulldb_summary.csv               AFM full_dbs MSAs (BFD + UniRef30)
AF3_server     af3_server_summary.csv               AlphaFold Server
AFM            afm_matched_summary.csv              AF 2.3.1 full_dbs (BFD + UniRef30)
AFM_reduced    afm_reduced_summary.csv              AF 2.3.1 reduced_dbs (small BFD)

AFM ran 5 models x 5 seeds = 25 predictions per complex. Only the top N_MODELS by
AFM's own ranking (iptm + ptm) are kept, so every method contributes five models.

AFM_screen is the original screen (False_positives.xlsx). It has iptm + ptm only, so
it appears only in Table 1 and the ipTM + pTM dot plot.
"""

N_MODELS = 5

# AlphaJudge methods, in display order (variants next to their parent method)
AJ_METHODS = ["Boltz2", "Boltz2_reduced", "AF3_cluster", "AF3_fulldb",
              "AF3_server", "AFM", "AFM_reduced"]

SUMMARY = {
    "Boltz2": "boltz2_matched_summary.csv",
    "Boltz2_reduced": "boltz2_reduced_summary.csv",
    "AF3_cluster": "af3_matched_summary.csv",
    "AF3_fulldb": "af3_fulldb_summary.csv",
    "AF3_server": "af3_server_summary.csv",
    "AFM": "afm_matched_summary.csv",
    "AFM_reduced": "afm_reduced_summary.csv",
}

# methods whose raw output has more models than N_MODELS -> keep the top-ranked
TOP_N = {"AFM", "AFM_reduced"}

# Figure 1 (ipTM + pTM): every AFM run -- the original screen and the full-DB and
# reduced-DB cluster runs (user, 2026-11-10) -- plus the other methods. Table 1 keeps
# its own column list in make_iptm_table.py (screen only for AFM).
IPTM_METHODS = ["AFM_screen", "AFM", "AFM_reduced", "Boltz2", "Boltz2_reduced",
                "AF3_cluster", "AF3_fulldb", "AF3_server"]

# parent hue = standard run; the darker/lighter partner = its database variant
COLOUR = {
    "AFM_screen": "#E69F00",
    "Boltz2": "#0072B2", "Boltz2_reduced": "#56B4E9",
    "AF3_cluster": "#009E73", "AF3_fulldb": "#004D38",
    "AF3_server": "#CC79A7",
    "AFM": "#D55E00", "AFM_reduced": "#8C510A",
}

PRETTY = {
    "AFM_screen": "AFM screen",
    "Boltz2": "Boltz-2", "Boltz2_reduced": "Boltz-2 reduced-DB",
    "AF3_cluster": "AF3 cluster", "AF3_fulldb": "AF3 cluster full-DB",
    "AF3_server": "AF3 server",
    "AFM": "AFM full-DB", "AFM_reduced": "AFM reduced-DB",
}

SLUG = {m: m.lower() for m in AJ_METHODS}
