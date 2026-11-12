#!/usr/bin/env python3
"""Single source of truth for the control labels.

Changed 2026-11-06. Previously only four complexes were labelled (two positive,
two negative) and the other six were treated as unlabelled screen hits. Every
complex is now labelled:

POSITIVE  interactions demonstrated in the lab.
NEGATIVE  everything else -- "could not be proven in the lab". Note this is
          absence of evidence, not demonstrated non-interaction, and a failed
          pulldown is weaker than a proven negative. Separation statistics
          computed against these labels inherit that weakness.

Import from here; do not re-declare these lists in individual scripts.
"""

POSITIVE = ["P69B_Epi1", "Pip1_EpiC2B", "Pip1_EpiC1", "P69B_SSP1"]

NEGATIVE = ["P69B_Ecp36", "P69B_EpiC1", "P69B_EpiC2B",
            "P69B_Six15", "P69B_TIL", "Pip1_Epi1"]

# P69B_SSP1 moved to POSITIVE on 2026-11-10 (user): SSP1 does interact with P69B,
# but the predicted interface is not the real binding site and the interaction does
# not affect P69B function.

ALL = POSITIVE + NEGATIVE

# Back-compat aliases used by the older scripts.
POS = POSITIVE
NEG = NEGATIVE


def control(complex_name):
    """'positive', 'negative', or 'unlabelled' for complexes outside the set."""
    if complex_name in POSITIVE:
        return "positive"
    if complex_name in NEGATIVE:
        return "negative"
    return "unlabelled"


def marker(complex_name):
    """The (+) / (-) suffix used in figure and table labels."""
    kind = control(complex_name)
    return " (+)" if kind == "positive" else " (−)" if kind == "negative" else ""
