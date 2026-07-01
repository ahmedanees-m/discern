"""Evidence-reliability flags (plan Phase 3, Step 3.1).

Documented fact: non-European-ancestry patients receive VUS results 1.5-2x more
often, driven by the allele-frequency (PM2/BS1) and in-silico (PP3/BP4) codes being
trained on European-skewed data; functional (PS3/BS3) evidence is ancestry-equitable.

This module computes a reliability r in [0, 1] for ancestry-biased codes as a
function of group-specific reference data density / predictor training composition.
`r` multiplies the contribution's effective points in the ledger (the
`EvidenceContribution.reliability` field), surfacing the uncertainty explicitly
rather than silently trusting biased evidence.
"""
from __future__ import annotations

import math

from core.schemas import PointsLedger

# Codes whose reliability is ancestry-sensitive (frequency + in-silico evidence).
ANCESTRY_BIASED = {"PM2", "BS1", "BA1", "PP3", "BP4"}


def reliability_pm2(gnomad_group_n: int) -> float:
    """Reliability of an 'absent/rare' frequency call given the group's reference N.

    Sparse reference data for an under-represented group -> a low-confidence PM2/BS1
    (we cannot be sure the variant is truly rare in that group). Saturates to 1.0 as
    the group reference size grows. r = clip(log10(N+1)/5, 0.3, 1.0).
    """
    return float(min(1.0, max(0.3, math.log10(gnomad_group_n + 1) / 5.0)))


def reliability_pp3(predictor_noneur_fraction: float) -> float:
    """Reliability of an in-silico call as a function of the predictor's training
    representation for the patient's group (more non-European training -> higher r)."""
    return float(min(1.0, max(0.3, 0.5 + predictor_noneur_fraction)))


def apply_reliability(ledger: PointsLedger, group_n: int | None = None,
                      predictor_noneur_fraction: float = 0.5) -> PointsLedger:
    """Down-weight ancestry-biased contributions in place; return the ledger.

    `group_n` is the gnomAD reference size for the patient's genetic-similarity group;
    None means unknown (treated as sparse -> conservative down-weighting).
    """
    r_freq = reliability_pm2(group_n if group_n is not None else 0)
    r_insilico = reliability_pp3(predictor_noneur_fraction)
    for c in ledger.contributions:
        base = c.code.split("_", 1)[0]
        if base in {"PM2", "BS1", "BA1"}:
            c.reliability = r_freq
        elif base in {"PP3", "BP4"}:
            c.reliability = r_insilico
    return ledger
