"""Equitable-evidence routing.

In action selection, the engine penalizes actions whose only yield is a
low-reliability (ancestry-biased) code, and prefers ancestry-robust actions
(functional assays, segregation) for under-represented-ancestry patients - turning
the documented equity gap into a routing rule. Concretely relevant: Glanzmann / BSS /
LAD cohorts are enriched in consanguineous populations where allele-frequency
evidence is least reliable.
"""
from __future__ import annotations

from core.schemas import Action, PatientContext

# Modalities whose evidence is ancestry-equitable (not skewed by reference databases).
EQUITABLE_MODALITIES = {"functional", "rna", "segregation"}
# Codes that are ancestry-biased (frequency + in-silico).
_BIASED_CODES = {"PM2", "BS1", "BA1", "PP3", "BP4"}


def is_equitable(action: Action) -> bool:
    """True if the action yields ancestry-robust evidence (not only biased codes)."""
    if action.modality in EQUITABLE_MODALITIES:
        return True
    return any(c not in _BIASED_CODES for c in action.yields_codes)


def equity_weight(action: Action, patient: PatientContext,
                  low_confidence_threshold: float = 0.7) -> float:
    """A multiplier in (0, 1] applied to an action's value when the patient is from an
    under-represented group (low ancestry confidence) and the action is non-equitable.

    1.0 = no penalty; <1.0 down-weights actions that lean on ancestry-biased evidence.
    """
    under_represented = (patient.ancestry_group is not None
                         and patient.ancestry_confidence < low_confidence_threshold)
    if not under_represented or is_equitable(action):
        return 1.0
    return 0.5            # halve the value of biased-only actions for these patients


def equity_filter(actions: list[Action], patient: PatientContext) -> list[Action]:
    """Drop actions whose ONLY yield is ancestry-biased for under-represented patients
    when ancestry-robust alternatives exist (route to equitable evidence)."""
    under_represented = (patient.ancestry_group is not None
                         and patient.ancestry_confidence < 0.7)
    if not under_represented:
        return actions
    robust = [a for a in actions if is_equitable(a)]
    return robust or actions      # never strand the case if only biased actions exist
