"""Whole-odyssey (cross-modality) action space (plan Phase 4, Step 4.1).

Beyond single-variant evidence acquisition, the case-level action space includes
modality escalations and meta-actions. Each carries cost/turnaround and a
conditional-yield prior (probability of surfacing a new causal finding), seeded from
the multimodal-escalation literature and updated by the learning loop.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ModalityAction:
    name: str
    cost_usd: float
    turnaround_days: int
    yield_prior: float          # P(surfaces a new causal finding) given the failure mode
    failure_mode: str           # which failure this addresses
    note: str = ""


# Literature-seeded priors (explicitly versioned; updated per-institution by Phase 5.2).
MODALITY_CATALOG: list[ModalityAction] = [
    ModalityAction("rna_seq", 900, 28, 0.16, "expression_splicing",
                   "splicing/expression aberrations missed by DNA"),
    ModalityAction("long_read_genome", 1800, 35, 0.12, "detection",
                   "SVs, repeats, hard-to-map regions; phasing"),
    ModalityAction("optical_genome_mapping", 1200, 30, 0.08, "detection",
                   "large structural variants / balanced rearrangements"),
    ModalityAction("methylation_episignature", 700, 21, 0.07, "imprinting_episignature",
                   "imprinting disorders / episignatures"),
    ModalityAction("proteomics", 1500, 42, 0.05, "functional", "protein-level consequence"),
    ModalityAction("reanalysis", 50, 180, 0.10, "knowledge_growth",
                   "low-cost 'wait' action; knowledge-base growth over time"),
    ModalityAction("matchmaking", 0, 120, 0.06, "novel_gene",
                   "GeneMatcher/Matchmaker Exchange for novel-gene candidates"),
]


def actions_for_failure(failure_mode: str) -> list[ModalityAction]:
    return [a for a in MODALITY_CATALOG if a.failure_mode == failure_mode]
