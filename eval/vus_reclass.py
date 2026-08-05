"""VUS-reclassification accuracy - the falsifiable headline metric.

% of input VUS reclassified by the joint disease model, **concordant with the VCEP/ClinVar
3-star truth**. This turns "mitigates VUS" from a slogan into a number (the v2 thesis,
echoing Ross 2021's 29%->20% on Glanzmann).
"""
from __future__ import annotations

from dataclasses import dataclass

from core.dx_schemas import VariantState
from diseases.ontology import route_clusters
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import reclassify
from jointdx.orchestrate import _best_cluster

_PATHOGENIC = {VariantState.PATH, VariantState.LP}
_BENIGN = {VariantState.LB, VariantState.BEN}


def _same_direction(a: VariantState, b: VariantState) -> bool:
    return (a in _PATHOGENIC and b in _PATHOGENIC) or (a in _BENIGN and b in _BENIGN) or a == b


@dataclass
class ReclassCase:
    ev: Evidence
    truth: VariantState                 # the VCEP/ClinVar 3-star classification


@dataclass
class ReclassReport:
    n: int
    n_reclassified: float
    reclassified_rate: float
    concordance: float                  # of reclassified, fraction matching truth direction


def evaluate(cases: list[ReclassCase]) -> ReclassReport:
    n = len(cases)
    reclassified = 0
    concordant = 0
    for c in cases:
        clusters = route_clusters([c.ev.variant_gene]) if c.ev.variant_gene else []
        if not clusters:
            continue
        j = joint(_best_cluster(clusters, c.ev), c.ev)
        old, new, _ = reclassify(j, VariantState.VUS)
        if new != old:                                  # moved out of VUS
            reclassified += 1
            if _same_direction(new, c.truth):
                concordant += 1
    return ReclassReport(
        n=n, n_reclassified=reclassified,
        reclassified_rate=reclassified / n if n else 0.0,
        concordance=concordant / reclassified if reclassified else 0.0,
    )
