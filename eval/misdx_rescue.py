"""Management-aware misdiagnosis-rescue - case-control.

On documented **mislabelled** cases ("2B"->PT-VWD; "ITP"->BSS), does DISCERN flag the
danger from the **pre-correction** evidence (corrected label hidden, no leakage)? Reported
as small / not powered, pre-registered.
"""
from __future__ import annotations

from dataclasses import dataclass

from diseases.ontology import route_clusters
from jointdx.factorgraph import Evidence, joint
from jointdx.orchestrate import _best_cluster
from safety.interlock import flags


@dataclass
class MisdxCase:
    ev: Evidence                 # PRE-correction evidence only (the corrected label is hidden)
    true_disease_id: str         # what it actually was (held out for scoring)
    planned_tx: str | None = None


@dataclass
class RescueReport:
    n: int
    rescued: float               # fraction where a flag named the true disease as a danger
    rescue_rate: float


def evaluate(cases: list[MisdxCase]) -> RescueReport:
    n = len(cases)
    rescued = 0
    for c in cases:
        clusters = route_clusters([c.ev.variant_gene]) if c.ev.variant_gene else []
        if not clusters:
            continue
        cl = _best_cluster(clusters, c.ev)
        fl = flags(cl, joint(cl, c.ev), planned_tx=c.planned_tx)
        # rescued if a safety flag named the true (dangerous) disease as a competitor
        if any(f.competitor_id == c.true_disease_id or f.leading_id == c.true_disease_id
               for f in fl):
            rescued += 1
    return RescueReport(n=n, rescued=rescued, rescue_rate=rescued / n if n else 0.0)
