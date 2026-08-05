"""Retrospective journey replay.

On cases where the *resolving action* is known (eRepo/EAHAD curated trails;
Solve-RD/UDN as managed access clears), test whether the engine's #1 recommended
action matches what actually resolved the case - and at lower predicted cost.

Managed-access cohorts (Solve-RD/UDN) are NOT on the v1 critical path (long DAC
timelines); this harness runs on any (recommended, resolving-action, cost) trail.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Journey:
    case_id: str
    resolving_action: str           # the action that actually resolved the case
    recommended_top: str            # the engine's #1 recommended action
    recommended_rank_of_resolver: int   # rank the engine gave the true resolver (1-based)
    predicted_cost: float           # engine's predicted cost to resolution
    actual_cost: float              # the realized cost of the actual journey


@dataclass
class RetrospectiveReport:
    n: int
    top1_concordance: float         # fraction where #1 == actual resolver
    top3_concordance: float
    mean_predicted_cost: float
    mean_actual_cost: float
    cheaper_fraction: float         # fraction where predicted < actual


def replay(journeys: list[Journey]) -> RetrospectiveReport:
    n = len(journeys)
    if n == 0:
        return RetrospectiveReport(0, 0.0, 0.0, 0.0, 0.0, 0.0)
    top1 = sum(j.recommended_top == j.resolving_action for j in journeys) / n
    top3 = sum(1 <= j.recommended_rank_of_resolver <= 3 for j in journeys) / n
    return RetrospectiveReport(
        n=n, top1_concordance=top1, top3_concordance=top3,
        mean_predicted_cost=sum(j.predicted_cost for j in journeys) / n,
        mean_actual_cost=sum(j.actual_cost for j in journeys) / n,
        cheaper_fraction=sum(j.predicted_cost < j.actual_cost for j in journeys) / n,
    )
