"""Core ablation - full system vs underlying-tools-only.

Paired comparison on identical cases: does the decision/orchestration layer add VUS
resolution and/or reduce cost/time over using the same tools WITHOUT it? Metrics +
paired statistics (McNemar on resolution, Wilcoxon on cost/time). **Pre-register
before running**; report negatives.
"""
from __future__ import annotations

from dataclasses import dataclass

from eval.stats import bootstrap_ci, mcnemar, wilcoxon_signed_rank


@dataclass
class PairedCase:
    full_resolved: bool          # full system resolved this case
    base_resolved: bool          # tools-only resolved it
    full_cost: float
    base_cost: float


@dataclass
class AblationReport:
    n: int
    full_resolution_rate: float
    base_resolution_rate: float
    mcnemar_p: float
    cost_wilcoxon_p: float
    mean_cost_delta: float
    cost_delta_ci: tuple[float, float, float]


def run_ablation(cases: list[PairedCase]) -> AblationReport:
    n = len(cases)
    b = sum(c.full_resolved and not c.base_resolved for c in cases)   # only full
    c_ = sum(c.base_resolved and not c.full_resolved for c in cases)  # only base
    mc = mcnemar(b, c_)
    full_costs = [c.full_cost for c in cases]
    base_costs = [c.base_cost for c in cases]
    wil = wilcoxon_signed_rank(full_costs, base_costs)
    deltas = [c.full_cost - c.base_cost for c in cases]
    return AblationReport(
        n=n,
        full_resolution_rate=sum(c.full_resolved for c in cases) / n if n else 0.0,
        base_resolution_rate=sum(c.base_resolved for c in cases) / n if n else 0.0,
        mcnemar_p=mc.p_value,
        cost_wilcoxon_p=wil.p_value,
        mean_cost_delta=sum(deltas) / n if n else 0.0,
        cost_delta_ci=bootstrap_ci(deltas),
    )
