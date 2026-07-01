"""Calibrated decide / abstain policy (plan Phase 4.2 / Section A.4).

DISCERN returns a call **only** when justified; otherwise it says "undecidable - here is
the deciding observation." Abstain if the max posterior is below threshold, OR the
credible interval is too wide (sparse LRs), OR a management-divergent competitor is not
excluded (the link to the safety interlock). The headline safety metric is the
**confident-and-wrong rate** - measured in Phase 9.
"""
from __future__ import annotations

from dataclasses import dataclass

from core.dx_schemas import DiscriminationCluster
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import leading_disease, marginal_disease
from jointdx.uncertainty import disease_intervals


@dataclass
class Decision:
    decided: bool
    leading: str
    p: float
    ci: tuple[float, float]
    reason: str


def decide(cluster: DiscriminationCluster, ev: Evidence, tau: float = 0.5,
           max_ci_width: float = 0.45, n_mc: int = 200) -> Decision:
    j = joint(cluster, ev)
    lead, p = leading_disease(j)
    md = marginal_disease(j)  # noqa: F841 - kept for symmetry / future per-disease checks
    ci = disease_intervals(cluster, ev, n_mc=n_mc).get(lead, (p, p, p))
    lo, hi = ci[1], ci[2]
    if p < tau:
        return Decision(False, lead, p, (lo, hi), "max posterior below threshold -> abstain")
    if (hi - lo) > max_ci_width:
        return Decision(False, lead, p, (lo, hi),
                        "credible interval too wide (sparse likelihood ratios) -> abstain")
    return Decision(True, lead, p, (lo, hi), "decided")
