"""Interactive what-if (plan Phase 6c).

Shows how the ranking shifts *before* ordering a test - e.g. "RIPA platelet-origin ->
PT-VWD 95%; plasma-origin -> 2B 90%" - so the clinician sees the value of each outcome
before committing.
"""
from __future__ import annotations

from core.dx_schemas import DiscriminationCluster
from jointdx.factorgraph import Evidence
from nextobs.recommend import predicted_shift


def whatif(cluster: DiscriminationCluster, ev: Evidence, observation_id: str) -> dict:
    """For the named observation, return {outcome: (leading_disease, prob)}."""
    obs = next((o for o in cluster.next_observations if o.id == observation_id), None)
    if obs is None:
        return {}
    return predicted_shift(cluster, ev, obs)


def whatif_text(cluster: DiscriminationCluster, ev: Evidence, observation_id: str) -> str:
    shifts = whatif(cluster, ev, observation_id)
    parts = [f"{outcome} -> {lead} {p:.0%}" for outcome, (lead, p) in shifts.items()]
    return "; ".join(parts)
