"""Management-aware misdiagnosis / treatment-safety interlock.

Fires on **treatment danger**, not posterior gap: a small probability of a
*treatment-changing* competitor fires; a large probability of a management-irrelevant one
does not. Plus a hard-stop interlock when a planned treatment is contraindicated by a
non-excluded competitor (DDAVP + 2B/PT-VWD; splenectomy + BSS; HSCT + LAD-III).
"""
from __future__ import annotations

from core.dx_schemas import DiscriminationCluster, Disease, SafetyFlag
from jointdx.infer import leading_disease, marginal_disease

# treatments whose mis-application is high-severity (transplant / irreversible / harmful).
_HIGH_SEVERITY_TX = {"hsct", "splenectomy", "corticosteroids", "corticosteroids_ivig_splenectomy",
                     "affected_related_donor_transplant"}


def treatment_changes(a: Disease, b: Disease) -> bool:
    return a.treatment != b.treatment


def severity(a: Disease, b: Disease) -> float:
    """1.0 when confusing a<->b risks a transplant-grade or harmful intervention."""
    if "hsct" in (a.treatment, b.treatment):
        return 1.0
    if any(t in _HIGH_SEVERITY_TX for t in (a.treatment, b.treatment)):
        return 0.9
    return 0.6


def _resolving_obs(cluster: DiscriminationCluster):
    for o in cluster.next_observations:
        if o.changes_management:
            return o
    return cluster.next_observations[0] if cluster.next_observations else None


def flags(cluster: DiscriminationCluster, joint: dict, planned_tx: str | None = None,
          tau: float = 0.05, lead_id: str | None = None) -> list[SafetyFlag]:
    """Safety flags for a leading call, judged against the competitors this joint leaves open.

    `lead_id` lets the caller pin the leading diagnosis to the engine's actual call while passing a
    gene-blind joint for the competitor probabilities, so the flags are always reported against the
    diagnosis the clinician was shown and a competitor is only retired by evidence that excludes it.
    """
    md = marginal_disease(joint)
    if lead_id is None:
        lead_id, _ = leading_disease(joint)
    by_id = {d.id: d for d in cluster.diseases}
    d_star = by_id.get(lead_id)
    out: list[SafetyFlag] = []
    if d_star is None:
        return out
    for d in cluster.diseases:
        p = md.get(d.id, 0.0)
        # management-aware divergence flag: competitor-only (the leading call is the baseline).
        if d.id != lead_id:
            sev = severity(d_star, d)
            if treatment_changes(d_star, d) and p * sev >= tau:
                out.append(SafetyFlag(
                    leading_id=lead_id, competitor_id=d.id,
                    management_divergence=f"{d_star.treatment} vs {d.treatment}",
                    p_competitor=round(p, 4), severity="high" if sev >= 0.9 else "moderate",
                    resolving_observation=_resolving_obs(cluster),
                    message=(f"If {d.name} (p={p:.2f}) rather than {d_star.name}, management changes "
                             f"({d_star.treatment} -> {d.treatment}). Resolve before treating.")))
        # A planned treatment contraindicated by ANY non-excluded disease fires the hard stop,
        # the leading call included: a contraindication does not stop applying because the disease
        # that carries it happens to be the one ranked first.
        if planned_tx and planned_tx in d.contraindications and p > 0:
            is_lead = d.id == lead_id
            out.append(SafetyFlag(
                leading_id=lead_id, competitor_id=d.id,
                management_divergence=f"{planned_tx} contraindicated if {d.name}",
                p_competitor=round(p, 4), severity="high",
                resolving_observation=_resolving_obs(cluster),
                message=(f"HARD STOP: {planned_tx} is contraindicated if {d.name} (p={p:.2f})"
                         + (" [leading diagnosis]" if is_lead else "")
                         + " - resolve with the deciding observation first.")))
    return out
