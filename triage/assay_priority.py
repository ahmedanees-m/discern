"""Scientist-facing VUS-triage.

Of all the VUS in a cohort, which - if functionally assayed - would resolve the **most**
diagnoses (highest expected information gain toward resolution)? This turns DISCERN into a
wet-lab prioritization engine. Most Tier-1 genes lack DMS (not growth-selectable), so
per-patient functional tests carry the load - making "which to assay next" genuinely
valuable.
"""
from __future__ import annotations

from dataclasses import dataclass

from diseases.ontology import route_clusters
from jointdx.factorgraph import Evidence, joint
from jointdx.orchestrate import _best_cluster
from nextobs.recommend import eig_joint


@dataclass
class TriageItem:
    variant_id: str
    gene: str
    eig: float
    best_assay: str
    changes_management: bool


def prioritize(cohort: list[Evidence]) -> list[TriageItem]:
    """Rank VUS by the information gain of their best functional/lab assay (desc)."""
    items: list[TriageItem] = []
    for ev in cohort:
        clusters = route_clusters([ev.variant_gene]) if ev.variant_gene else []
        if not clusters:
            continue
        cl = _best_cluster(clusters, ev)
        j = joint(cl, ev)
        assays = [o for o in cl.next_observations if o.kind in ("functional", "lab")]
        scored = [(o, eig_joint(j, o)) for o in assays]
        if not scored:
            continue
        best, eig = max(scored, key=lambda t: t[1])
        items.append(TriageItem(ev.variant_id, ev.variant_gene, round(eig, 4),
                                best.name, best.changes_management))
    return sorted(items, key=lambda x: x.eig, reverse=True)
