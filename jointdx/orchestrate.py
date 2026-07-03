"""End-to-end orchestration - one auditable diagnose() call (plan Phase 8.2).

route -> joint -> abstain -> safety flags -> next observation -> assemble. Returns a
``DxRecommendation`` with the ranked diagnosis (credible intervals), the measured VUS
reclassification, the management-aware safety flags, and the cheapest decisive next
observation - all deterministic, with a templated explanation and an audit dict.
"""
from __future__ import annotations

from core.dx_schemas import DxRecommendation, JointPosterior, VariantState
from diseases.ontology import route_clusters
from jointdx.abstain import decide
from jointdx.explain import explain
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import leading_disease, marginal_variant, reclassify
from jointdx.uncertainty import disease_intervals
from nextobs.recommend import best_observation
from safety.interlock import flags


def _best_cluster(clusters, ev: Evidence):
    """When a gene maps to more than one cluster, pick the one that best fits the evidence.

    Routing keys on the PRESENT (positive) findings the patient actually exhibits, not on
    every feature the user was asked about. A pertinent negative (an explicitly *absent*
    feature) must not pull the case toward the cluster that happens to define it: ruling out
    the macrothrombocytopenia findings for a GP1BA case whose phenotype is enhanced-RIPA
    must route to the platelet-type-VWD cluster, not to the macrothrombocytopenia/ITP
    cluster. Ties break on the leading-disease posterior, which already reflects the full
    evidence (present and absent)."""
    present_ids = ({f.id for f in ev.clinical if f.observed}
                   | {f.id for f in ev.functional if f.observed})

    def score(c):
        known = {fid for d in c.diseases for fid in d.feature_lr}
        return (len(present_ids & known), leading_disease(joint(c, ev))[1])

    return max(clusters, key=score)


def diagnose(ev: Evidence, planned_tx: str | None = None,
             old_state: VariantState = VariantState.VUS, decide_tau: float = 0.5,
             n_mc: int = 150) -> DxRecommendation | None:
    clusters = route_clusters([ev.variant_gene]) if ev.variant_gene else []
    if not clusters:
        return None
    cluster = _best_cluster(clusters, ev)      # gene may map to >1 cluster -> best evidence fit

    j = joint(cluster, ev)
    dec = decide(cluster, ev, tau=decide_tau, n_mc=n_mc)
    intervals = disease_intervals(cluster, ev, n_mc=n_mc)
    p_disease = {d.id: (round(intervals[d.id][0], 4), round(intervals[d.id][1], 4),
                        round(intervals[d.id][2], 4)) for d in cluster.diseases}
    p_variant = {ev.variant_id: {s.name: round(p, 4) for s, p in marginal_variant(j).items()}}

    posterior = JointPosterior(
        cluster=cluster, p_disease=p_disease, p_variant=p_variant,
        leading=dec.leading, confidence=round(dec.p, 4), decided=dec.decided)

    safety = flags(cluster, j, planned_tx=planned_tx)
    next_obs = best_observation(cluster, ev)
    old, new, drivers = reclassify(j, old_state)
    text = explain(cluster, p_disease, (old, new, drivers), safety, next_obs, dec.decided)

    return DxRecommendation(
        posterior=posterior,
        reclassified_variants={ev.variant_id: (old.name, new.name, drivers)},
        safety_flags=safety, next_observation=next_obs, explanation=text,
        audit={"cluster": cluster.id, "decided": dec.decided, "reason": dec.reason,
               "vcep_specs": {d.id: d.vcep_spec for d in cluster.diseases},
               "evidence": {"gene": ev.variant_gene, "codes": ev.genetic_codes,
                            "n_clinical": len(ev.clinical), "n_functional": len(ev.functional)}})
