"""Engine API.

Endpoints: POST /classify, POST /recommend, POST /case. The deterministic engine
computes everything; the LLM (if wired) only phrases explanations.

The core handler functions are pure-Python and unit-testable. ``build_app()`` imports
FastAPI lazily and wires the handlers, so this module imports fine in CI without
FastAPI installed (it is installed in the API image on the VM).
"""
from __future__ import annotations

from typing import Any

from core.schemas import (
    EvidenceContribution,
    PatientContext,
    PointsLedger,
    Strength,
    Variant,
)
from engine.recommend import recommend_next_action
from rules.acmg_codes import contributions_from_codes
from rules.point_engine import classify
from rules.posterior import posterior
from rules.vcep_loader import get_spec


def _ledger_from_payload(p: dict[str, Any]) -> tuple[PointsLedger, object]:
    v = Variant(p.get("chrom", "0"), int(p.get("pos", 0)), p.get("ref", "N"),
                p.get("alt", "N"), gene=p.get("gene", ""),
                hgvs_p=p.get("hgvs_p"), hgvs_c=p.get("hgvs_c"))
    spec = get_spec(v.gene)
    contribs: list[EvidenceContribution] = []
    if "codes" in p:                                  # explicit applied codes
        contribs, _ = contributions_from_codes(p["codes"])
    elif "contributions" in p:                        # explicit code/strength pairs
        for c in p["contributions"]:
            contribs.append(EvidenceContribution(c["code"], Strength[c["strength"]],
                                                 True, c.get("source", "api")))
    return PointsLedger(v, contributions=contribs, spec_version=spec.version), spec


def handle_classify(payload: dict[str, Any]) -> dict[str, Any]:
    ledger, spec = _ledger_from_payload(payload)
    cls = classify(ledger, spec)
    return {"gene": ledger.variant.gene, "variant": ledger.variant.key(),
            "points": ledger.points, "classification": cls.name,
            "posterior": round(posterior(ledger.points, spec.prior_p), 4),
            "spec_version": spec.version,
            "applied_codes": [c.code for c in ledger.contributions if c.applied]}


def handle_recommend(payload: dict[str, Any]) -> dict[str, Any]:
    ledger, spec = _ledger_from_payload(payload)
    patient = PatientContext(**payload.get("patient", {}))
    prefs = payload.get("prefs", {})
    rec = recommend_next_action(
        ledger, patient, mechanism=payload.get("mechanism", spec.mechanism),
        domain=payload.get("domain", "bleeding"), inheritance=spec.inheritance,
        prior_p=spec.prior_p, **{k: prefs[k] for k in
                                 ("objective", "alpha", "beta", "gamma", "tau") if k in prefs})
    return {
        "current_class": rec.current_class.name, "current_posterior": rec.current_posterior,
        "in_conflict": rec.in_conflict, "gap_to_lp": rec.gap_to_lp,
        "actions": rec.ranked, "pareto_frontier": rec.pareto_frontier,
        "case_note": rec.case_note, "explanation": rec.explanation,
    }


def handle_case(payload: dict[str, Any]) -> dict[str, Any]:
    """Whole-case routing decision (variant vs modality vs wait)."""
    from engine.case_policy import CaseState, decide_case
    st = CaseState(**{k: payload[k] for k in payload
                     if k in CaseState.__dataclass_fields__})
    plan = decide_case(st)
    return {"route": plan.route, "rationale": plan.rationale,
            "modality": plan.modality.name if plan.modality else None}


def handle_diagnose(payload: dict[str, Any]) -> dict[str, Any]:
    """DISCERN coupled disease x variant diagnosis."""
    from core.dx_schemas import Feature, FeatureKind
    from jointdx.factorgraph import Evidence
    from jointdx.orchestrate import diagnose

    def _feats(items, kind):
        return [Feature(f["id"], kind, f.get("value", True), observed=f.get("observed", True))
                for f in items]

    ev = Evidence(
        variant_gene=payload.get("gene", ""), variant_id=payload.get("variant_id", "V"),
        genetic_codes=payload.get("codes", []),
        clinical=_feats(payload.get("clinical", []), FeatureKind.CLINICAL),
        functional=_feats(payload.get("functional", []), FeatureKind.FUNCTIONAL))
    rec = diagnose(ev, planned_tx=payload.get("planned_tx"))
    if rec is None:
        return {"error": f"no discrimination cluster for gene {ev.variant_gene!r}"}
    return {
        "leading": rec.posterior.leading, "decided": rec.posterior.decided,
        "p_disease": rec.posterior.p_disease, "p_variant": rec.posterior.p_variant,
        "reclassified": rec.reclassified_variants,
        "safety_flags": [{"competitor": f.competitor_id, "severity": f.severity,
                          "message": f.message} for f in rec.safety_flags],
        "next_observation": rec.next_observation.name if rec.next_observation else None,
        "explanation": rec.explanation, "audit": rec.audit,
    }


def build_app():  # pragma: no cover - requires fastapi (API image)
    from fastapi import FastAPI
    app = FastAPI(title="DISCERN", version="0.0.1")

    @app.post("/classify")
    def _classify(payload: dict):
        return handle_classify(payload)

    @app.post("/recommend")
    def _recommend(payload: dict):
        return handle_recommend(payload)

    @app.post("/case")
    def _case(payload: dict):
        return handle_case(payload)

    @app.post("/diagnose")
    def _diagnose(payload: dict):
        return handle_diagnose(payload)

    @app.get("/health")
    def _health():
        return {"status": "ok"}

    return app
