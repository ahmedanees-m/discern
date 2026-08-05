"""Recommendation assembly.

Ties the decision core together: from a points ledger + patient context + the gene's
disease mechanism, it (1) finds attainable evidence codes, (2) maps them to concrete
mechanism-correct actions, (3) ranks by value-of-information per cost - escalating to
lookahead on ties, applying a risk gate, and (for VUS-by-conflict) preferring an
orthogonal disambiguating assay - and (4) emits an explained, audited recommendation
with the *expected post-action classification* and a cost-EIG Pareto frontier.

The explanation is templated from the deterministic plan - an LLM (via llm/gateway)
may phrase it, but never decides it.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from core.schemas import Action, Classification, PatientContext, PointsLedger
from engine.action_map import expand_actions
from engine.gap import (
    ORTHOGONAL_MODALITIES,
    attainable_codes,
    case_note,
    detect_conflict,
    gap_to_target,
)
from engine.voi import lookahead, pareto_frontier, rank_actions, risk_gate
from equity.routing import equity_filter
from rules.point_engine import BANDS, classify
from rules.posterior import posterior


@dataclass
class Recommendation:
    current_class: Classification
    current_posterior: float
    in_conflict: bool
    gap_to_lp: float
    ranked: list[dict]
    pareto_frontier: list[dict] = field(default_factory=list)
    case_note: str = ""
    explanation: str = ""
    audit_index: int | None = None


def _class_for_points(points: float) -> Classification:
    for threshold, cls in BANDS:
        if points >= threshold:
            return cls
    return Classification.B


def _expected_post_class(points: float, action: Action) -> tuple[Classification, float]:
    """Best-case (positive-result) class + posterior if this action yields its top code."""
    k = 4.0 if any(c.startswith("PS") for c in action.yields_codes) else 2.0
    pts = points + k
    return _class_for_points(pts), posterior(pts)


def recommend_next_action(ledger: PointsLedger, patient: PatientContext, mechanism: str,
                          domain: str = "bleeding", prior_p: float = 0.10,
                          objective: str = "delta_utility", alpha: float = 1.0,
                          beta: float = 50.0, gamma: float = 0.0, tau: float = 0.0,
                          inheritance: str = "", target: Classification = Classification.LP,
                          tie_eps: float = 0.05, equity_routing: bool = True,
                          audit=None) -> Recommendation:
    points = ledger.points
    current_class = classify(ledger)
    conflict = detect_conflict(ledger)

    opps = attainable_codes(ledger, patient, mechanism, domain)
    codes = [o.code for o in opps]
    actions = expand_actions(codes, mechanism, domain)

    # VUS-by-conflict: prefer orthogonal disambiguating evidence, not point-stacking (Section 6).
    if conflict.in_conflict:
        orth = [a for a in actions if a.modality in ORTHOGONAL_MODALITIES]
        if orth:
            actions = orth

    # Equity routing: for under-represented-ancestry patients, route to
    # ancestry-robust evidence (functional/segregation) over ancestry-biased codes.
    if equity_routing:
        actions = equity_filter(actions, patient)

    ranked = rank_actions(points, actions, prior_p=prior_p, objective=objective,
                          alpha=alpha, beta=beta, gamma=gamma)
    ranked = risk_gate(ranked, tau)

    # Lookahead escalation when the top two are within epsilon (Section 3.8).
    escalated = False
    if len(ranked) >= 2 and ranked[0][1] > 0:
        if abs(ranked[0][1] - ranked[1][1]) <= tie_eps * ranked[0][1]:
            top_actions = [a for a, _, _ in ranked[: min(6, len(ranked))]]
            _, seq = lookahead(points, top_actions, depth=2, prior_p=prior_p)
            if seq and ranked[0][0] is not seq[0]:
                ranked.sort(key=lambda t: (t[0] is not seq[0], -t[1]))
                escalated = True

    rows = []
    for a, vd, s in ranked:
        post_cls, post_p = _expected_post_class(points, a)
        rows.append({
            "action": a.name, "modality": a.modality, "yields": a.yields_codes,
            "cost_usd": a.cost_usd, "turnaround_days": a.turnaround_days,
            "value_density": round(vd, 6), "eig_bits": round(s["eig"], 4),
            "delta_utility": round(s["delta_utility"], 4),
            "expected_post_class": post_cls.name, "expected_post_posterior": round(post_p, 3),
        })

    pareto = pareto_frontier(points, actions, prior_p=prior_p)
    note = case_note(inheritance)
    explanation = _explain(points, current_class, rows, gap_to_target(ledger, target),
                           target, conflict.in_conflict, escalated, note)

    audit_index = None
    if audit is not None:
        entry = audit.append("recommend", {
            "gene": ledger.variant.gene, "variant": ledger.variant.key(),
            "points": points, "current_class": current_class.name,
            "in_conflict": conflict.in_conflict, "top_action": rows[0]["action"] if rows else None,
            "spec_version": ledger.spec_version,
        })
        audit_index = entry.index

    return Recommendation(
        current_class=current_class,
        current_posterior=round(posterior(points, prior_p), 3),
        in_conflict=conflict.in_conflict,
        gap_to_lp=gap_to_target(ledger, target),
        ranked=rows,
        pareto_frontier=pareto,
        case_note=note,
        explanation=explanation,
        audit_index=audit_index,
    )


def _explain(points, current_class, rows, gap, target, in_conflict, escalated, note) -> str:
    if not rows:
        return (f"{current_class.value} at {points:+.0f} points. No attainable evidence "
                "for this patient/mechanism - recommend reanalysis or matchmaking.")
    top = rows[0]
    head = "VUS-by-conflict" if in_conflict else current_class.value
    msg = (f"{head} at {points:+.0f} points ({gap:.0f} to {target.name}). "
           f"Recommended: {top['action']} (~${top['cost_usd']:.0f}, "
           f"~{top['turnaround_days']}d) -> expected {'/'.join(top['yields'])} -> "
           f"posterior {top['expected_post_posterior']:.2f} -> {top['expected_post_class']}.")
    if in_conflict:
        msg += " (Orthogonal disambiguating assay chosen to resolve conflicting evidence.)"
    if escalated:
        msg += " (Lookahead used to break a near-tie.)"
    if len(rows) > 1:
        msg += f" Alternative: {rows[1]['action']} (~${rows[1]['cost_usd']:.0f})."
    if note:
        msg += f" Note: {note}"
    return msg
