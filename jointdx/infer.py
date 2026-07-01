"""Inference over the joint posterior (plan Phase 3.1 / 3.2).

Marginals P(D|E)=Sum_V and P(V|E)=Sum_D from the joint, the leading disease with its
confidence, and VUS reclassification (old->new state with drivers) - the measured,
falsifiable VUS-mitigation output.
"""
from __future__ import annotations

from core.dx_schemas import VariantState


def marginal_disease(joint: dict) -> dict[str, float]:
    out: dict[str, float] = {}
    for (d_id, _v), p in joint.items():
        out[d_id] = out.get(d_id, 0.0) + p
    return out


def marginal_variant(joint: dict) -> dict[VariantState, float]:
    out: dict[VariantState, float] = {v: 0.0 for v in VariantState}
    for (_d, v), p in joint.items():
        out[v] += p
    return out


def leading_disease(joint: dict) -> tuple[str, float]:
    md = marginal_disease(joint)
    if not md:
        return "", 0.0
    d_id = max(md, key=md.get)
    return d_id, md[d_id]


def variant_call(joint: dict) -> tuple[VariantState, float]:
    mv = marginal_variant(joint)
    v = max(mv, key=mv.get)
    return v, mv[v]


# Collapse the 5-state marginal to a P/LP/VUS/LB/B call by the dominant + neighbor mass.
_ORDER = [VariantState.BEN, VariantState.LB, VariantState.VUS, VariantState.LP, VariantState.PATH]


def reclassify(joint: dict, old_state: VariantState = VariantState.VUS,
               decide_threshold: float = 0.5) -> tuple[VariantState, VariantState, dict]:
    """Return (old, new, drivers). `new` is the argmax variant state if it clears the
    threshold, else unchanged (stays the old state - non-reclassification)."""
    mv = marginal_variant(joint)
    new_state, p_new = max(mv.items(), key=lambda kv: kv[1])
    drivers = {"p_variant": {s.name: round(p, 4) for s, p in mv.items()},
               "leading_disease": leading_disease(joint)}
    if p_new < decide_threshold:
        return old_state, old_state, drivers          # not confident enough to reclassify
    return old_state, new_state, drivers
