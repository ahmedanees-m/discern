"""Cheapest decisive next observation - usefulness-load-bearing.

Ranks candidate observations (lab / functional **+ segregation + phasing**) by expected
information gain over the JOINT P(D,V|E) - one score that values disease-resolving and
variant-upgrading steps together - then by `changes_management`, with accessibility as a
tiebreak. No currency in the primary ranking (the equity choice).
"""
from __future__ import annotations

import math

from core.dx_schemas import DiscriminationCluster, Observation, VariantState
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import leading_disease

_PATH_LIKE = {VariantState.PATH: 1.0, VariantState.LP: 0.8, VariantState.VUS: 0.5,
              VariantState.LB: 0.2, VariantState.BEN: 0.0}
_ACCESS_RANK = {"high": 2, "moderate": 1, "low": 0}


def _entropy(dist: dict) -> float:
    return -sum(p * math.log(p) for p in dist.values() if p > 0)


def _cell_factor(d_id: str, v: VariantState, lr_map: dict) -> float:
    f = 1.0
    if d_id in lr_map:
        f *= float(lr_map[d_id])
    if "V" in lr_map:
        f *= float(lr_map["V"]) ** (2 * _PATH_LIKE[v] - 1)
    return f


def eig_joint(j: dict, obs: Observation) -> float:
    """Expected information gain of `obs` over the joint posterior `j`."""
    if not obs.outcome_lr:
        return 0.0
    h0 = _entropy(j)
    outcomes = list(obs.outcome_lr.keys())
    # P(o | cell) proportional-to factor(cell, o), normalized over outcomes per cell.
    p_o = dict.fromkeys(outcomes, 0.0)
    post = {o: {} for o in outcomes}
    for cell, pcell in j.items():
        d_id, v = cell
        facs = {o: _cell_factor(d_id, v, obs.outcome_lr[o]) for o in outcomes}
        z = sum(facs.values()) or 1.0
        for o in outcomes:
            w = pcell * facs[o] / z
            p_o[o] += w
            post[o][cell] = post[o].get(cell, 0.0) + w
    eig = h0
    for o in outcomes:
        if p_o[o] <= 0:
            continue
        norm = {c: w / p_o[o] for c, w in post[o].items()}
        eig -= p_o[o] * _entropy(norm)
    return max(eig, 0.0)


def recommend(cluster: DiscriminationCluster, ev: Evidence) -> list[tuple[Observation, float]]:
    """Ranked observations: (changes_management, EIG, accessibility) desc."""
    j = joint(cluster, ev)
    scored = [(o, eig_joint(j, o)) for o in cluster.next_observations]
    scored.sort(key=lambda t: (t[0].changes_management, round(t[1], 6),
                               _ACCESS_RANK.get(t[0].accessibility, 1)), reverse=True)
    return scored


def best_observation(cluster: DiscriminationCluster, ev: Evidence) -> Observation | None:
    r = recommend(cluster, ev)
    return r[0][0] if r else None


def predicted_shift(cluster: DiscriminationCluster, ev: Evidence, obs: Observation) -> dict:
    """For each outcome, the resulting leading disease + its probability (for the explanation)."""
    j = joint(cluster, ev)
    out = {}
    for o, lr_map in obs.outcome_lr.items():
        post = {}
        for cell, pcell in j.items():
            d_id, v = cell
            post[cell] = pcell * _cell_factor(d_id, v, lr_map)
        z = sum(post.values()) or 1.0
        post = {c: p / z for c, p in post.items()}
        lead, p = leading_disease(post)
        out[o] = (lead, round(p, 3))
    return out
