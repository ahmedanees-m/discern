"""Mondrian split-conformal selective prediction over the diagnosis posterior (v3.1 Track B2).

Gives a per-class (Mondrian) coverage guarantee: when DISCERN commits, the true disease is in
the prediction set with probability >= 1 - alpha PER CLASS. Selective mode commits to the
leading disease only when the conformal set is the singleton {leading}; otherwise it abstains
and defers to the next-observation module (composing with jointdx/abstain.py).

Nonconformity score = 1 - P(disease | evidence). Calibration uses the standard split-conformal
quantile k = ceil((n+1)(1-alpha)). Coverage on REAL labels (H5 diagnosis half / Gate G11) is
cohort-gated; the synthetic-sampling sanity test verifies the machinery's guarantee, NOT
clinical accuracy (no synthetic result is reported as a headline).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from core.dx_schemas import DiscriminationCluster
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import leading_disease, marginal_disease


def nonconformity(cluster: DiscriminationCluster, ev: Evidence, label: str) -> float:
    return 1.0 - marginal_disease(joint(cluster, ev)).get(label, 0.0)


def calibrate(cluster: DiscriminationCluster, labeled: list[tuple[Evidence, str]],
              alpha: float = 0.1) -> dict[str, float]:
    """Mondrian per-class nonconformity thresholds from a calibration set of (Evidence, true_id)."""
    by_class: dict[str, list[float]] = {d.id: [] for d in cluster.diseases}
    for ev, y in labeled:
        by_class.setdefault(y, []).append(nonconformity(cluster, ev, y))
    q: dict[str, float] = {}
    for y, scores in by_class.items():
        if not scores:
            q[y] = 1.0
            continue
        scores.sort()
        n = len(scores)
        k = math.ceil((n + 1) * (1 - alpha))
        q[y] = scores[min(k, n) - 1]
    return q


def prediction_set(cluster: DiscriminationCluster, ev: Evidence, q: dict[str, float]) -> set[str]:
    md = marginal_disease(joint(cluster, ev))
    return {d.id for d in cluster.diseases if (1.0 - md.get(d.id, 0.0)) <= q.get(d.id, 1.0)}


@dataclass
class SelectiveResult:
    committed: bool
    leading: str
    pred_set: set
    reason: str


def selective_predict(cluster: DiscriminationCluster, ev: Evidence, q: dict[str, float]) -> SelectiveResult:
    lead, _ = leading_disease(joint(cluster, ev))
    ps = prediction_set(cluster, ev, q)
    if ps == {lead}:
        return SelectiveResult(True, lead, ps, "singleton conformal set -> commit")
    return SelectiveResult(False, lead, ps, f"conformal set size {len(ps)} -> abstain (defer to next-obs)")


def empirical_coverage(cluster: DiscriminationCluster, labeled: list[tuple[Evidence, str]],
                       q: dict[str, float]) -> dict:
    """Fraction of (Evidence, true_id) whose true label is in the prediction set, overall + per class."""
    hit = 0
    per: dict[str, list[int]] = {}
    for ev, y in labeled:
        ok = y in prediction_set(cluster, ev, q)
        hit += ok
        per.setdefault(y, [0, 0])
        per[y][0] += ok
        per[y][1] += 1
    n = len(labeled)
    return {"overall": hit / n if n else 0.0,
            "per_class": {y: (c / t if t else 0.0) for y, (c, t) in per.items()}}
