"""Equity reporting (plan Phase 3, Step 3.2).

Per-case and cohort-level metrics quantifying how much of a classification rests on
ancestry-biased vs ancestry-equitable evidence, and the resolution-rate gap across
genetic-similarity groups (target: the equity-aware config narrows it).
"""
from __future__ import annotations

from dataclasses import dataclass

from core.schemas import PointsLedger

_BIASED = {"PM2", "BS1", "BA1", "PP3", "BP4"}


@dataclass
class CaseEquity:
    biased_points: float
    equitable_points: float
    biased_fraction: float       # share of |points| resting on ancestry-biased evidence


def case_equity(ledger: PointsLedger) -> CaseEquity:
    biased = equit = 0.0
    for c in ledger.contributions:
        if not c.applied:
            continue
        mag = abs(c.strength.value * c.reliability)
        if c.code.split("_", 1)[0] in _BIASED:
            biased += mag
        else:
            equit += mag
    total = biased + equit
    return CaseEquity(biased, equit, (biased / total) if total else 0.0)


def cohort_resolution_gap(resolved_by_group: dict[str, tuple[int, int]]) -> dict:
    """resolved_by_group: group -> (n_resolved, n_total). Returns per-group rates and the
    max-min gap (the disparity the equity-aware routing aims to reduce)."""
    rates = {g: (r / t if t else 0.0) for g, (r, t) in resolved_by_group.items()}
    gap = (max(rates.values()) - min(rates.values())) if rates else 0.0
    return {"rates": rates, "resolution_gap": gap}
