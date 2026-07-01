"""Calibration - reliability diagram + Expected Calibration Error (plan Section A.5).

The engine's *promises* (predicted posterior / predicted yield) must be calibrated,
not just its classifications. Bins predicted probabilities and compares to observed
frequencies; ECE is the weighted mean |confidence - accuracy| across bins.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class CalibrationBin:
    lo: float
    hi: float
    n: int
    mean_pred: float
    observed: float


def reliability_curve(pred: list[float], outcome: list[int], n_bins: int = 10):
    bins: list[CalibrationBin] = []
    for i in range(n_bins):
        lo, hi = i / n_bins, (i + 1) / n_bins
        idx = [j for j, p in enumerate(pred) if (lo <= p < hi) or (i == n_bins - 1 and p == 1.0)]
        if not idx:
            bins.append(CalibrationBin(lo, hi, 0, 0.0, 0.0))
            continue
        mp = sum(pred[j] for j in idx) / len(idx)
        obs = sum(outcome[j] for j in idx) / len(idx)
        bins.append(CalibrationBin(lo, hi, len(idx), mp, obs))
    return bins


def expected_calibration_error(pred: list[float], outcome: list[int], n_bins: int = 10) -> float:
    n = len(pred)
    if n == 0:
        return 0.0
    return sum(b.n / n * abs(b.mean_pred - b.observed)
               for b in reliability_curve(pred, outcome, n_bins))
