"""The statistical primitives behind the reported numbers, in one tested place.

These functions compute the paper's headline statistics: the Wilson interval on the partition's
inflation rate, DeLong's test on the AUROC difference, expected calibration error, percentile
bootstrap intervals, exact McNemar on paired case correctness, and rank metrics for the external
comparison. They previously lived as private helpers inside the analysis harnesses, where they
could not be exercised by continuous integration - the harnesses need scikit-learn, pandas and
several gigabytes of third-party data, none of which CI has. A silent arithmetic error in any of
them would have propagated straight into the manuscript.

Everything here depends only on numpy and scipy, so it runs everywhere, and every function is
covered by tests/test_stats.py against closed-form or hand-computed values.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.stats import binomtest, norm


# ---------------------------------------------------------------------------- proportions
def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion.

    Preferred over the normal approximation because the rates reported here sit far from 0.5 with
    large n, where Wald intervals misbehave.
    """
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, centre - half), min(1.0, centre + half))


# ---------------------------------------------------------------------------- calibration
def ece(probs, labels, bins: int = 10) -> float:
    """Expected calibration error: mean |confidence - accuracy| weighted by bin occupancy."""
    probs, labels = np.asarray(probs, float), np.asarray(labels, float)
    edges = np.linspace(0, 1, bins + 1)
    total = 0.0
    for i in range(bins):
        hi = edges[i + 1] if i < bins - 1 else 1.0001
        m = (probs >= edges[i]) & (probs < hi)
        if m.sum():
            total += (m.sum() / len(probs)) * abs(labels[m].mean() - probs[m].mean())
    return float(total)


# ---------------------------------------------------------------------------- resampling
def bootstrap_indices(n: int, rng, n_boot: int = 1000):
    """Row indices for `n_boot` bootstrap resamples of `n` observations."""
    return rng.integers(0, n, size=(n_boot, n))


def percentile_ci(values, lo: float = 2.5, hi: float = 97.5, ndigits: int = 4):
    """Percentile interval over a bootstrap distribution, ignoring non-finite draws."""
    v = np.asarray([x for x in values if x is not None and np.isfinite(x)], float)
    if not len(v):
        return None
    return [round(float(np.percentile(v, lo)), ndigits),
            round(float(np.percentile(v, hi)), ndigits)]


# ---------------------------------------------------------------------------- ROC comparison
def midrank(x):
    """Ranks with ties averaged - the tie handling DeLong's covariance estimate requires."""
    order = np.argsort(x)
    z = np.asarray(x, float)[order]
    n = len(z)
    t = np.zeros(n, float)
    i = 0
    while i < n:
        j = i
        while j < n and z[j] == z[i]:
            j += 1
        t[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    out = np.empty(n, float)
    out[order] = t
    return out


def delong(labels, score_a, score_b) -> dict:
    """DeLong's test for two correlated ROC curves on the same observations.

    Fast form of Sun and Xu (2014). Returns both AUROCs, their difference, the z statistic and a
    two-sided p value. `z` and `p_value` are None when the estimated variance of the difference is
    not positive, which happens when the two scores are identical.
    """
    labels = np.asarray(labels, int)
    pos = labels == 1
    if pos.sum() == 0 or (~pos).sum() == 0:
        raise ValueError("DeLong needs both classes present")
    mat = np.vstack([np.concatenate([np.asarray(s, float)[pos], np.asarray(s, float)[~pos]])
                     for s in (score_a, score_b)])
    m = int(pos.sum())
    n = mat.shape[1] - m
    tx = np.vstack([midrank(mat[r, :m]) for r in range(2)])
    ty = np.vstack([midrank(mat[r, m:]) for r in range(2)])
    tz = np.vstack([midrank(mat[r, :]) for r in range(2)])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    cov = np.cov(v01) / m + np.cov(v10) / n
    var_diff = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    base = {"auroc_a": round(float(aucs[0]), 4), "auroc_b": round(float(aucs[1]), 4),
            "delta": round(float(aucs[0] - aucs[1]), 4)}
    if var_diff <= 0:
        return {**base, "z": None, "p_value": None}
    z = float((aucs[0] - aucs[1]) / np.sqrt(var_diff))
    return {**base, "z": round(z, 4), "p_value": round(float(2 * norm.sf(abs(z))), 6)}


# ---------------------------------------------------------------------------- paired outcomes
def mcnemar_exact(a_correct, b_correct) -> dict:
    """Exact McNemar on paired binary correctness. `a` is the system under test."""
    a_only = sum(1 for a, b in zip(a_correct, b_correct, strict=True) if a and not b)
    b_only = sum(1 for a, b in zip(a_correct, b_correct, strict=True) if b and not a)
    n = a_only + b_only
    p = float(binomtest(a_only, n, 0.5).pvalue) if n else 1.0
    return {"a_only_correct": a_only, "b_only_correct": b_only,
            "discordant_pairs": n, "p_value_exact": round(p, 6)}


# ---------------------------------------------------------------------------- rank metrics
def rank_metrics(hit_ranks, ks=(1, 3, 5, 10)) -> dict:
    """Recall@k and mean reciprocal rank from 1-based hit positions; None means no hit."""
    n = len(hit_ranks)
    if not n:
        return {"n": 0}
    out = {"n": n}
    for k in ks:
        out[f"recall@{k}"] = round(sum(1 for h in hit_ranks if h is not None and h <= k) / n, 4)
    out["mrr"] = round(sum(1.0 / h for h in hit_ranks if h is not None) / n, 4)
    return out
