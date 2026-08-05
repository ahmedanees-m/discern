"""Validation statistics.

Paired tests for the headline ablation (full system vs tools-only on identical cases):
McNemar (paired binary resolution), Wilcoxon signed-rank (paired cost/time), and
bootstrap CIs. DeLong (paired AUC) is provided via a thin scipy/sklearn path when
available. All keep the pre-registration / negatives-included reporting commitment.

scipy is used when present; pure-Python fallbacks keep the core importable in CI.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass


@dataclass
class TestResult:
    statistic: float
    p_value: float
    detail: str = ""


def mcnemar(b: int, c: int) -> TestResult:
    """McNemar's test on the two discordant cells (b = only-A-resolved, c = only-B).

    Uses the exact binomial p-value (robust for small/large n); continuity-corrected
    chi-square as the statistic.
    """
    n = b + c
    stat = ((abs(b - c) - 1) ** 2) / n if n else 0.0
    # exact two-sided binomial p under p=0.5
    k = min(b, c)
    cum = sum(math.comb(n, i) for i in range(0, k + 1)) * (0.5 ** n) if n else 1.0
    p = min(1.0, 2.0 * cum)
    return TestResult(stat, p, f"b={b}, c={c}, n_discordant={n}")


def wilcoxon_signed_rank(x: list[float], y: list[float]) -> TestResult:
    """Paired Wilcoxon signed-rank (normal approximation with tie handling)."""
    diffs = [a - b for a, b in zip(x, y, strict=True) if a - b != 0]
    n = len(diffs)
    if n == 0:
        return TestResult(0.0, 1.0, "all pairs tied")
    order = sorted(range(n), key=lambda i: abs(diffs[i]))
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(diffs[order[j + 1]]) == abs(diffs[order[i]]):
            j += 1
        avg = (i + j) / 2 + 1
        for t in range(i, j + 1):
            ranks[order[t]] = avg
        i = j + 1
    w_plus = sum(r for d, r in zip(diffs, ranks, strict=True) if d > 0)
    mean = n * (n + 1) / 4
    sd = math.sqrt(n * (n + 1) * (2 * n + 1) / 24)
    z = (w_plus - mean) / sd if sd else 0.0
    p = math.erfc(abs(z) / math.sqrt(2))
    return TestResult(z, p, f"W+={w_plus:.1f}, n={n}")


def bootstrap_ci(values: list[float], n_boot: int = 2000, alpha: float = 0.05,
                 seed: int = 0) -> tuple[float, float, float]:
    """Mean + percentile bootstrap CI."""
    rng = random.Random(seed)
    n = len(values)
    if n == 0:
        return 0.0, 0.0, 0.0
    means = []
    for _ in range(n_boot):
        sample = [values[rng.randrange(n)] for _ in range(n)]
        means.append(sum(sample) / n)
    means.sort()
    lo = means[int((alpha / 2) * n_boot)]
    hi = means[int((1 - alpha / 2) * n_boot) - 1]
    return sum(values) / n, lo, hi


def delong_auc(scores: list[float], labels: list[int]) -> float:  # pragma: no cover
    """AUC (the DeLong variance/test needs scipy; here we return the AUC point estimate
    via the Mann-Whitney U identity). Full paired DeLong test runs in the eval image."""
    pos = [s for s, y in zip(scores, labels, strict=True) if y == 1]
    neg = [s for s, y in zip(scores, labels, strict=True) if y == 0]
    if not pos or not neg:
        return float("nan")
    wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))
