"""Tests for the statistics that produce the manuscript's numbers.

These functions used to live as private helpers inside analysis harnesses that continuous
integration cannot execute - they need scikit-learn and several gigabytes of third-party data. That
left the arithmetic behind every reported confidence interval, p value and calibration figure
unexercised. Each is checked here against a closed form, a hand-computed value, or an invariant
that must hold regardless of implementation.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from core.stats import (
    bootstrap_indices,
    delong,
    ece,
    mcnemar_exact,
    midrank,
    percentile_ci,
    rank_metrics,
    wilson,
)


# ---------------------------------------------------------------- wilson
def test_wilson_matches_the_closed_form():
    """Hand-computed against the published formula for a mid-range proportion."""
    lo, hi = wilson(50, 100)
    z, n, p = 1.96, 100, 0.5
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    assert lo == pytest.approx(centre - half, abs=1e-12)
    assert hi == pytest.approx(centre + half, abs=1e-12)


def test_wilson_reproduces_the_reported_partition_interval():
    """The 33.2% inflation figure and its interval, exactly as the manuscript states them."""
    lo, hi = wilson(4068, 12240)
    assert round(4068 / 12240, 3) == 0.332
    assert (round(lo, 3), round(hi, 3)) == (0.324, 0.341)


def test_wilson_is_bounded_and_brackets_the_estimate():
    for k, n in ((0, 10), (10, 10), (1, 1000), (999, 1000), (3, 7)):
        lo, hi = wilson(k, n)
        assert 0.0 <= lo <= hi <= 1.0
        assert lo <= k / n <= hi
    assert wilson(0, 0) == (0.0, 0.0)


def test_wilson_interval_narrows_with_n():
    wide = wilson(5, 10)
    narrow = wilson(500, 1000)
    assert (narrow[1] - narrow[0]) < (wide[1] - wide[0])


# ---------------------------------------------------------------- ece
def test_ece_is_zero_for_a_perfectly_calibrated_predictor():
    # in each bin the predicted probability equals the observed frequency
    probs = [0.25] * 100 + [0.75] * 100
    labels = [1] * 25 + [0] * 75 + [1] * 75 + [0] * 25
    assert ece(probs, labels) == pytest.approx(0.0, abs=1e-12)


def test_ece_equals_the_gap_for_a_maximally_wrong_predictor():
    assert ece([1.0] * 50, [0] * 50) == pytest.approx(1.0)
    assert ece([0.0] * 50, [1] * 50) == pytest.approx(1.0)


def test_ece_weights_bins_by_occupancy():
    # 90 well-calibrated observations at 0.5, 10 badly wrong at 0.95 -> 0.10 * 0.95
    probs = [0.5] * 90 + [0.95] * 10
    labels = [1] * 45 + [0] * 45 + [0] * 10
    assert ece(probs, labels) == pytest.approx(0.10 * 0.95, abs=1e-9)


def test_ece_includes_probability_one_in_the_last_bin():
    """A boundary the two earlier private copies handled differently; both must include 1.0."""
    assert ece([1.0, 1.0], [1, 1]) == pytest.approx(0.0)
    assert ece([1.0, 1.0], [0, 0]) == pytest.approx(1.0)


# ---------------------------------------------------------------- midrank
def test_midrank_without_ties_is_the_plain_rank():
    assert list(midrank([10.0, 20.0, 30.0])) == [1.0, 2.0, 3.0]


def test_midrank_averages_tied_positions():
    # three-way tie occupying ranks 2,3,4 -> every member gets 3.0
    assert list(midrank([1.0, 5.0, 5.0, 5.0, 9.0])) == [1.0, 3.0, 3.0, 3.0, 5.0]


def test_midrank_is_order_independent():
    x = np.array([3.0, 1.0, 2.0, 2.0])
    order = np.argsort(x)
    assert sorted(midrank(x)) == sorted(midrank(x[order]))


# ---------------------------------------------------------------- delong
def _labels_scores():
    rng = np.random.default_rng(0)
    y = np.array([1] * 60 + [0] * 60)
    good = np.concatenate([rng.normal(1.0, 1.0, 60), rng.normal(-1.0, 1.0, 60)])
    return y, good


def test_delong_reports_no_difference_for_identical_scores():
    y, s = _labels_scores()
    r = delong(y, s, s)
    assert r["delta"] == 0.0
    assert r["auroc_a"] == r["auroc_b"]
    assert r["z"] is None and r["p_value"] is None      # variance of the difference is zero


def test_delong_auroc_agrees_with_sklearn():
    pytest.importorskip("sklearn")
    from sklearn.metrics import roc_auc_score
    y, s = _labels_scores()
    r = delong(y, s, s)
    assert r["auroc_a"] == pytest.approx(round(roc_auc_score(y, s), 4), abs=1e-4)


def test_delong_is_antisymmetric_in_its_arguments():
    y, good = _labels_scores()
    weak = good + np.random.default_rng(1).normal(0, 3.0, len(good))
    fwd, rev = delong(y, good, weak), delong(y, weak, good)
    assert fwd["delta"] == pytest.approx(-rev["delta"], abs=1e-9)
    assert fwd["p_value"] == pytest.approx(rev["p_value"], abs=1e-9)


def test_delong_detects_a_large_real_difference():
    y, good = _labels_scores()
    noise = np.random.default_rng(2).normal(0, 1.0, len(good))
    r = delong(y, good, noise)
    assert r["delta"] > 0.3
    assert r["p_value"] < 0.001


def test_delong_requires_both_classes():
    with pytest.raises(ValueError):
        delong([1, 1, 1], [0.1, 0.2, 0.3], [0.3, 0.2, 0.1])


# ---------------------------------------------------------------- bootstrap
def test_bootstrap_indices_shape_and_range():
    rng = np.random.default_rng(0)
    idx = bootstrap_indices(25, rng, n_boot=100)
    assert idx.shape == (100, 25)
    assert idx.min() >= 0 and idx.max() < 25


def test_bootstrap_indices_are_reproducible_from_a_seed():
    a = bootstrap_indices(10, np.random.default_rng(7), 20)
    b = bootstrap_indices(10, np.random.default_rng(7), 20)
    assert np.array_equal(a, b)


def test_percentile_ci_brackets_the_bulk_of_the_distribution():
    lo, hi = percentile_ci(list(range(101)))
    assert lo == pytest.approx(2.5, abs=0.5)
    assert hi == pytest.approx(97.5, abs=0.5)


def test_percentile_ci_ignores_nan_and_none():
    assert percentile_ci([1.0, 2.0, float("nan"), None, 3.0]) is not None
    assert percentile_ci([]) is None
    assert percentile_ci([None, float("nan")]) is None


# ---------------------------------------------------------------- mcnemar
def test_mcnemar_counts_only_discordant_pairs():
    a = [1, 1, 0, 0, 1, 1]
    b = [1, 0, 1, 0, 0, 0]
    m = mcnemar_exact(a, b)
    assert m["a_only_correct"] == 3          # positions 1, 4, 5
    assert m["b_only_correct"] == 1          # position 2
    assert m["discordant_pairs"] == 4


def test_mcnemar_is_one_when_nothing_is_discordant():
    m = mcnemar_exact([1, 0, 1], [1, 0, 1])
    assert m["discordant_pairs"] == 0
    assert m["p_value_exact"] == 1.0


def test_mcnemar_matches_the_exact_binomial():
    # 6 vs 2 discordant is the reported matched-evidence LIRICAL comparison
    m = mcnemar_exact([1] * 6 + [0] * 2, [0] * 6 + [1] * 2)
    from scipy.stats import binomtest
    assert m["p_value_exact"] == pytest.approx(round(binomtest(6, 8, 0.5).pvalue, 6))
    assert m["p_value_exact"] > 0.05          # the manuscript reports this as not significant


def test_mcnemar_reproduces_the_significant_full_evidence_arm():
    m = mcnemar_exact([1] * 9 + [0], [0] * 9 + [1])
    assert m["p_value_exact"] < 0.05
    assert round(m["p_value_exact"], 2) == 0.02


def test_mcnemar_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        mcnemar_exact([1, 0], [1, 0, 1])


# ---------------------------------------------------------------- rank metrics
def test_rank_metrics_on_hand_computed_ranks():
    m = rank_metrics([1, 2, 3, None], ks=(1, 3))
    assert m["n"] == 4
    assert m["recall@1"] == 0.25
    assert m["recall@3"] == 0.75
    assert m["mrr"] == pytest.approx(round((1 + 0.5 + 1 / 3) / 4, 4))


def test_rank_metrics_are_monotone_in_k():
    m = rank_metrics([1, 4, 9, None, 2], ks=(1, 3, 5, 10))
    assert m["recall@1"] <= m["recall@3"] <= m["recall@5"] <= m["recall@10"]


def test_rank_metrics_handle_no_hits_and_empty_input():
    m = rank_metrics([None, None])
    assert m["mrr"] == 0.0 and m["recall@1"] == 0.0
    assert rank_metrics([]) == {"n": 0}


def test_rank_metrics_perfect_ranking():
    m = rank_metrics([1, 1, 1])
    assert m["recall@1"] == 1.0 and m["mrr"] == 1.0
