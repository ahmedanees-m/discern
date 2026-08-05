"""Tavtigian posterior bridge checks."""
from rules.posterior import posterior


def test_anchor_points():
    # prior 0.10; the four canonical anchors from the plan.
    assert abs(posterior(10) - 0.99) < 0.02      # Pathogenic
    assert abs(posterior(6) - 0.90) < 0.03       # Likely pathogenic
    assert abs(posterior(0) - 0.10) < 1e-6        # prior
    assert posterior(-7) < 0.001                  # Benign


def test_monotonic():
    pts = [-8, -4, 0, 2, 4, 6, 8, 10]
    vals = [posterior(p) for p in pts]
    assert all(b > a for a, b in zip(vals, vals[1:], strict=False))


def test_prior_shifts_curve():
    assert posterior(0, prior_p=0.25) == 0.25
    assert posterior(4, prior_p=0.25) > posterior(4, prior_p=0.10)
