"""Odyssey simulator harness (plan Phase 2.8 / 6.2).

These assertions validate the harness and the *defensible* findings only - VOI beats
the naive (random) baseline, graceful-degradation saves cost, and VOI uses no more
tests than greedy. We deliberately do NOT assert VOI dollar-dominance over a strong
cost-aware greedy (it does not hold on a homogeneous cohort - reported as such).
"""
from sim.odyssey_sim import compare


def test_harness_runs_all_policies():
    stats = compare(n_cases=800, seed=7)
    assert set(stats) == {"voi", "voi_stop", "greedy", "fixed", "random"}
    for s in stats.values():
        assert 0.0 <= s.resolved_rate <= 1.0
        assert s.mean_cost > 0 and s.mean_tests >= 1.0


def test_voi_beats_random_baseline():
    stats = compare(n_cases=800, seed=7)
    # VOI is more accurate AND cheaper than picking actions at random
    assert stats["voi"].accuracy > stats["random"].accuracy
    assert stats["voi"].mean_cost < stats["random"].mean_cost


def test_graceful_degradation_saves_cost():
    stats = compare(n_cases=800, seed=7)
    # stopping when no action improves resolution avoids futile spend
    assert stats["voi_stop"].mean_cost < stats["greedy"].mean_cost


def test_voi_no_more_tests_than_greedy():
    stats = compare(n_cases=800, seed=7)
    assert stats["voi"].mean_tests <= stats["greedy"].mean_tests + 0.05
