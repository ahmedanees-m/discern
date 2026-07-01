"""v3.1 Track C1: the synthetic-coupling sanity harness runs and its circularity guard holds."""
from eval.synthetic_coupling_harness import run


def test_harness_runs_and_guard_holds():
    r = run("integrin", n=120, seed=1)
    assert r.n == 120
    # circularity guard: genetics-only variant marginal is invariant to phenotype (no double-count)
    assert r.intrinsic_invariant_to_phenotype is True
    # both arms produce a valid probability mass
    assert 0.0 <= r.intrinsic_only_path_mass <= 1.0
    assert 0.0 <= r.coupled_path_mass <= 1.0


def test_coupling_channel_is_active():
    # on synthetic disease-specific phenotype, the coupled arm should move P(PATH) vs genetics-only
    # (a wiring sanity check, NOT evidence the coupling is correct - that is cohort-gated, G13)
    r = run("integrin", n=200, seed=2)
    assert r.coupled_path_mass != r.intrinsic_only_path_mass
