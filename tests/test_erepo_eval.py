"""CI guard for the ERepo per-code eval (plan validation cleanup).

Runs eval/erepo_reconstruction.run on a tiny hand-built ERepo fixture (including
bundled-code variants). Asserts the partition vocabulary covers every applied code
(0 unknown) so a future VCEP code we have not mapped fails loudly in CI rather than
landing silently in `unknown`, and that the at-risk (owned-code) population is detected.
"""
import os

from eval.erepo_reconstruction import run

_FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "erepo_sample.tab")


def test_partition_covers_all_codes():
    s = run(_FIXTURE)
    assert s["partition_coverage"] == 1.0, f"unmapped codes: {s['unknown_codes']}"


def test_at_risk_population_detected():
    s = run(_FIXTURE)
    assert s["bleeding_variants"] >= 5
    assert s["variants_with_owned_codes"] > 0          # bundled PP4/PS3/PM3/PP1 variants present
    assert s["points_routed_out"].get("disease_pp4", 0) > 0   # PP4 routed to the coupling


def test_combining_arithmetic_runs():
    s = run(_FIXTURE)
    assert 0.0 <= s["acmg_combining_exact"] <= 1.0
    assert s["inflation_prevented"] >= 0
