"""Regression guards for the benchmark tracks (DISCERN_Benchmark_Execution_Plan_v1).

Track 1 assertions read the committed metrics snapshot (bench/track1_metrics.json, derived from the
committed GeneBe cache) so they are network-free. Track 3 safety is recomputed live (deterministic:
the hard-stop fires iff the planned tx is a contraindication of a non-excluded disease).
"""
from __future__ import annotations

import json
import os

import pytest

HERE = os.path.dirname(__file__)
T1 = os.path.join(HERE, "..", "bench", "track1_metrics.json")


@pytest.mark.skipif(not os.path.exists(T1), reason="track1_metrics.json not generated")
def test_track1_discrimination_and_calibration():
    m = json.load(open(T1, encoding="utf-8"))
    disc = m["discrimination"]
    # DISCERN tracks REVEL (it ingests REVEL as PP3) and clears the legacy InterVar anchor.
    assert disc["DISCERN_points"]["auroc"] >= 0.90
    assert disc["DISCERN_points"]["auroc"] > disc["InterVar_full_DB_prior"]["auroc"]
    # Calibration is the differentiator: DISCERN's calibrated prob beats the raw predictor scores,
    # and the categorical tools are not calibratable at all.
    cal = m["calibration"]
    assert cal["DISCERN_isotonic_oof"]["ece"] < cal["REVEL_raw_score"]["ece"]
    assert cal["DISCERN_isotonic_oof"]["ece"] < cal["AlphaMissense_raw_score"]["ece"]
    assert "not calibratable" in cal["GeneBe"]["note"]
    # ClinVar circularity: GeneBe reproduces ClinVar -> not a fair comparison surface.
    cv = m["clinvar_circularity"]
    assert cv["genebe_auroc_overall"] >= 0.99


def test_track3_safety_interlock_both_directions():
    from bench.track3_trustworthiness import safety_interlock
    s = safety_interlock()
    # fires on every contraindicated plan AND stays silent on every harmless plan
    assert s["hardstop_sensitivity"] == 1.0
    assert s["hardstop_specificity"] == 1.0


def test_track3_no_confidently_wrong_and_monotone_payoff():
    from bench.track3_trustworthiness import _case_predictions, calibration_diagnosis, risk_coverage
    preds = _case_predictions()
    cal = calibration_diagnosis(preds)
    assert cal["confidently_wrong_at_0.8"] == 0          # no high-confidence diagnosis errors
    rc = risk_coverage(preds)
    assert rc["monotone_payoff"] is True                 # most-confident half is at least as accurate


T1B = os.path.join(HERE, "..", "bench", "track1b_erepo_metrics.json")


@pytest.mark.skipif(not os.path.exists(T1B), reason="track1b_erepo_metrics.json not generated")
def test_track1b_calibration_edge_holds_on_independent_surface():
    m = json.load(open(T1B, encoding="utf-8"))
    # G-T1': the calibration edge must hold on BOTH the eRepo-primary and the time-split surface.
    for surf in ("eRepo_primary", "time_split"):
        cal = m[surf]["calibration"]
        assert cal["DISCERN_isotonic_oof"]["ece"] < cal["REVEL_raw"]["ece"]
        assert cal["DISCERN_isotonic_oof"]["ece"] < cal["AlphaMissense_raw"]["ece"]
        assert m[surf]["discrimination"]["DISCERN"]["auroc"] >= 0.90
    # partition made visible: DISCERN agrees on the in-silico code it derives (PP3) and applies NONE
    # of the evidence-stream codes it is not given inputs for (PS3 functional, PS4 case-control).
    k = m["per_code_kappa_vs_erepo"]
    assert k["PP3"]["kappa"] >= 0.7
    assert k["PS3"]["discern_applied"] == 0
    assert k["PS4"]["discern_applied"] == 0
    # Strength-modified codes must be counted. PM2 is recorded as PM2_Supporting by most panels, so
    # a parser that drops the suffix collapses its expert count to a couple of dozen; the true
    # count is in the hundreds. This is the number that exposed the parser defect.
    assert k["PM2"]["erepo_applied"] > 500, "strength-modified codes are being dropped again"
    assert k["PVS1"]["erepo_applied"] > 150
    # Agreement is reported at two levels, and both must be present for every applied criterion.
    for v in k.values():
        if v["discern_applied"] and v["erepo_applied"]:
            assert v["both_applied"] is not None
            assert v["strength_agreement"] is not None


def test_criterion_parsing_matches_the_partition_vocabulary():
    """The kappa path and the coverage path must consume one normalized code vocabulary.

    An earlier parser anchored each criterion with \b. Because "_" is a word character, "\bPM2\b"
    does not match "PM2_Supporting", so every strength-modified code was silently dropped from the
    agreement analysis while the coverage analysis, which normalizes with base_code, kept them.
    Two paths over one input, one normalized and one not, is invisible until someone reads the
    supplementary table. This asserts they agree.
    """
    from bench.track1b_erepo_headtohead import _erepo_codes
    from rules.vcep.partition import base_code

    applied = ("PM2_Supporting, PP3_Moderate, PVS1_Strong, PM3_Very Strong, "
               "PS3_Supporting, BS1_Supporting, PP1_Strong, PM2, PP3")
    parsed = _erepo_codes(applied)

    # Every criterion present in the raw string is recovered, in base form.
    assert {"PM2", "PP3", "PVS1", "PM3", "PS3", "BS1", "PP1"} == parsed
    # No strength suffix survives into the criterion vocabulary.
    assert not any("_" in c for c in parsed)
    # The two paths agree token for token.
    raw = [t.strip() for t in applied.replace(";", ",").split(",") if t.strip()]
    assert {base_code(t) for t in raw} == parsed


def test_strength_modified_codes_are_not_dropped():
    """A strength-modified code must be counted, not discarded."""
    from bench.track1b_erepo_headtohead import _erepo_codes, _erepo_codes_with_strength

    assert _erepo_codes("PM2_Supporting") == {"PM2"}
    assert _erepo_codes_with_strength("PM2_Supporting") == {"PM2": "Supporting"}
    # An unmodified code carries its criterion's default strength, so these two agree.
    assert _erepo_codes_with_strength("PP3")["PP3"] == "Supporting"
    assert _erepo_codes_with_strength("PP3_Supporting")["PP3"] == "Supporting"
    assert _erepo_codes_with_strength("PVS1")["PVS1"] == "Very Strong"
