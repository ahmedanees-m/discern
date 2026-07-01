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
