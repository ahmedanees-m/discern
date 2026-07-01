"""Validation modules: ablation, calibration, equity-eval, retrospective (Phase 6)."""
from eval.ablation import PairedCase, run_ablation
from eval.calibration import expected_calibration_error, reliability_curve
from eval.equity_eval import compare_equity_configs
from eval.retrospective import Journey, replay


def test_ablation_detects_full_system_advantage():
    # full system resolves 20 extra cases that tools-only misses, at lower cost
    cases = [PairedCase(True, False, 500, 800) for _ in range(20)]
    cases += [PairedCase(True, True, 500, 800) for _ in range(80)]
    rep = run_ablation(cases)
    assert rep.full_resolution_rate == 1.0 and rep.base_resolution_rate == 0.8
    assert rep.mcnemar_p < 0.05          # significant resolution advantage
    assert rep.mean_cost_delta < 0       # full system cheaper on average


def test_calibration_perfect_is_zero_ece():
    pred = [0.05, 0.25, 0.55, 0.85, 0.95]
    # outcomes matching the bin means -> low ECE
    assert expected_calibration_error([0.0, 0.0, 1.0, 1.0, 1.0], [0, 0, 1, 1, 1]) == 0.0
    curve = reliability_curve(pred, [0, 0, 1, 1, 1], n_bins=5)
    assert len(curve) == 5


def test_equity_routing_narrows_gap():
    off = {"nfe": (90, 100), "sas": (60, 100)}      # 0.30 gap
    on = {"nfe": (88, 100), "sas": (78, 100)}       # 0.10 gap
    cmp = compare_equity_configs(off, on)
    assert cmp.narrowed and cmp.gap_equity_on < cmp.gap_equity_off


def test_retrospective_concordance():
    js = [Journey("c1", "factor_activity", "factor_activity", 1, 250, 900),
          Journey("c2", "rna_seq", "flow_cytometry", 2, 450, 600)]
    rep = replay(js)
    assert rep.top1_concordance == 0.5 and rep.top3_concordance == 1.0
    assert rep.cheaper_fraction == 1.0
