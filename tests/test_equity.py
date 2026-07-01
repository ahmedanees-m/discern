"""Phase 3 equity: reliability down-weighting, routing, ancestry, dashboard."""
from core.schemas import (
    Action,
    EvidenceContribution,
    PatientContext,
    PointsLedger,
    Strength,
    Variant,
)
from engine.recommend import recommend_next_action
from equity.ancestry import assign_group
from equity.dashboard import case_equity, cohort_resolution_gap
from equity.reliability import apply_reliability, reliability_pm2
from equity.routing import equity_filter, is_equitable


def _led(codes):
    v = Variant("17", 1, "C", "T", gene="ITGB3")
    return PointsLedger(v, contributions=[EvidenceContribution(c, s, True, "t") for c, s in codes])


# ---- reliability ----
def test_reliability_monotonic_in_group_size():
    assert reliability_pm2(0) == 0.3                 # no reference data -> floor
    assert reliability_pm2(10) < reliability_pm2(100000)
    assert reliability_pm2(10**6) == 1.0             # ample data -> full reliability


def test_apply_reliability_downweights_biased_codes():
    led = _led([("PM2", Strength.PM), ("PS3", Strength.PS)])
    assert led.points == 6.0
    apply_reliability(led, group_n=0)                # sparse group
    # PM2 down-weighted to 0.3 -> 2*0.3 + 4*1.0 = 4.6; PS3 (functional) unchanged
    assert abs(led.points - 4.6) < 1e-9


# ---- routing ----
def _act(name, mod, yields):
    return Action(name, yields, 500, 7, 0.9, 0.9, mod)


def test_is_equitable():
    assert is_equitable(_act("flow", "functional", ["PS3"]))
    assert not is_equitable(_act("freq_only", "frequency", ["PM2"]))


def test_equity_filter_routes_to_robust_for_underrepresented():
    patient = PatientContext(ancestry_group="sas", ancestry_confidence=0.4)
    actions = [_act("freq_only", "frequency", ["PM2"]), _act("flow", "functional", ["PS3"])]
    routed = equity_filter(actions, patient)
    assert [a.name for a in routed] == ["flow"]       # biased-only dropped
    # well-characterized ancestry -> no filtering
    assert len(equity_filter(actions, PatientContext(ancestry_group="nfe",
                                                     ancestry_confidence=0.95))) == 2


# ---- ancestry ----
def test_assign_group_nearest_centroid():
    centroids = {"nfe": [0.0, 0.0], "sas": [1.0, 1.0]}
    r = assign_group([0.9, 0.95], centroids)
    assert r.group == "sas" and 0.0 < r.confidence <= 1.0


# ---- dashboard ----
def test_case_equity_fraction():
    led = _led([("PM2", Strength.PM), ("PS3", Strength.PS)])   # 2 biased / 4 equitable
    ce = case_equity(led)
    assert abs(ce.biased_fraction - (2.0 / 6.0)) < 1e-9


def test_cohort_gap():
    out = cohort_resolution_gap({"nfe": (80, 100), "sas": (50, 100)})
    assert abs(out["resolution_gap"] - 0.30) < 1e-9


# ---- integration: equity routing in recommend ----
def test_recommend_equity_routing_runs():
    led = _led([("PP3", Strength.PP)])
    rec = recommend_next_action(led, PatientContext(ancestry_group="sas", ancestry_confidence=0.3),
                                mechanism="integrin_expression", equity_routing=True)
    assert all(r["modality"] in {"functional", "rna", "segregation"} for r in rec.ranked)
