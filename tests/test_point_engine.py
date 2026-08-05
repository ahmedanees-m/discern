"""Deterministic classification band checks."""
from core.schemas import Classification, EvidenceContribution, PointsLedger, Strength, Variant
from rules.point_engine import classify


def _ledger(*contribs) -> PointsLedger:
    v = Variant("X", 100, "A", "T", gene="F8")
    return PointsLedger(variant=v, contributions=list(contribs))


def _c(code, strength, applied=True, reliability=1.0):
    return EvidenceContribution(code, strength, applied, "test", reliability=reliability)


def test_bands():
    assert classify(_ledger(_c("PVS1", Strength.PVS), _c("PM2", Strength.PM))) is Classification.P   # 10
    assert classify(_ledger(_c("PS3", Strength.PS), _c("PM2", Strength.PM))) is Classification.LP    # 6
    assert classify(_ledger(_c("PP3", Strength.PP))) is Classification.VUS                            # 1
    assert classify(_ledger(_c("BS1", Strength.BS))) is Classification.LB                             # -4
    assert classify(_ledger(_c("BVS1", Strength.BVS))) is Classification.B                            # -8


def test_unapplied_codes_excluded():
    led = _ledger(_c("PS3", Strength.PS), _c("PVS1", Strength.PVS, applied=False))
    assert led.points == 4
    assert classify(led) is Classification.VUS


def test_reliability_downweight():
    # An equity down-weighted PM2 contributes less than its nominal points.
    led = _ledger(_c("PM2", Strength.PM, reliability=0.5))
    assert led.points == 1.0


def test_ba1_prefilter():
    class Spec:
        version = "test-1"

        def ba1_met(self, ledger):
            return True

    led = _ledger(_c("PVS1", Strength.PVS), _c("PS3", Strength.PS))  # would be Pathogenic
    assert classify(led, Spec()) is Classification.B                  # BA1 overrides
