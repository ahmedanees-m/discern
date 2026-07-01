"""gnomAD PM2/BS1/BA1 adapter + VCEP loader (plan Step 1.1 / 0.2)."""
from adapters.gnomad import GnomadAdapter
from core.adapter import EvidenceAdapter
from core.schemas import EvidenceContribution, PatientContext, PointsLedger, Strength, Variant
from rules.vcep_loader import get_spec, load_spec  # noqa: F401


def _adapter(af):
    return GnomadAdapter(get_spec("F8"), af_lookup=lambda _vid: af)


def _v():
    return Variant("X", 154861759, "C", "T", gene="F8")


def test_is_adapter():
    a = _adapter(0.0)
    assert isinstance(a, EvidenceAdapter)
    assert a.code_group == ("PM2", "BS1", "BA1")
    assert a.health_check()


def test_ba1_high_frequency():
    out = _adapter(0.01).evaluate(_v(), PatientContext())
    assert out[0].code == "BA1" and out[0].strength is Strength.BVS


def test_bs1_moderate_frequency():
    out = _adapter(0.0005).evaluate(_v(), PatientContext())   # > bs1 (1e-4), < ba1 (1.5e-3)
    assert out[0].code == "BS1" and out[0].strength is Strength.BS


def test_pm2_absent():
    assert _adapter(None).evaluate(_v(), PatientContext())[0].code == "PM2"
    assert _adapter(0.0).evaluate(_v(), PatientContext())[0].code == "PM2"


def test_intermediate_no_call():
    # between pm2 (1e-5) and bs1 (1e-4) -> no frequency code
    assert _adapter(5e-5).evaluate(_v(), PatientContext()) == []


# ---- VCEP loader ----
def test_spec_f8_loaded():
    s = get_spec("F8")
    assert s.gene == "F8" and s.inheritance == "XLR" and s.mechanism == "coag_factor"
    assert s.ba1_threshold() > s.bs1_threshold() > s.pm2_threshold()


def test_spec_fallback_to_base():
    s = get_spec("NOSUCHGENE")
    assert s.gene == "*"        # base_acmg fallback


def test_ba1_met_prefilter():
    led = PointsLedger(_v(), contributions=[
        EvidenceContribution("BA1", Strength.BVS, True, "gnomAD")])
    assert get_spec("F8").ba1_met(led)
