"""Remaining Phase 1 adapters: insilico, splice, autopvs1, mave, litmine, prioritizer."""
from adapters.autopvs1 import AutoPVS1Adapter
from adapters.insilico import InSilicoAdapter, revel_strength
from adapters.litmine_pm3 import LitMineAdapter, MinedEvidence
from adapters.mave import MaveAdapter, oddspath_to_strength
from adapters.prioritizer import normalize
from adapters.splice import SpliceAdapter, SpliceResult
from core.adapter import EvidenceAdapter
from core.schemas import PatientContext, Strength, Variant


def _v(gene="BRCA1"):
    return Variant("17", 43000000, "C", "T", gene=gene, hgvs_p="p.Arg100Gln")


def test_insilico_revel_calibration():
    # Pejaver 2022 verified thresholds
    assert revel_strength(0.95) == ("PP3", Strength.PS)
    assert revel_strength(0.80) == ("PP3", Strength.PM)
    assert revel_strength(0.70) == ("PP3", Strength.PP)
    assert revel_strength(0.20) == ("BP4", Strength.BP)
    assert revel_strength(0.10) == ("BP4", Strength.BM)
    assert revel_strength(0.50) is None          # intermediate -> no call
    a = InSilicoAdapter(score_lookup=lambda v: 0.95)
    assert isinstance(a, EvidenceAdapter) and a.health_check()
    out = a.evaluate(_v(), PatientContext())
    assert out[0].code == "PP3" and out[0].strength is Strength.PS


def test_splice_pp3_and_rna_flag():
    a = SpliceAdapter(splice_lookup=lambda v: SpliceResult(0.8, True))
    assert a.evaluate(_v(), PatientContext())[0].code == "PP3"
    assert a.rna_resolvable(_v())                 # high score at splice site -> RNA-resolvable
    b = SpliceAdapter(splice_lookup=lambda v: SpliceResult(0.05, True))
    assert b.evaluate(_v(), PatientContext())[0].code == "BP4"


def test_mave_oddspath():
    assert oddspath_to_strength(20) == ("PS3", Strength.PS)
    assert oddspath_to_strength(5) == ("PS3", Strength.PM)
    assert oddspath_to_strength(0.04) == ("BS3", Strength.BS)
    a = MaveAdapter(oddspath_lookup=lambda v: 25.0)
    assert a.evaluate(_v(), PatientContext())[0].code == "PS3"


def test_autopvs1():
    a = AutoPVS1Adapter(pvs1_lookup=lambda v: "Very Strong")
    out = a.evaluate(_v(), PatientContext())
    assert out[0].code == "PVS1" and out[0].strength is Strength.PVS


def test_litmine_requires_citations():
    grounded = LitMineAdapter(miner=lambda v, p: [
        MinedEvidence("PM3", Strength.PM, ["12345678"], "in trans with pathogenic")])
    assert grounded.evaluate(_v(), PatientContext())[0].code == "PM3"
    ungrounded = LitMineAdapter(miner=lambda v, p: [
        MinedEvidence("PM3", Strength.PM, [], "no citation")])
    assert ungrounded.evaluate(_v(), PatientContext()) == []     # ungrounded -> rejected


def test_prioritizer_disagreement_flag():
    cs = normalize({"exomiser": [("F8", 0.9)], "ai_marrvel": [("VWF", 0.8)]})
    assert cs.disagreement_flag                  # different top genes -> human review
    cs2 = normalize({"exomiser": [("F8", 0.9)], "ai_marrvel": [("F8", 0.85)]})
    assert not cs2.disagreement_flag
