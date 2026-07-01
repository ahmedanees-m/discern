"""v3.1 Track A2: novel-variant intrinsic scoring + PVS1/PS4 strength trees + spec-aware adapters."""
from adapters.insilico import InSilicoAdapter
from adapters.splice import SpliceAdapter
from core.schemas import Strength
from rules.point_engine import Classification
from rules.variant_scoring import Annotations, ps4_strength, pvs1_strength, score_variant
from rules.vcep.loader import get_spec


def test_high_frequency_variant_is_benign():
    sv = score_variant("F8", "x", Annotations(af=0.01))   # F8 BA1 = 0.000333
    assert "BA1" in sv.codes
    assert sv.classification == Classification.B


def test_rare_lof_with_nmd_is_pathogenic_leaning():
    sv = score_variant("F8", "x", Annotations(af=0.0, consequence="frameshift", nmd_predicted=True))
    assert any(c.startswith("PVS1_VeryStrong") for c in sv.codes)
    assert "PM2_Supporting" in sv.codes                    # F8 applies PM2 at Supporting
    assert sv.classification in (Classification.P, Classification.LP)


def test_missing_predictors_gives_reduced_confidence():
    sv = score_variant("F8", "x", Annotations(af=0.0))     # no REVEL/splice provided
    assert sv.confidence.startswith("reduced")


def test_gene_without_spec_is_reduced_confidence():
    sv = score_variant("SERPINF2", "x", Annotations(af=0.0))   # no VCEP spec -> Gate G2
    assert sv.covered is False
    assert sv.confidence.startswith("reduced")


def test_computational_pp3_and_bp4():
    hi = score_variant("F8", "x", Annotations(af=0.0, revel=0.9, splice=0.0))
    assert "PP3_Supporting" in hi.codes
    lo = score_variant("VWF", "x", Annotations(af=0.005, revel=0.1, splice=0.0))  # VWF BP4 REVEL<=0.290
    assert "BP4_Supporting" in lo.codes


def test_pvs1_tree_levels():
    assert pvs1_strength(Annotations(consequence="frameshift", nmd_predicted=True)) == Strength.PVS
    assert pvs1_strength(Annotations(consequence="frameshift", nmd_predicted=False,
                                     in_functional_domain=True)) == Strength.PS
    assert pvs1_strength(Annotations(consequence="frameshift", nmd_predicted=False)) == Strength.PM
    assert pvs1_strength(Annotations(consequence="missense")) is None


def test_ps4_proband_ratio_tree():
    assert ps4_strength(Annotations(proband_count=20)) == Strength.PS
    assert ps4_strength(Annotations(proband_count=6)) == Strength.PM
    assert ps4_strength(Annotations(proband_count=2)) == Strength.PP
    assert ps4_strength(Annotations(proband_count=None)) is None


# ---- v3.1 refinements: full PVS1 tree, OR-based PS4, RUNX1 thresholds, spec computational + adapters ----

def test_pvs1_full_tree_branches():
    # canonical splice undergoing NMD -> VeryStrong
    assert pvs1_strength(Annotations(consequence="canonical_splice", nmd_predicted=True)) == Strength.PVS
    # initiation-codon variant -> Moderate
    assert pvs1_strength(Annotations(consequence="start_lost")) == Strength.PM
    # NMD-escaping truncation that removes >10% of the protein -> Strong
    assert pvs1_strength(Annotations(consequence="frameshift", nmd_predicted=False,
                                     removes_gt10pct=True)) == Strength.PS
    # in-frame deletion in a critical domain -> Strong
    assert pvs1_strength(Annotations(consequence="inframe_deletion", in_functional_domain=True)) == Strength.PS
    # not in a biologically-relevant transcript -> downgraded to Supporting
    assert pvs1_strength(Annotations(consequence="frameshift", nmd_predicted=True,
                                     exon_in_biorelevant_transcript=False)) == Strength.PP


def test_ps4_odds_ratio_tree():
    assert ps4_strength(Annotations(odds_ratio=6.0, or_ci_low=2.0)) == Strength.PS
    assert ps4_strength(Annotations(odds_ratio=3.5, or_ci_low=1.5)) == Strength.PM
    assert ps4_strength(Annotations(odds_ratio=2.2, or_ci_low=1.1)) == Strength.PP
    assert ps4_strength(Annotations(odds_ratio=10.0, or_ci_low=0.8)) is None   # CI includes 1


def test_runx1_thresholds_verified():
    s = get_spec("RUNX1")
    assert s.af["ba1"] == 0.0015 and s.af["bs1"] == 0.00015     # MM-VCEP (Luo 2019)
    assert s.computational["pp3_revel"] == 0.88                  # MM-VCEP v2 revision (2022)
    assert s.computational["bp4_revel"] == 0.50


def test_computational_cutoffs_are_gene_specific():
    # REVEL 0.62: above F8's PP3 cut-off (0.6) but below VWF's (0.644)
    f8 = score_variant("F8", "x", Annotations(af=0.0, revel=0.62))
    vwf = score_variant("VWF", "x", Annotations(af=0.005, revel=0.62))
    assert "PP3_Supporting" in f8.codes
    assert not any(c.startswith("PP3") for c in vwf.codes)


def test_gp1ba_graded_pp3_moderate():
    sv = score_variant("GP1BA", "x", Annotations(af=0.00005, revel=0.8))   # >= 0.773 moderate cut-off
    assert "PP3_Moderate" in sv.codes


def test_insilico_adapter_for_spec_uses_vcep_cutoffs():
    assert InSilicoAdapter.for_spec(get_spec("F8"))._spec_strength(0.65) == ("PP3", Strength.PP)
    assert InSilicoAdapter.for_spec(get_spec("GP1BA"))._spec_strength(0.8) == ("PP3", Strength.PM)
    assert InSilicoAdapter.for_spec(get_spec("VWF"))._spec_strength(0.2) == ("BP4", Strength.BP)
    # no spec computational -> falls back to Pejaver-calibrated bands (gene-agnostic)
    assert InSilicoAdapter(predictor="REVEL", score_lookup=lambda v: 0.5)._spec_aware is False


def test_splice_adapter_for_spec_uses_vcep_cutoff():
    assert SpliceAdapter.for_spec(get_spec("F8"))._pp3 == 0.2     # CSpec GN071 SpliceAI cut-off
    assert SpliceAdapter.for_spec(get_spec("VWF"))._pp3 == 0.5
