"""DISCERN Phases 9-10: VUS-reclassification, misdx-rescue, reader study, VUS-triage."""
from core.dx_schemas import Feature, FeatureKind, VariantState
from eval.misdx_rescue import MisdxCase
from eval.misdx_rescue import evaluate as eval_rescue
from eval.reader_study import Vignette, score_discern
from eval.vus_reclass import ReclassCase
from eval.vus_reclass import evaluate as eval_reclass
from jointdx.factorgraph import Evidence
from triage.assay_priority import prioritize


def _clin(fid, present=True):
    return Feature(fid, FeatureKind.CLINICAL, present, observed=present)


def _func(fid, val="absent"):
    return Feature(fid, FeatureKind.FUNCTIONAL, val, observed=True)


def test_vus_reclassification_metric():
    cases = [
        ReclassCase(Evidence(variant_gene="ITGB3", genetic_codes=["PVS1", "PM2"],
                             clinical=[_clin("glanzmann_type_bleeding")],
                             functional=[_func("aiib3_expression_absent")]), VariantState.PATH),
        ReclassCase(Evidence(variant_gene="ITGB3", genetic_codes=["BA1"],
                             clinical=[_clin("glanzmann_type_bleeding")]), VariantState.BEN),
    ]
    rep = eval_reclass(cases)
    assert rep.n == 2 and 0.0 <= rep.reclassified_rate <= 1.0
    assert rep.concordance >= 0.5          # reclassified ones agree with the 3-star direction


def test_misdx_rescue_fires_for_true_disease():
    # looks like ITP (thrombocytopenia) but is BSS (giant platelets, reduced CD42) -> planned splenectomy
    ev = Evidence(variant_gene="GP9",
                  clinical=[_clin("macrothrombocytopenia"), _clin("giant_platelets")],
                  functional=[_func("flow_cd42_reduced", "reduced")])
    rep = eval_rescue([MisdxCase(ev, true_disease_id="bss", planned_tx="splenectomy")])
    assert rep.rescue_rate == 1.0


def test_reader_study_scoring():
    vig = Vignette(
        id="v1",
        ev=Evidence(variant_gene="ITGB3", genetic_codes=["PVS1", "PM2"],
                    clinical=[_clin("glanzmann_type_bleeding"), _clin("leukocytosis", False)]),
        correct_disease_id="gt", correct_next_obs_id="father_segregation")
    score = score_discern([vig])
    assert score.n == 1 and score.disease_accuracy == 1.0


def test_vus_triage_ranks_by_information_gain():
    cohort = [
        Evidence(variant_gene="ITGB3", variant_id="v_gt", genetic_codes=["PM2"],
                 clinical=[_clin("glanzmann_type_bleeding")]),
        Evidence(variant_gene="VWF", variant_id="v_vwf", genetic_codes=["PM2"],
                 clinical=[_clin("ripa_low_dose_enhanced")]),
    ]
    ranked = prioritize(cohort)
    assert ranked and ranked == sorted(ranked, key=lambda x: x.eig, reverse=True)
    assert all(r.best_assay for r in ranked)
