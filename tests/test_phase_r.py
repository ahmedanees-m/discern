"""Phase R guards: the gene term, the gene-blind safety view, and the committed headline numbers.

The metric assertions read committed JSON snapshots rather than recomputing, because the Phase R
harnesses need scikit-learn and CI installs only the light dependency set. The engine assertions
are recomputed live, since they are deterministic and cheap.
"""
from __future__ import annotations

import json
import os

import pytest

from core.dx_schemas import Feature, FeatureKind
from diseases.ontology import cluster_for
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import leading_disease
from jointdx.orchestrate import diagnose
from nextobs.whatif import whatif

HERE = os.path.dirname(__file__)
PR_VARIANT = os.path.join(HERE, "..", "bench", "phase_r_variant_metrics.json")
PR_SENS = os.path.join(HERE, "..", "bench", "phase_r_gene_term_sensitivity.json")
PR_GENE = os.path.join(HERE, "..", "eval", "gene_only_baseline.json")


def _clin(fid, present):
    return Feature(fid, FeatureKind.CLINICAL, present, observed=present)


# ---- the gene actually reaches the disease posterior (Phase R R4) ----
def test_gene_separates_haemophilia_a_from_b():
    """F8 and F9 give near-identical phenotypes; only the gene tells them apart."""
    cluster = cluster_for("coag_factor")
    shared = [_clin("prolonged_aptt", True), _clin("delayed_bleeding", True)]
    lead_a, _ = leading_disease(joint(cluster, Evidence(variant_gene="F8", clinical=shared)))
    lead_b, _ = leading_disease(joint(cluster, Evidence(variant_gene="F9", clinical=shared)))
    assert lead_a == "hemophilia_a"
    assert lead_b == "hemophilia_b"


def test_gene_term_is_evidence_not_a_filter():
    """An off-gene disease keeps real posterior mass, so it can still be argued back."""
    cluster = cluster_for("coag_factor")
    j = joint(cluster, Evidence(variant_gene="F8", clinical=[_clin("prolonged_aptt", True)]))
    from jointdx.infer import marginal_disease
    md = marginal_disease(j)
    assert md["hemophilia_b"] > 0.01           # not vetoed
    assert md["hemophilia_a"] > md["hemophilia_b"]


def test_deciding_assay_still_overturns_the_gene():
    """The value-of-information layer is only meaningful if a test can beat the gene."""
    cluster = cluster_for("vwf_gpib")
    ev = Evidence(variant_gene="VWF",
                  clinical=[Feature("ripa_low_dose_enhanced", FeatureKind.LAB, True, observed=True)])
    shifts = whatif(cluster, ev, "ripa_mixing")
    assert shifts["plasma_origin"][0] == "vwd2b"
    assert shifts["platelet_origin"][0] == "ptvwd"      # GP1BA disease, despite the VWF variant


# ---- safety is judged gene-blind (Phase R) ----
def test_offgene_contraindication_still_hard_stops():
    """A GP1BA case must still stop desmopressin: type 2B is a VWF disease but is not excluded."""
    ev = Evidence(variant_gene="GP1BA", clinical=[
        Feature("ripa_low_dose_enhanced", FeatureKind.LAB, True, observed=True),
        Feature("ripa_mixing_platelet_origin", FeatureKind.LAB, True, observed=True)])
    rec = diagnose(ev, planned_tx="ddavp", n_mc=40)
    assert rec is not None
    assert rec.posterior.leading == "ptvwd"
    assert any("HARD STOP" in f.message and "2B" in f.message for f in rec.safety_flags)


def test_safety_flags_are_reported_against_the_engines_own_leading_call():
    ev = Evidence(variant_gene="GP9", clinical=[
        Feature("flow_cd42_reduced", FeatureKind.LAB, True, observed=True),
        Feature("giant_platelets", FeatureKind.LAB, True, observed=True)])
    rec = diagnose(ev, planned_tx="splenectomy", n_mc=40)
    assert rec is not None
    assert all(f.leading_id == rec.posterior.leading for f in rec.safety_flags)


def test_safety_uses_the_gene_blind_view_not_the_gene_aware_one():
    """The re-validation the gene term made necessary, asserted rather than assumed.

    An ITGB3 case with recurrent infections must still flag LAD-III, whose gene is FERMT3. Under a
    gene-aware posterior LAD-III falls below the divergence threshold; under the gene-blind view the
    interlock still sees it. This pins the adjudication path, so the reported sensitivity and
    specificity cannot silently revert to the pre-Phase-R logic.
    """
    from jointdx.infer import marginal_disease
    ev = Evidence(variant_gene="ITGB3", genetic_codes=["PVS1", "PM2"],
                  clinical=[_clin("glanzmann_type_bleeding", True),
                            _clin("recurrent_infections", True)])
    cluster = cluster_for("integrin")
    p_aware = marginal_disease(joint(cluster, ev))["lad3"]
    p_blind = marginal_disease(joint(cluster, ev, gene_evidence=False))["lad3"]
    assert p_blind > p_aware                       # the gene would otherwise shrink the competitor
    rec = diagnose(ev, n_mc=40)
    assert rec is not None
    assert any(f.competitor_id == "lad3" for f in rec.safety_flags)


def test_safety_metrics_are_regenerated_under_current_logic():
    """Recompute rather than trust the committed snapshot."""
    from bench.track3_trustworthiness import safety_interlock
    s = safety_interlock()
    assert s["n_scenarios"] == 5
    assert s["hardstop_sensitivity"] == 1.0
    assert s["hardstop_specificity"] == 1.0
    snap = os.path.join(HERE, "..", "bench", "track3_metrics.json")
    if os.path.exists(snap):
        committed = json.load(open(snap, encoding="utf-8"))["safety_interlock"]
        assert committed["hardstop_sensitivity"] == s["hardstop_sensitivity"]
        assert committed["hardstop_specificity"] == s["hardstop_specificity"]
        assert committed["n_scenarios"] == s["n_scenarios"]


def test_variant_arm_is_independent_of_the_gene_term():
    """R1/R2/R3/R6 must be untouched by a disease-model change - proved, not assumed."""
    from jointdx import factorgraph
    from rules.variant_scoring import Annotations, score_variant

    probe = [("ITGB3", Annotations(af=1e-6, revel=0.9, consequence="missense")),
             ("F8", Annotations(af=None, revel=0.1, consequence="missense")),
             ("VWF", Annotations(af=1e-3, revel=0.5, consequence="nonsense"))]
    before = [(score_variant(g, g, a).points, score_variant(g, g, a).classification.name)
              for g, a in probe]
    on, off = factorgraph.ON_GENE, factorgraph.OFF_GENE
    try:
        factorgraph.ON_GENE, factorgraph.OFF_GENE = 0.99, 0.01
        after = [(score_variant(g, g, a).points, score_variant(g, g, a).classification.name)
                 for g, a in probe]
    finally:
        factorgraph.ON_GENE, factorgraph.OFF_GENE = on, off
    assert before == after


# ---- committed Phase R numbers ----
@pytest.mark.skipif(not os.path.exists(PR_VARIANT), reason="phase_r_variant_metrics.json not generated")
def test_calibration_is_out_of_sample_and_beaten_by_nobody():
    m = json.load(open(PR_VARIANT, encoding="utf-8"))
    for surf in ("eRepo_primary", "time_split"):
        s = m[surf]
        p = s["calibration_protocol"]
        assert p["reported_on"].startswith("held-out folds only")
        assert "training folds only" in p["fit_on"]
        cal = s["calibration"]
        # every ECE carries a CI, and post-hoc calibration is applied to the comparators too
        for entry in cal.values():
            assert entry["ece_ci95"] is not None
        assert "REVEL_isotonic_oof" in cal and "AlphaMissense_isotonic_oof" in cal
        # DISCERN stays best on the point estimate; the CIs overlap, which the paper must say
        assert cal["DISCERN_isotonic_oof"]["ece"] <= cal["REVEL_isotonic_oof"]["ece"]
        assert cal["DISCERN_isotonic_oof"]["ece"] <= cal["AlphaMissense_isotonic_oof"]["ece"]
        # post-hoc calibration helps the raw scores a lot - that is the retired claim
        assert cal["REVEL_isotonic_oof"]["ece"] < cal["REVEL_raw"]["ece"]


@pytest.mark.skipif(not os.path.exists(PR_VARIANT), reason="phase_r_variant_metrics.json not generated")
def test_discern_tracks_revel_rather_than_beating_it():
    m = json.load(open(PR_VARIANT, encoding="utf-8"))
    for surf in ("eRepo_primary", "time_split"):
        d = m[surf]["discern_vs_revel"]
        assert d["delong"]["p_value"] > 0.05          # the gap is not statistically significant
        lo, hi = d["paired_bootstrap"]["delta_ci95"]
        assert lo < 0 < hi                            # the CI on the difference spans zero


@pytest.mark.skipif(not os.path.exists(PR_VARIANT), reason="phase_r_variant_metrics.json not generated")
def test_intrinsic_evidence_cannot_reach_a_pathogenic_missense_band():
    """The most consequential finding, and the argument for the coupling."""
    m = json.load(open(PR_VARIANT, encoding="utf-8"))
    c = m["added_value_over_ranking_score"]["intrinsic_only_ceiling"]
    assert c["n_reaching_lp_or_p"] == 0
    assert c["max_discern_points_on_missense"] < c["lp_threshold_points"]


def test_pm1_and_pm5_are_variant_intrinsic_not_routed_away():
    """Guards the corrected claim: the partition does not move the hotspot and same-residue codes.

    The ceiling is not produced by routing PM1/PM5/PS1/PS4 elsewhere - they are variant-intrinsic
    and simply have no input in this pipeline. Saying otherwise was wrong in v1 of the manuscript
    and this test exists so it cannot be reintroduced.
    """
    from rules.vcep.partition import owner
    for code in ("PM1", "PM5", "PS1", "PS4"):
        assert owner(code) == "variant_intrinsic", code
    assert owner("PS3") == "functional"
    assert owner("PP4") == "disease_pp4"


@pytest.mark.skipif(not os.path.exists(PR_VARIANT), reason="phase_r_variant_metrics.json not generated")
def test_ceiling_attribution_separates_partition_from_missing_inputs():
    m = json.load(open(PR_VARIANT, encoding="utf-8"))["ceiling_attribution"]
    assert m["reach_lp_as_scored"] == 0
    # both streams contribute, and neither explains the ceiling on its own
    assert m["reach_lp_if_intrinsic_codes_were_available"] > 0
    assert m["reach_lp_if_routed_codes_were_re_added"] > 0
    assert m["reach_lp_with_both"] < m["n"] / 2       # still nowhere near resolving the set


PR_LIRICAL = os.path.join(HERE, "..", "eval", "lirical_arm.json")


@pytest.mark.skipif(not os.path.exists(PR_LIRICAL), reason="lirical_arm.json not generated")
def test_lirical_comparison_separates_reasoning_from_encoding():
    """The evidence-matched arm must exist, and its weaker result must not be quietly dropped."""
    m = json.load(open(PR_LIRICAL, encoding="utf-8"))
    assert m["headline_arm"] == "DISCERN_hpo_representable_only"
    matched = m["paired_tests_vs_lirical_restricted"]["hpo_representable_only"]
    full = m["paired_tests_vs_lirical_restricted"]["phenotype_only"]
    # on identical evidence DISCERN leads but not significantly - the claim must stay hedged
    assert matched["delta"] > 0
    assert matched["mcnemar"]["p_value_exact"] > 0.05
    assert matched["delta_ci95"][0] < 0 < matched["delta_ci95"][1]
    # the significant lead belongs to the arm that also encodes non-HPO findings
    assert full["mcnemar"]["p_value_exact"] < 0.05
    assert full["delta"] > matched["delta"]
    assert "two_distinct_claims" in m and "ordering_note" in m


@pytest.mark.skipif(not os.path.exists(PR_VARIANT), reason="phase_r_variant_metrics.json not generated")
def test_every_unavailable_intrinsic_code_has_a_stated_reason():
    m = json.load(open(PR_VARIANT, encoding="utf-8"))["ceiling_attribution"]
    reasons = m["why_each_intrinsic_code_has_no_input"]
    for code in m["intrinsic_but_no_input_here"]:
        assert code in reasons and len(reasons[code]) > 40, code


@pytest.mark.skipif(not os.path.exists(PR_GENE), reason="gene_only_baseline.json not generated")
def test_diagnosis_arm_is_reported_against_the_gene_only_baseline():
    m = json.load(open(PR_GENE, encoding="utf-8"))
    pooled = m["pooled"]
    # the baseline is strong, which is the finding: the pooled number cannot carry a claim
    assert pooled["gene_lookup"]["top1"] >= 0.9
    assert pooled["mcnemar_top1"]["p_value_exact"] > 0.05
    # and the subset where the gene does not settle it is tiny, which must be stated
    assert m["strata_sizes"].get("shared_gene", 0) <= 5


@pytest.mark.skipif(not os.path.exists(PR_SENS), reason="phase_r_gene_term_sensitivity.json not generated")
def test_committed_gene_term_is_the_largest_that_keeps_tests_decisive():
    m = json.load(open(PR_SENS, encoding="utf-8"))
    grid = sorted(m["grid"], key=lambda r: r["likelihood_ratio"])
    committed = next(r for r in grid if r["committed"])
    assert committed["deciding_assay_can_overturn_gene"] is True
    stronger = [r for r in grid if r["likelihood_ratio"] > committed["likelihood_ratio"]]
    assert stronger and all(not r["deciding_assay_can_overturn_gene"] for r in stronger)
    # safety never depends on the choice
    assert all(r["hardstop_sensitivity"] == 1.0 and r["hardstop_specificity"] == 1.0 for r in grid)
