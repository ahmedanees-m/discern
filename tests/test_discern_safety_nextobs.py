"""DISCERN Phases 4-6: abstention, management-aware safety, next observation, what-if."""
from core.dx_schemas import Feature, FeatureKind
from diseases.ontology import cluster_for
from jointdx.abstain import decide
from jointdx.factorgraph import Evidence, joint
from nextobs.partial import diagnose_partial
from nextobs.recommend import best_observation, recommend
from nextobs.whatif import whatif
from safety.interlock import flags


def _clin(fid, present):
    return Feature(fid, FeatureKind.CLINICAL, present, observed=present)


# ---- Phase 4: abstention ----
def test_decide_confident_case():
    cluster = cluster_for("integrin")
    ev = Evidence(variant_gene="ITGB3", genetic_codes=["PVS1", "PS1", "PM2"],
                  clinical=[_clin("glanzmann_type_bleeding", True), _clin("leukocytosis", False)])
    d = decide(cluster, ev, n_mc=60)
    assert d.leading == "gt"


def test_abstain_when_below_confidence_threshold():
    cluster = cluster_for("integrin")
    ev = Evidence(variant_gene="", genetic_codes=[], clinical=[])   # prior-only, low confidence
    d = decide(cluster, ev, tau=0.85, n_mc=60)
    assert not d.decided and "below threshold" in d.reason    # not confident enough -> abstain


# ---- Phase 5: management-aware safety ----
def test_safety_flag_fires_on_treatment_danger():
    cluster = cluster_for("integrin")
    # leaning GT, but LAD-III (HSCT) retains some probability (infections present)
    ev = Evidence(variant_gene="ITGB3", genetic_codes=["PVS1", "PM2"],
                  clinical=[_clin("glanzmann_type_bleeding", True), _clin("recurrent_infections", True)])
    fl = flags(cluster, joint(cluster, ev))
    assert any(f.competitor_id == "lad3" for f in fl)      # HSCT divergence flagged


def test_ddavp_interlock_hard_stop():
    cluster = cluster_for("vwf_gpib")
    # platelet-origin leaning (PT-VWD), but 2B retains probability; planned DDAVP -> hard stop
    ev = Evidence(variant_gene="GP1BA",
                  clinical=[Feature("ripa_low_dose_enhanced", FeatureKind.LAB, True, observed=True),
                            Feature("ripa_mixing_platelet_origin", FeatureKind.LAB, True, observed=True)])
    fl = flags(cluster, joint(cluster, ev), planned_tx="ddavp")
    assert any("HARD STOP" in f.message and f.competitor_id == "vwd2b" for f in fl)


def test_ddavp_hard_stop_fires_when_contraindicated_disease_is_leading():
    """v3.1 B3 regression: the hard-stop must fire even when the contraindicated disease is the
    LEADING call. Previously the loop `continue`d past the leading disease, so "DDAVP planned +
    2B leading" emitted only a divergence flag (or nothing) instead of a hard stop."""
    from jointdx.infer import leading_disease
    cluster = cluster_for("vwf_gpib")
    # plasma-origin RIPA + a VWF variant -> type 2B leads; DDAVP is contraindicated in 2B.
    ev = Evidence(variant_gene="VWF",
                  clinical=[Feature("ripa_low_dose_enhanced", FeatureKind.LAB, True, observed=True),
                            Feature("ripa_mixing_plasma_origin", FeatureKind.LAB, True, observed=True)])
    j = joint(cluster, ev)
    lead_id, _ = leading_disease(j)
    assert lead_id == "vwd2b"                                   # 2B is the leading call
    fl = flags(cluster, j, planned_tx="ddavp")
    hard = [f for f in fl if "HARD STOP" in f.message and f.competitor_id == "vwd2b"]
    assert hard, "hard-stop must fire on the contraindicated LEADING diagnosis"
    assert "[leading diagnosis]" in hard[0].message


# ---- Phase 6: next observation + what-if + partial ----
def test_next_observation_prefers_management_changing():
    cluster = cluster_for("integrin")
    ev = Evidence(variant_gene="ITGB3", genetic_codes=["PM2"],
                  clinical=[_clin("glanzmann_type_bleeding", True)])
    ranked = recommend(cluster, ev)
    assert ranked
    assert ranked[0][0].changes_management        # top recommendation changes management


def test_whatif_ripa_mixing_separates_diseases():
    cluster = cluster_for("vwf_gpib")
    ev = Evidence(variant_gene="VWF",
                  clinical=[Feature("ripa_low_dose_enhanced", FeatureKind.LAB, True, observed=True)])
    shifts = whatif(cluster, ev, "ripa_mixing")
    assert shifts["plasma_origin"][0] == "vwd2b"
    assert shifts["platelet_origin"][0] == "ptvwd"


def test_partial_input_runs_and_reports_missing():
    cluster = cluster_for("integrin")
    geno_only = Evidence(variant_gene="ITGB3", genetic_codes=["PVS1"])
    out = diagnose_partial(cluster, geno_only)
    assert "clinical" in out["missing"] and "functional" in out["missing"]
    assert out["leading"]
    assert best_observation(cluster, geno_only) is not None
