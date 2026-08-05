"""end-to-end diagnose() orchestration, intake, API."""
from api.main import handle_diagnose
from intake.extract import extract_features
from jointdx.factorgraph import Evidence
from jointdx.orchestrate import diagnose


def test_diagnose_end_to_end_gt():
    ev = Evidence(variant_gene="ITGB3", variant_id="ITGB3:c.X", genetic_codes=["PVS1", "PM2"],
                  clinical=[{"id": "glanzmann_type_bleeding"}])  # built below properly
    # build proper Features
    from core.dx_schemas import Feature, FeatureKind
    ev.clinical = [Feature("glanzmann_type_bleeding", FeatureKind.CLINICAL, True, observed=True),
                   Feature("recurrent_infections", FeatureKind.CLINICAL, True, observed=True)]
    rec = diagnose(ev, n_mc=40)
    assert rec is not None
    assert rec.posterior.leading in ("gt", "lad3")
    assert rec.next_observation is not None
    assert "Leading" in rec.explanation or "Undecidable" in rec.explanation
    # LAD-III (HSCT) competitor should be flagged given the infections
    assert any(f.competitor_id == "lad3" for f in rec.safety_flags)
    assert "ITGB3:c.X" in rec.reclassified_variants


def test_diagnose_unknown_gene_returns_none():
    assert diagnose(Evidence(variant_gene="NOTAGENE")) is None


def test_api_handle_diagnose_ddavp_interlock():
    out = handle_diagnose({
        "gene": "GP1BA", "variant_id": "GP1BA:c.Y", "planned_tx": "ddavp",
        "clinical": [{"id": "ripa_low_dose_enhanced"},
                     {"id": "ripa_mixing_platelet_origin"}],
    })
    assert out["leading"]
    assert any("HARD STOP" in f["message"] for f in out["safety_flags"])
    assert "explanation" in out and "audit" in out


def test_gp1ba_routes_to_ripa_cluster_when_bss_features_ruled_out():
    """Regression: GP1BA maps to both the macrothrombocytopenia and the enhanced-RIPA
    clusters. A platelet-type-VWD phenotype (enhanced RIPA present) with the BSS findings
    explicitly ruled out must still route to the enhanced-RIPA cluster and fire the DDAVP
    hard-stop; routing must key on present findings, not on which panels were filled in."""
    from core.dx_schemas import Feature, FeatureKind
    present = ["ripa_low_dose_enhanced", "ripa_mixing_platelet_origin", "normal_platelet_size"]
    absent = ["ripa_mixing_plasma_origin", "macrothrombocytopenia", "giant_platelets",
              "flow_cd42_reduced", "ripa_absent"]
    clin = ([Feature(f, FeatureKind.LAB, True, observed=True) for f in present]
            + [Feature(f, FeatureKind.LAB, False, observed=False) for f in absent])
    rec = diagnose(Evidence(variant_gene="GP1BA", clinical=clin), planned_tx="ddavp", n_mc=60)
    assert rec.posterior.cluster.id == "vwf_gpib"
    assert rec.posterior.leading == "ptvwd"
    assert any("HARD STOP" in f.message and "2B" in f.message for f in rec.safety_flags)

    # And the genuine BSS phenotype (macro findings present) must still route to BSS.
    bss = [Feature(f, FeatureKind.LAB, True, observed=True)
           for f in ("macrothrombocytopenia", "giant_platelets", "flow_cd42_reduced", "ripa_absent")]
    rec2 = diagnose(Evidence(variant_gene="GP1BA", clinical=bss), planned_tx="splenectomy", n_mc=60)
    assert rec2.posterior.cluster.id == "macrothrombocytopenia"
    assert rec2.posterior.leading == "bss"


def test_lad1_vs_lad3_activation_assay_discriminates():
    """The platelet aIIbb3 activation assay must separate LAD-I (activation intact, because
    ITGB2 is a leukocyte beta2 defect) from LAD-III (inside-out activation impaired). It was
    previously inert (declared in next_observations only, absent from the LR model)."""
    from core.dx_schemas import Feature, FeatureKind

    def post(gene, feats):
        clin = [Feature(f, FeatureKind.LAB, bool(p), observed=bool(p)) for f, p in feats.items()]
        return diagnose(Evidence(variant_gene=gene, clinical=clin), n_mc=80).posterior

    lad1 = post("ITGB2", {"leukocytosis": True, "recurrent_infections": True,
                          "aiib3_activation": True, "glanzmann_type_bleeding": False})
    assert lad1.leading == "lad1"
    lad3 = post("FERMT3", {"leukocytosis": True, "recurrent_infections": True,
                           "aiib3_activation": False, "glanzmann_type_bleeding": True})
    assert lad3.leading == "lad3"

    # the assay is a live input: normal vs reduced activation must move LAD-I's share.
    def share(gene, feats, did):
        return {d: v[0] for d, v in post(gene, feats).p_disease.items()}[did]
    normal = share("ITGB2", {"aiib3_activation": True}, "lad1")
    reduced = share("ITGB2", {"aiib3_activation": False}, "lad1")
    assert normal > reduced + 0.2


def test_intake_extracts_present_and_pertinent_negatives():
    note = "Mucocutaneous bleeding since birth. No leukocytosis. No recurrent infections."
    def stub(_note):
        return [("glanzmann_type_bleeding", True, "HP:0000001"),
                ("leukocytosis", False, "HP:0001974"),
                ("recurrent_infections", False, "HP:0002719")]
    feats = extract_features(note, extractor=stub)
    assert len(feats) == 3
    neg = [f for f in feats if not f.observed]
    assert {f.id for f in neg} == {"leukocytosis", "recurrent_infections"}   # pertinent negatives
