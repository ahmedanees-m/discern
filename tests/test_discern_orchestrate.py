"""DISCERN Phase 8: end-to-end diagnose() orchestration, intake, API."""
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
