"""Domain-agnosticism (Phase 7 packs) + API handlers (Section 9)."""
from api.main import handle_case, handle_classify, handle_recommend
from core.schemas import EvidenceContribution, PatientContext, PointsLedger, Strength, Variant
from engine.recommend import recommend_next_action
from rules.vcep_loader import get_spec


def _led(gene, codes):
    v = Variant("2", 1, "C", "T", gene=gene)
    return PointsLedger(v, contributions=[EvidenceContribution(c, s, True, "t") for c, s in codes])


# ---- domain-agnosticism: same engine, new packs only ----
def test_epilepsy_pack_routes_to_patch_clamp():
    led = _led("SCN1A", [("PM2", Strength.PM), ("PP3", Strength.PP)])
    rec = recommend_next_action(led, PatientContext(), mechanism="ion_channel_loss",
                                domain="epilepsy")
    assert rec.ranked
    assert any("patch_clamp" in r["action"] for r in rec.ranked)     # electrophysiology assay
    # top action should be high-yield (PS3-class functional or RNA), not weak segregation
    assert rec.ranked[0]["modality"] in {"functional", "rna"}


def test_cancer_pack_routes_to_mave():
    led = _led("BRCA1", [("PM2", Strength.PM), ("PP3", Strength.PP)])
    rec = recommend_next_action(led, PatientContext(), mechanism="brca", domain="cancer")
    assert rec.ranked
    assert any(r["action"] in ("sge_saturation", "mave_functional", "hdr_assay")
               for r in rec.ranked)


def test_specs_loaded_for_new_domains():
    assert get_spec("SCN1A").mechanism == "ion_channel_loss"
    assert get_spec("BRCA1").mechanism == "brca"


# ---- API handlers ----
def test_api_classify():
    out = handle_classify({"gene": "F8", "chrom": "X", "pos": 1, "ref": "C", "alt": "T",
                           "codes": ["PVS1", "PM2"]})              # 8 + 2 = 10 -> Pathogenic
    assert out["classification"] == "P"
    assert out["points"] == 10.0


def test_api_recommend():
    out = handle_recommend({"gene": "ITGB3", "codes": ["PM2", "PP3"],
                            "mechanism": "integrin_expression", "domain": "bleeding"})
    assert out["current_class"] == "VUS"
    assert out["actions"] and "explanation" in out


def test_api_case():
    out = handle_case({"has_vus_candidate": True, "has_any_candidate": True})
    assert out["route"] == "variant_voi"
