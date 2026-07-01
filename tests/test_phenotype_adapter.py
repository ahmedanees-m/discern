"""Phenotype PP4 adapter (plan Step 1.7)."""
from adapters.phenotype import PhenotypeAdapter
from core.adapter import EvidenceAdapter
from core.schemas import PatientContext, Strength, Variant


def _v(gene="VWF"):
    return Variant("12", 6000000, "C", "T", gene=gene, hgvs_p="p.Arg854Gln")


def test_is_adapter():
    a = PhenotypeAdapter(similarity_lookup=lambda terms, gene: 0.6)
    assert isinstance(a, EvidenceAdapter)
    assert a.code_group == ("PP4",) and a.health_check()


def test_pp4_supporting_and_moderate():
    p = PatientContext(hpo_terms=["HP:0011890", "HP:0000967"])      # bleeding phenotype
    sup = PhenotypeAdapter(similarity_lookup=lambda t, g: 0.6).evaluate(_v(), p)
    assert sup[0].code == "PP4" and sup[0].strength is Strength.PP
    mod = PhenotypeAdapter(similarity_lookup=lambda t, g: 0.85).evaluate(_v(), p)   # 2N-strong match
    assert mod[0].strength is Strength.PM


def test_no_call_below_threshold_or_no_terms():
    p = PatientContext(hpo_terms=["HP:0011890"])
    assert PhenotypeAdapter(similarity_lookup=lambda t, g: 0.3).evaluate(_v(), p) == []
    # no HPO terms -> nothing
    assert PhenotypeAdapter(similarity_lookup=lambda t, g: 0.9).evaluate(_v(), PatientContext()) == []
