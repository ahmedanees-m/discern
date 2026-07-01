"""ClinVar PS1/PM5 adapter (plan Step 1.5)."""
from adapters.clinvar import ClinVarAdapter
from core.adapter import EvidenceAdapter
from core.schemas import PatientContext, Variant
from data.sources.build_clinvar_index import parse_protein_change

_INDEX = {"F8": {"Arg2228": {"Gln": "Pathogenic", "Trp": "Likely pathogenic"}}}


def _v(hgvs_p, gene="F8"):
    return Variant("X", 154861759, "C", "T", gene=gene, hgvs_p=hgvs_p)


def test_is_adapter_and_health():
    a = ClinVarAdapter(index=_INDEX)
    assert isinstance(a, EvidenceAdapter)
    assert a.code_group == ("PS1", "PM5")
    assert a.health_check()


def test_ps1_same_change():
    a = ClinVarAdapter(index=_INDEX)
    out = a.evaluate(_v("p.Arg2228Gln"), PatientContext())
    assert len(out) == 1 and out[0].code == "PS1"
    assert out[0].strength.name == "PS"


def test_pm5_same_residue_different_change():
    a = ClinVarAdapter(index=_INDEX)
    out = a.evaluate(_v("p.Arg2228Leu"), PatientContext())   # Leu not in index, residue is
    assert len(out) == 1 and out[0].code == "PM5"
    assert out[0].strength.name == "PM"


def test_no_hit_returns_empty():
    a = ClinVarAdapter(index=_INDEX)
    assert a.evaluate(_v("p.Gly100Ser"), PatientContext()) == []          # unknown residue
    assert a.evaluate(_v("p.Arg2326Ter"), PatientContext()) == []         # nonsense -> not missense
    assert a.evaluate(Variant("X", 1, "C", "T", gene="F8"), PatientContext()) == []  # no hgvs_p


def test_protein_parser():
    assert parse_protein_change("NM_000132.4(F8):c.6682C>T (p.Arg2228Gln)") == ("Arg", 2228, "Gln")
    assert parse_protein_change("p.Arg2326Ter") is None          # nonsense
    assert parse_protein_change("p.Gly12=") is None              # synonymous
