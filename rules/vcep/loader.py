"""VCEP spec loader.

Loads a machine-readable ClinGen gene-specific VCEP specification: per-code calibrated
strengths, allele-frequency thresholds, the disease/inheritance/mechanism, and a source
provenance string. Genes without a spec resolve to a ``covered=False`` reduced-confidence
record. The spec is consumed **per code** (via ``rules/vcep/partition.py``),
never as a bottom-line label.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from functools import cache

import yaml

from core.schemas import Strength
from rules.vcep.partition import is_variant_intrinsic

_SPEC_DIR = os.path.join(os.path.dirname(__file__), "specs")


@dataclass
class VcepSpec:
    spec_id: str
    genes: list[str]
    disease: str
    inheritance: str
    mechanism: str
    source: str
    prior_p: float = 0.10
    af: dict = field(default_factory=lambda: {"ba1": 0.05, "bs1": 0.0015, "pm2": 0.0001})
    code_strengths: dict = field(default_factory=dict)   # base code -> strength name
    computational: dict = field(default_factory=dict)    # pp3_revel/pp3_splice/bp4_revel/bp4_splice cut-offs
    covered: bool = True
    notes: str = ""

    def strength_for(self, code: str) -> Strength:
        """VCEP-specified strength for a base code, falling back to the ACMG default."""
        base = code.split("_", 1)[0]
        name = self.code_strengths.get(base)
        if name:
            return Strength[name]
        # ACMG defaults by prefix
        prefix = base[:3] if base[:3] in {"PVS", "BVS"} else base[:2]
        return {"PVS": Strength.PVS, "PS": Strength.PS, "PM": Strength.PM, "PP": Strength.PP,
                "BVS": Strength.BVS, "BS": Strength.BS, "BM": Strength.BM,
                "BA": Strength.BVS, "BP": Strength.BP}.get(prefix, Strength.PP)

    def af_threshold(self, kind: str) -> float:
        return float(self.af[kind])

    def variant_intrinsic_codes(self, codes: list[str]) -> list[str]:
        """Subset of applied codes owned by the variant-intrinsic factor."""
        return [c for c in codes if is_variant_intrinsic(c)]


def _reduced_confidence(gene: str) -> VcepSpec:
    return VcepSpec(
        spec_id=f"generic-acmg:{gene}", genes=[gene], disease="", inheritance="",
        mechanism="", source="no gene-specific VCEP specification - reduced confidence",
        covered=False, notes="Generic ACMG applied; output tagged reduced-confidence.",
    )


def load_spec_file(path: str) -> VcepSpec:
    with open(path, encoding="utf-8") as fh:
        d = yaml.safe_load(fh)
    return VcepSpec(
        spec_id=d["spec_id"], genes=d["genes"], disease=d.get("disease", ""),
        inheritance=d.get("inheritance", ""), mechanism=d.get("mechanism", ""),
        source=d.get("source", ""), prior_p=float(d.get("prior_p", 0.10)),
        af=d.get("af", {"ba1": 0.05, "bs1": 0.0015, "pm2": 0.0001}),
        code_strengths=d.get("code_strengths", {}),
        computational=d.get("computational", {}),
        covered=bool(d.get("covered", True)), notes=d.get("notes", ""),
    )


@cache
def get_spec(gene: str) -> VcepSpec:
    """Return the VCEP spec covering ``gene``, or a reduced-confidence record (G2)."""
    for fn in os.listdir(_SPEC_DIR):
        if not fn.endswith(".yaml"):
            continue
        spec = load_spec_file(os.path.join(_SPEC_DIR, fn))
        if gene in spec.genes:
            return spec
    return _reduced_confidence(gene)
