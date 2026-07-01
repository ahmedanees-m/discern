"""MaveDB adapter - PS3 / BS3 functional (plan Step 1.6).

Looks up the variant in any MaveDB score set for the gene and maps the functional
score to ACMG strength via OddsPath calibration (Brnich et al. 2020, Genome Medicine):
  OddsPath >= 18.7 -> Strong, >= 4.3 -> Moderate, >= 2.1 -> Supporting (benign inverse:
  <= 0.053 -> Strong, <= 0.23 -> Moderate, <= 0.48 -> Supporting).

Backed by the MaveDB API (`pip install mavedb`); the lookup is injectable for tests.
Where saturation data exists (e.g. BRCA1 SGE) this can resolve a VUS directly.
"""
from __future__ import annotations

from collections.abc import Callable

from core.adapter import EvidenceAdapter
from core.schemas import EvidenceContribution, PatientContext, Strength, Variant

# Brnich 2020 OddsPath -> evidence strength.
_PATH = [(18.7, Strength.PS), (4.3, Strength.PM), (2.1, Strength.PP)]
_BENIGN = [(0.053, Strength.BS), (0.23, Strength.BM), (0.48, Strength.BP)]


def oddspath_to_strength(oddspath: float):
    for thr, st in _PATH:
        if oddspath >= thr:
            return "PS3", st
    for thr, st in _BENIGN:
        if oddspath <= thr:
            return "BS3", st
    return None


class MaveAdapter(EvidenceAdapter):
    code_group = ("PS3", "BS3")
    version = "MaveDB-2024"

    def __init__(self, oddspath_lookup: Callable[[Variant], float | None] | None = None):
        self._lookup = oddspath_lookup

    def health_check(self) -> bool:
        return self._lookup is not None

    def evaluate(self, v: Variant, p: PatientContext) -> list[EvidenceContribution]:
        if self._lookup is None:
            return []
        op = self._lookup(v)
        if op is None:
            return []
        hit = oddspath_to_strength(op)
        if not hit:
            return []
        code, strength = hit
        return [EvidenceContribution(code, strength, True, "MaveDB",
                rationale=f"MAVE OddsPath={op:.2f} -> {code} {strength.name} (Brnich2020)")]
