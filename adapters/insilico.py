"""Calibrated in-silico adapter - PP3 / BP4 (plan Step 1.2, Section A.3).

Maps a single predefined predictor's score to ACMG strength via the *calibrated* local
posteriors of Pejaver et al. 2022 (AJHG 109:2163) - not a flat "Supporting". One tool
is chosen per VCEP rule, per ClinGen PP3/BP4 guidance.

Verified REVEL thresholds (Pejaver 2022, Table 2):
  PP3  Supporting >= 0.644, Moderate >= 0.773, Strong >= 0.932
  BP4  Supporting <= 0.290, Moderate <= 0.183, Strong <= 0.016 (Very Strong <= 0.003)

The score lookup is injectable (a dbNSFP/AlphaMissense reader runs in a tool image);
the calibration logic here is exact and unit-tested.
"""
from __future__ import annotations

from collections.abc import Callable

from core.adapter import EvidenceAdapter
from core.schemas import EvidenceContribution, PatientContext, Strength, Variant

# Calibrated thresholds per predictor. Each is an ordered list of
# (low, high, code, strength); a score in [low, high) maps to that code/strength.
THRESHOLDS = {
    "REVEL": {
        "pp3": [(0.644, 0.773, Strength.PP), (0.773, 0.932, Strength.PM), (0.932, 1.01, Strength.PS)],
        "bp4": [(0.183, 0.290, Strength.BP), (0.016, 0.183, Strength.BM), (-0.01, 0.016, Strength.BS)],
    },
}


def revel_strength(score: float) -> tuple[str, Strength] | None:
    t = THRESHOLDS["REVEL"]
    for lo, hi, st in t["pp3"]:
        if lo <= score < hi:
            return "PP3", st
    for lo, hi, st in t["bp4"]:
        if lo <= score < hi:
            return "BP4", st
    return None


class InSilicoAdapter(EvidenceAdapter):
    code_group = ("PP3", "BP4")
    version = "REVEL/Pejaver2022"

    def __init__(self, predictor: str = "REVEL",
                 score_lookup: Callable[[Variant], float | None] | None = None,
                 pp3_revel: float | None = None, pp3_revel_moderate: float | None = None,
                 bp4_revel: float | None = None):
        self.predictor = predictor
        self._lookup = score_lookup
        # When the gene-specific VCEP cut-offs are supplied, use them; else fall back to the
        # Pejaver2022 calibrated bands. This is the forward-path wiring of the spec computational block.
        self._pp3 = pp3_revel
        self._pp3_mod = pp3_revel_moderate
        self._bp4 = bp4_revel

    @classmethod
    def for_spec(cls, spec, score_lookup: Callable[[Variant], float | None] | None = None,
                 predictor: str = "REVEL") -> InSilicoAdapter:
        """Build a spec-aware in-silico adapter from a VcepSpec's `computational` cut-offs."""
        c = getattr(spec, "computational", {}) or {}
        return cls(predictor, score_lookup, c.get("pp3_revel"), c.get("pp3_revel_moderate"),
                   c.get("bp4_revel"))

    def _spec_strength(self, score: float) -> tuple[str, Strength] | None:
        if self._pp3_mod is not None and score >= self._pp3_mod:
            return "PP3", Strength.PM
        if self._pp3 is not None and score >= self._pp3:
            return "PP3", Strength.PP
        if self._bp4 is not None and score <= self._bp4:
            return "BP4", Strength.BP
        return None

    @property
    def _spec_aware(self) -> bool:
        return self._pp3 is not None or self._bp4 is not None

    def health_check(self) -> bool:
        return self._lookup is not None and (self._spec_aware or self.predictor in THRESHOLDS)

    def evaluate(self, v: Variant, p: PatientContext) -> list[EvidenceContribution]:
        if self._lookup is None:
            return []
        score = self._lookup(v)
        if score is None:
            return []
        if self._spec_aware:
            hit = self._spec_strength(score)
            src = f"{self.predictor} (VCEP cut-off)"
        else:
            hit = revel_strength(score) if self.predictor == "REVEL" else None
            src = f"{self.predictor} (Pejaver2022 calibrated)"
        if not hit:
            return []
        code, strength = hit
        return [EvidenceContribution(code, strength, True, src,
                rationale=f"{self.predictor}={score:.3f} -> {code} {strength.name}")]
