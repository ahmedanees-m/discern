"""autoPVS1 adapter - PVS1.

Wraps autoPVS1 (github.com/JiguangPeng/autopvs1, GPL-3.0) to assign PVS1 strength for
null variants per the ClinGen SVI null-variant decision tree, honoring gene-specific
VCEP overrides where defined. autoPVS1 needs a VEP cache and runs in
docker/Dockerfile.tools; the strength lookup is injectable for tests.
"""
from __future__ import annotations

from collections.abc import Callable

from core.adapter import EvidenceAdapter
from core.schemas import EvidenceContribution, PatientContext, Strength, Variant

# autoPVS1 strength label -> ACMG points.
_STRENGTH = {"Very Strong": Strength.PVS, "Strong": Strength.PS,
             "Moderate": Strength.PM, "Supporting": Strength.PP}


class AutoPVS1Adapter(EvidenceAdapter):
    code_group = ("PVS1",)
    version = "JiguangPeng@autopvs1-v2.0"

    def __init__(self, pvs1_lookup: Callable[[Variant], str | None] | None = None):
        # pvs1_lookup(variant) -> autoPVS1 strength label, or None if not a null variant
        self._lookup = pvs1_lookup

    def health_check(self) -> bool:
        return self._lookup is not None

    def evaluate(self, v: Variant, p: PatientContext) -> list[EvidenceContribution]:
        if self._lookup is None:
            return []
        label = self._lookup(v)
        if not label or label not in _STRENGTH:
            return []
        return [EvidenceContribution("PVS1", _STRENGTH[label], True, "autoPVS1",
                rationale=f"null variant; ClinGen SVI tree -> PVS1 {label}")]
