"""ClinVar adapter - PS1 / PM5 (plan Step 1.5).

PS1: the variant's amino-acid change is already classified Pathogenic/Likely
     pathogenic in ClinVar (same AA change, established as damaging).
PM5: a *different* missense at the *same residue* is Pathogenic/Likely pathogenic,
     and the query change itself is not already established (novel missense at a
     known pathogenic residue).

Backed by a JSON index built from the kept ClinVar variant_summary
(``data/sources/build_clinvar_index.py``). The index path is configurable so the
adapter stays swappable and testable with a tiny in-memory fixture.
"""
from __future__ import annotations

import json
import os
from functools import lru_cache

from core.adapter import EvidenceAdapter
from core.schemas import EvidenceContribution, PatientContext, Strength, Variant
from data.sources.build_clinvar_index import parse_protein_change

_DEFAULT_INDEX = os.path.join(
    os.path.dirname(__file__), "..", "data", "processed", "clinvar_ps1_pm5.json"
)


@lru_cache(maxsize=4)
def _load_index(path: str) -> dict:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


class ClinVarAdapter(EvidenceAdapter):
    code_group = ("PS1", "PM5")
    version = "clinvar_20260503"

    def __init__(self, index_path: str | None = None, index: dict | None = None):
        self._path = index_path or _DEFAULT_INDEX
        self._index = index            # allow direct injection (tests)

    def _idx(self) -> dict:
        return self._index if self._index is not None else _load_index(self._path)

    def health_check(self) -> bool:
        if self._index is not None:
            return True
        return os.path.exists(self._path)

    def evaluate(self, v: Variant, p: PatientContext) -> list[EvidenceContribution]:
        if not v.hgvs_p:
            return []
        pc = parse_protein_change(v.hgvs_p)
        if not pc:
            return []
        ref, pos, alt = pc
        gene_idx = self._idx().get(v.gene, {})
        residue = gene_idx.get(f"{ref}{pos}", {})
        if not residue:
            return []
        src = f"ClinVar {self.version}"
        if alt in residue:
            return [EvidenceContribution(
                "PS1", Strength.PS, True, src,
                rationale=f"same AA change p.{ref}{pos}{alt} is {residue[alt]} in ClinVar")]
        # different pathogenic missense at the same residue -> PM5
        others = {a: s for a, s in residue.items() if a != alt}
        if others:
            example = next(iter(others))
            return [EvidenceContribution(
                "PM5", Strength.PM, True, src,
                rationale=f"different pathogenic missense p.{ref}{pos}{example} "
                          f"({residue[example]}) at the same residue")]
        return []
