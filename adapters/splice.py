"""Splice adapter - PP3 / BP4 / BP7 + RNA-resolvable flag (plan Step 1.3).

Scores splice impact with Pangolin (github.com/tkzeng/Pangolin, GPL-3.0) against the
VCEP cutoffs, and - crucially for the synonymous-but-not-silent ITGB3 class of case -
flags variants where an RNA assay would yield PS3-splice (the gap engine routes to
RNA-seq/RT-PCR). Pangolin runs in docker/Dockerfile.tools (NOT `pip install pangolin`,
which is the SARS-CoV-2 lineage tool); the score lookup is injectable for tests.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from core.adapter import EvidenceAdapter
from core.schemas import EvidenceContribution, PatientContext, Strength, Variant

# Pangolin splice-score cutoffs (per current VCEP guidance; re-verify at execution).
PP3_CUTOFF = 0.5      # high splice-altering score -> PP3
BP4_CUTOFF = 0.1      # low score at a splice-relevant site -> BP4/BP7
RNA_RESOLVABLE_CUTOFF = 0.2   # ambiguous score near a splice region -> recommend RNA


@dataclass
class SpliceResult:
    score: float
    splice_relevant: bool         # in/near a splice region


class SpliceAdapter(EvidenceAdapter):
    code_group = ("PP3", "BP4", "BP7")
    version = "Pangolin"

    def __init__(self, splice_lookup: Callable[[Variant], SpliceResult | None] | None = None,
                 pp3_cutoff: float = PP3_CUTOFF, bp4_cutoff: float = BP4_CUTOFF,
                 rna_cutoff: float = RNA_RESOLVABLE_CUTOFF):
        self._lookup = splice_lookup
        self._pp3 = pp3_cutoff
        self._bp4 = bp4_cutoff
        self._rna = rna_cutoff

    @classmethod
    def for_spec(cls, spec, splice_lookup: Callable[[Variant], SpliceResult | None] | None = None) -> SpliceAdapter:
        """Build a spec-aware splice adapter from a VcepSpec's `computational` cut-offs."""
        c = getattr(spec, "computational", {}) or {}
        return cls(splice_lookup, c.get("pp3_splice", PP3_CUTOFF), c.get("bp4_splice", BP4_CUTOFF))

    def health_check(self) -> bool:
        return self._lookup is not None

    def evaluate(self, v: Variant, p: PatientContext) -> list[EvidenceContribution]:
        if self._lookup is None:
            return []
        r = self._lookup(v)
        if r is None:
            return []
        out = []
        if r.score >= self._pp3:
            out.append(EvidenceContribution("PP3", Strength.PP, True, "Pangolin",
                       rationale=f"splice score {r.score:.2f} >= {self._pp3}"))
        elif r.splice_relevant and r.score <= self._bp4:
            out.append(EvidenceContribution("BP4", Strength.BP, True, "Pangolin",
                       rationale=f"splice score {r.score:.2f} <= {self._bp4} at splice site"))
        return out

    def rna_resolvable(self, v: Variant) -> bool:
        """Flag for the gap engine: an RNA assay (-> PS3-splice) would be informative."""
        if self._lookup is None:
            return False
        r = self._lookup(v)
        return bool(r and r.splice_relevant and r.score >= self._rna)
