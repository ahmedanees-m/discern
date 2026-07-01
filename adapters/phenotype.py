"""Phenotype adapter - PP4 (plan Step 1.7).

Computes semantic similarity between the patient's HPO terms and the gene-disease HPO
profile (Resnik IC over the HPO DAG, best-match-average between the two term sets), then
applies PP4 at the VCEP-specified strength. Recent ClinGen guidance allows PP4 up to +5
points and subtype-specific phenotypes (e.g. VWD type 2N) - so the strength is tiered.

The similarity backend runs against the kept HPO data (hp.obo / phenotype.hpoa /
genes_to_phenotype) - heavy graph work that lives behind an injectable callable, so the
adapter stays testable without loading the ontology. Notes-> HPO term extraction is an
LLM soft-task (llm/gateway, cloud Nemotron) validated against the ontology - the rule
engine, not the LLM, assigns the points.
"""
from __future__ import annotations

from collections.abc import Callable

from core.adapter import EvidenceAdapter
from core.schemas import EvidenceContribution, PatientContext, Strength, Variant

# Best-match-average similarity (0..1) -> PP4 strength tier.
PP4_SUPPORTING = 0.55
PP4_MODERATE = 0.75      # strong gene-phenotype / subtype match (e.g. VWD 2N)


class PhenotypeAdapter(EvidenceAdapter):
    code_group = ("PP4",)
    version = "HPO-Resnik-BMA"

    def __init__(self,
                 similarity_lookup: Callable[[list[str], str], float | None] | None = None,
                 moderate_threshold: float = PP4_MODERATE):
        # similarity_lookup(patient_hpo_terms, gene) -> similarity in [0,1], or None
        self._lookup = similarity_lookup
        self._moderate = moderate_threshold

    def health_check(self) -> bool:
        return self._lookup is not None

    def evaluate(self, v: Variant, p: PatientContext) -> list[EvidenceContribution]:
        if self._lookup is None or not p.hpo_terms:
            return []
        sim = self._lookup(p.hpo_terms, v.gene)
        if sim is None or sim < PP4_SUPPORTING:
            return []
        strength = Strength.PM if sim >= self._moderate else Strength.PP
        return [EvidenceContribution("PP4", strength, True, "HPO/Resnik",
                rationale=f"phenotype similarity {sim:.2f} -> PP4 {strength.name}")]
