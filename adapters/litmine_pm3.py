"""Literature-mining adapter - PM3 / PS4.

AutoPM3 (Bioinformatics 2025, PMC12263107) extracts PM3 evidence (in trans with a
pathogenic variant) from literature via LLM + RAG; RAG over literature also supports
PS4 (case enrichment). **All LLM output is a proposal carrying citations; the rule
engine assigns the points - never the LLM** (the project's core safety property).

The mining backend is injectable; in production it routes through llm/gateway (cloud
Nemotron) with required citations, logged to the audit ledger.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from core.adapter import EvidenceAdapter
from core.schemas import EvidenceContribution, PatientContext, Strength, Variant


@dataclass
class MinedEvidence:
    code: str                 # "PM3" | "PS4"
    strength: Strength
    citations: list[str]      # PMIDs - REQUIRED; no citation -> not applied
    rationale: str = ""


class LitMineAdapter(EvidenceAdapter):
    code_group = ("PM3", "PS4")
    version = "AutoPM3-2025+RAG"

    def __init__(self, miner: Callable[[Variant, PatientContext], list[MinedEvidence]] | None = None):
        self._miner = miner

    def health_check(self) -> bool:
        return self._miner is not None

    def evaluate(self, v: Variant, p: PatientContext) -> list[EvidenceContribution]:
        if self._miner is None:
            return []
        out = []
        for m in self._miner(v, p):
            if not m.citations:                # ungrounded LLM proposal -> reject
                continue
            out.append(EvidenceContribution(
                m.code, m.strength, True, "AutoPM3/RAG (LLM-proposed, rule-assigned)",
                rationale=f"{m.rationale} [PMIDs: {','.join(m.citations)}]"))
        return out
