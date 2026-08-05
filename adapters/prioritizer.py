"""Prioritizer adapter - candidate belief state.

Normalizes the output of upstream gene/variant prioritizers (Exomiser, AI-MARRVEL,
optionally DeepRare) into a common ``CandidateSet``, and records cross-tool
disagreement - which is a *feature*: disagreement raises a human-review flag rather
than being silently accepted. These tools run in their own images / services; the
per-tool ranking is injected here.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Candidate:
    gene: str
    score: float
    source: str


@dataclass
class CandidateSet:
    candidates: list[Candidate] = field(default_factory=list)
    disagreement_flag: bool = False
    notes: str = ""


def normalize(tool_rankings: dict[str, list[tuple[str, float]]],
              top_k: int = 10) -> CandidateSet:
    """Merge per-tool (gene, score) rankings; flag cross-tool top-gene disagreement."""
    cands: list[Candidate] = []
    top_genes = []
    for tool, ranking in tool_rankings.items():
        for gene, score in ranking[:top_k]:
            cands.append(Candidate(gene, score, tool))
        if ranking:
            top_genes.append(ranking[0][0])
    disagree = len(set(top_genes)) > 1
    note = "cross-tool top-gene disagreement -> human review" if disagree else ""
    return CandidateSet(cands, disagreement_flag=disagree, notes=note)
