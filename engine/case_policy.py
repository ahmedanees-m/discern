"""Whole-odyssey case policy (plan Phase 4 / design Section 4.5).

Lifts the engine from single-variant to whole-case decisions. Routes by failure mode:
  * interpretation failure (a VUS candidate exists) -> variant-level VOI (Phase 2);
  * detection failure (no candidate / missing second hit / hard region) -> modality
    escalation (long-read / optical / RNA / methylation);
and decides whether to stop/report, wait (reanalysis), or matchmake.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from engine.odyssey_actions import MODALITY_CATALOG, ModalityAction
from rules.point_engine import Classification


@dataclass
class CaseState:
    has_vus_candidate: bool          # an interpretable VUS exists (interpretation failure)
    has_any_candidate: bool          # any plausible candidate at all
    suspected_splicing: bool = False
    suspected_structural: bool = False
    suspected_imprinting: bool = False
    novel_gene_candidate: bool = False
    modalities_done: set[str] = field(default_factory=set)


@dataclass
class CasePlan:
    route: str                       # "variant_voi" | "modality_escalation" | "reanalysis" | "matchmaking" | "report"
    modality: ModalityAction | None
    rationale: str


def _best_modality(names: list[str], done: set[str]) -> ModalityAction | None:
    cands = [a for a in MODALITY_CATALOG if a.name in names and a.name not in done]
    return max(cands, key=lambda a: a.yield_prior / max(a.cost_usd, 1)) if cands else None


def decide_case(state: CaseState) -> CasePlan:
    """Pick the next whole-case step (variant resolution vs modality escalation vs wait)."""
    # Interpretation failure: an interpretable VUS exists -> resolve it (Phase 2 VOI).
    if state.has_vus_candidate:
        return CasePlan("variant_voi", None,
                        "Interpretable VUS present -> resolve via variant-level value-of-information.")
    # Detection failure: escalate the modality matched to the suspected mechanism.
    if state.suspected_splicing:
        m = _best_modality(["rna_seq"], state.modalities_done)
        if m:
            return CasePlan("modality_escalation", m, "Suspected splicing -> RNA-seq.")
    if state.suspected_structural:
        m = _best_modality(["long_read_genome", "optical_genome_mapping"], state.modalities_done)
        if m:
            return CasePlan("modality_escalation", m, "Suspected structural/hard region -> long-read/optical.")
    if state.suspected_imprinting:
        m = _best_modality(["methylation_episignature"], state.modalities_done)
        if m:
            return CasePlan("modality_escalation", m, "Suspected imprinting -> methylation episignature.")
    if not state.has_any_candidate:
        # generic detection escalation, cheapest-yielding modality not yet done
        m = _best_modality(["rna_seq", "long_read_genome", "optical_genome_mapping",
                            "methylation_episignature"], state.modalities_done)
        if m:
            return CasePlan("modality_escalation", m,
                            "No candidate yet -> escalate the highest-yield-per-cost modality.")
    if state.novel_gene_candidate:
        return CasePlan("matchmaking", _best_modality(["matchmaking"], set()),
                        "Novel-gene candidate -> matchmaking (GeneMatcher / MME).")
    # Nothing actionable now -> 'wait/reanalyze'.
    return CasePlan("reanalysis", _best_modality(["reanalysis"], set()),
                    "Not currently resolvable -> scheduled reanalysis as the knowledge base grows.")


def already_actionable(current_class: Classification) -> bool:
    return current_class in (Classification.P, Classification.LP, Classification.B, Classification.LB)
