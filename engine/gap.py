"""Evidence-gap & conflict analysis (plan Phase 2.1 / design Section 4.1).

A VUS is an evidence deficit (not enough points either way) or an evidence conflict
(opposing codes). This module:
  * enumerates the *attainable* not-yet-applied codes (context- and mechanism-aware),
  * computes the points gap to a target classification, and
  * flags VUS-by-conflict (opposing applied codes) so the recommender targets a
    disambiguating action rather than just any point-adder.

Attainability is mechanism- and context-gated: PS2 only if parents are available;
PP1/BS4 only with an informative family; PM3 (recessive) needs phasing; PS3/BS3
only if a validated assay exists for the gene's mechanism.
"""
from __future__ import annotations

from dataclasses import dataclass

from core.schemas import CodeOpportunity, PatientContext, PointsLedger, Strength
from engine.action_map import actions_for_code
from rules.point_engine import BANDS, Classification

# Acquirable codes an action can produce, with default strength and prerequisite flags.
# prereq names are attributes on PatientContext (or the sentinel "assay_exists",
# resolved against the domain action catalog for the variant's mechanism).
ACQUIRABLE = {
    "PS3": (Strength.PS, ["assay_exists"], "Functional assay supports damaging effect"),
    "BS3": (Strength.BS, ["assay_exists"], "Functional assay shows no damaging effect"),
    "PS2": (Strength.PS, ["parents_available"], "De novo (parental/trio testing)"),
    "PM3": (Strength.PM, ["parents_available"], "In trans with pathogenic (phasing)"),
    "PP1": (Strength.PP, ["informative_family"], "Co-segregation with disease"),
    "BS4": (Strength.BS, ["informative_family"], "Non-segregation with disease"),
    "PP4": (Strength.PP, [], "Phenotype highly specific for the gene"),
    "PS4": (Strength.PS, [], "Case enrichment / prevalence in affected"),
}


@dataclass
class ConflictReport:
    in_conflict: bool
    pathogenic_codes: list[str]
    benign_codes: list[str]


def _applied(ledger: PointsLedger) -> set[str]:
    # normalize to base code (strip strength suffix)
    return {c.code.split("_", 1)[0] for c in ledger.contributions if c.applied}


def attainable_codes(ledger: PointsLedger, patient: PatientContext, mechanism: str,
                     domain: str = "bleeding") -> list[CodeOpportunity]:
    applied = _applied(ledger)
    out: list[CodeOpportunity] = []
    for code, (strength, prereqs, rationale) in ACQUIRABLE.items():
        if code in applied:
            continue
        missing = []
        ok = True
        for pr in prereqs:
            if pr == "assay_exists":
                if not actions_for_code(code, mechanism, domain):
                    ok = False
            elif not getattr(patient, pr, False):
                ok = False
                missing.append(pr)
        if not ok and not missing:
            continue          # hard-blocked (e.g. no assay for this mechanism)
        out.append(CodeOpportunity(
            code=code,
            strength_if_pos=strength,
            delta_points=float(strength.value),
            prerequisites=missing,
            rationale=rationale,
        ))
    return out


def gap_to_target(ledger: PointsLedger, target: Classification = Classification.LP) -> float:
    """Points still needed to reach `target` (0 if already there)."""
    threshold = next(t for t, c in BANDS if c is target)
    return max(0.0, threshold - ledger.points)


def detect_conflict(ledger: PointsLedger) -> ConflictReport:
    path = [c.code for c in ledger.contributions if c.applied and c.strength.value > 0]
    ben = [c.code for c in ledger.contributions if c.applied and c.strength.value < 0]
    return ConflictReport(bool(path and ben), path, ben)


# Modalities that provide *orthogonal* evidence - used to break VUS-by-conflict by
# acquiring an independent readout rather than stacking more of the same (design Section 6).
ORTHOGONAL_MODALITIES = {"functional", "rna", "segregation"}


def case_note(inheritance: str) -> str:
    """Inheritance-specific guidance (recessive diplotype / X-linked hemizygote, Section 3.1/Section 6)."""
    inh = (inheritance or "").upper()
    if inh.startswith("AR"):
        return ("Recessive: disease causality is a property of the *genotype* - also pursue "
                "the in-trans relationship (PM3 via phasing/parental testing) and functional "
                "protein loss, which together support the diplotype even before the second "
                "allele is fully resolved.")
    if inh.startswith("XL"):
        return ("X-linked: a hemizygous male needs one pathogenic allele; the CFD-VCEP "
                "hemizygote PS4 rule applies (F8/F9).")
    return ""
