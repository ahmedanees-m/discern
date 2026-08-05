"""Core data model for DISCERN.

Defined once and shared by the evidence adapters, the deterministic rule engine, and the
value-of-information decision core.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class Strength(Enum):
    """Signed ACMG/AMP evidence points (Tavtigian 2020, naturally scaled)."""

    PVS = 8
    PS = 4
    PM = 2
    PP = 1
    BP = -1
    BM = -2
    BS = -4
    BVS = -8


class Classification(Enum):
    P = "Pathogenic"
    LP = "Likely pathogenic"
    VUS = "Uncertain"
    LB = "Likely benign"
    B = "Benign"


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    build: str = "GRCh38"
    transcript: str | None = None

    def key(self) -> str:
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"


@dataclass
class PatientContext:
    hpo_terms: list[str] = field(default_factory=list)
    ancestry_group: str | None = None          # inferred (equity module)
    ancestry_confidence: float = 0.0
    parents_available: bool = False
    informative_family: bool = False              # for segregation (PP1/BS4)
    existing_assays: dict = field(default_factory=dict)   # already-done results
    sex: str | None = None


@dataclass
class EvidenceContribution:
    code: str                  # e.g. "PS3", "PM2", "PP3"
    strength: Strength
    applied: bool              # True = currently met; False = a potential opportunity
    source: str                # provenance: which adapter/DB/version
    reliability: float = 1.0   # equity down-weighting (0..1)
    rationale: str = ""


@dataclass
class PointsLedger:
    variant: Variant
    contributions: list[EvidenceContribution] = field(default_factory=list)
    spec_version: str = ""

    @property
    def points(self) -> float:
        """Signed points sum, with equity reliability weighting.

        Only currently-applied codes count. BA1 is handled as a standalone
        benign pre-filter by the rule engine and is excluded from this sum.
        """
        return sum(
            c.strength.value * c.reliability
            for c in self.contributions
            if c.applied
        )


@dataclass
class CodeOpportunity:
    """An attainable, not-yet-met evidence code (output of the gap engine)."""

    code: str
    strength_if_pos: Strength
    delta_points: float
    prerequisites: list[str]     # e.g. ["parents_available"]
    rationale: str


@dataclass
class Action:
    """A real-world orderable test/assay/modality."""

    name: str
    yields_codes: list[str]
    cost_usd: float
    turnaround_days: int
    sensitivity: float           # priors (literature -> learning loop)
    specificity: float
    modality: str
    burden: float = 0.0          # patient burden / invasiveness penalty


@dataclass
class RankedPlan:
    actions: list[tuple[Action, float]]      # (action, value_per_cost)
    explanation: str
    current_class: Classification
    current_posterior: float
