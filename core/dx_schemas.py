"""DISCERN data model - the coupled disease x variant vocabulary (plan Section 3).

Extends the base core data model (`core/schemas.py`) with the disease-discrimination
types. The engine's output is named ``DxRecommendation`` to avoid colliding with the
single-variant ``engine.recommend.Recommendation``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class VariantState(Enum):
    PATH = "pathogenic"
    LP = "likely_path"
    VUS = "vus"
    LB = "likely_benign"
    BEN = "benign"


class FeatureKind(Enum):
    CLINICAL = "clinical"
    LAB = "lab"
    FUNCTIONAL = "functional"
    GENETIC = "genetic"


@dataclass
class Feature:
    id: str
    kind: FeatureKind
    value: object
    observed: bool = True        # observed=False => an explicit *pertinent negative*
    source: str = ""


@dataclass
class Observation:
    """A candidate next step - a test OR segregation OR phasing OR a functional assay."""

    id: str
    name: str
    kind: str                    # "lab" | "functional" | "segregation" | "phasing"
    informs: list[str]           # feature ids and/or "V" (variant) and/or "D" (disease)
    outcome_lr: dict             # outcome -> {disease_id | "V": LR}
    changes_management: bool = False
    accessibility: str = "moderate"   # high | moderate | low  (NO currency)


@dataclass
class Disease:
    id: str
    name: str
    genes: list[str]
    inheritance: str
    mechanism: str
    prior: float                                 # prevalence prior
    feature_lr: dict = field(default_factory=dict)   # feature_id -> (LR, n_cases, pmid)
    p_path_given_disease: float = 0.5            # P(variant pathogenic | this disease) = the PP4 coupling
    treatment: str = ""
    contraindications: list[str] = field(default_factory=list)
    vcep_spec: str | None = None              # rules/vcep/... or None (-> reduced-confidence tag)


@dataclass
class DiscriminationCluster:
    id: str
    name: str
    diseases: list[Disease]
    discriminating_features: list[str] = field(default_factory=list)
    next_observations: list[Observation] = field(default_factory=list)


@dataclass
class JointPosterior:
    cluster: DiscriminationCluster
    p_disease: dict                              # disease_id -> (prob, lo, hi)  credible interval
    p_variant: dict                              # variant_id -> {VariantState: prob}
    leading: str
    confidence: float
    decided: bool                                # abstain if not decided


@dataclass
class SafetyFlag:
    leading_id: str
    competitor_id: str
    management_divergence: str                   # what changes (e.g. "DDAVP caution if 2B")
    p_competitor: float
    severity: str
    resolving_observation: Observation | None
    message: str


@dataclass
class DxRecommendation:
    """The engine's output (plan Section 3 `Recommendation`)."""

    posterior: JointPosterior
    reclassified_variants: dict = field(default_factory=dict)   # variant_id -> (old, new, drivers)
    safety_flags: list[SafetyFlag] = field(default_factory=list)
    next_observation: Observation | None = None
    explanation: str = ""
    audit: dict = field(default_factory=dict)
