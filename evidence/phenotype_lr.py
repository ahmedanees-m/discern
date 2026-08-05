"""Phenotype likelihood ratios with pertinent negatives.

For each observed HPO feature, LR(term | disease) = freq(term | disease) / freq(term | bg).
A **pertinent negative** (an explicitly *absent* term, ``observed=False``) contributes
LR = (1 - freq(term | disease)) / (1 - freq(term | bg)) - so "no leukocytosis" argues
against LAD-III. Phenotype enters the joint model **once** (informing the disease); its
pull on the variant flows only through the disease->variant coupling (calibrated PP4).

``Disease.feature_lr[feature_id]`` holds ``(freq_in_disease, n_cases, pmid)`` - the
disease-conditional frequency (with provenance + sample size for the uncertainty layer).
"""
from __future__ import annotations

import math

from core.dx_schemas import Disease, Feature, FeatureKind

_BG_DEFAULT = 0.02      # background population frequency for an unspecified term


def feature_lr(freq_in_disease: float, observed: bool, bg: float = _BG_DEFAULT) -> float:
    """Present/absent likelihood ratio from a disease-conditional frequency."""
    freq = min(max(freq_in_disease, 1e-4), 1 - 1e-4)
    bg = min(max(bg, 1e-4), 1 - 1e-4)
    return (freq / bg) if observed else ((1 - freq) / (1 - bg))


def phenotype_loglik(features: list[Feature], disease: Disease) -> float:
    """log P(E_pheno | D) over clinical features the disease characterizes (present + absent)."""
    ll = 0.0
    for f in features:
        if f.kind not in (FeatureKind.CLINICAL, FeatureKind.LAB):
            continue
        entry = disease.feature_lr.get(f.id)
        if not entry:
            continue
        lr = feature_lr(float(entry[0]), f.observed)
        ll += math.log(max(lr, 1e-9))
    return ll
