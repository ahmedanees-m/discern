"""Lab / functional evidence - the one factor that touches BOTH disease and variant
(plan Phase 1.3 / Section 3.1).

A functional result (flow %, RIPA pattern, aggregometry, factor activity, smear)
informs the disease (it fits one disease's predicted pattern) AND the variant (a
damaging functional result is PS3-grade evidence the variant is pathogenic). It is one
factor `P(E_func | D, V)`, so it is counted **once** - never added separately as PS3 and
as a disease feature. The variant side uses Brnich-2020 OddsPath-style strength.
"""
from __future__ import annotations

import math

from core.dx_schemas import Disease, Feature, FeatureKind, VariantState
from evidence.phenotype_lr import feature_lr

# Functional-result strength on the variant (Brnich 2020 OddsPath, log-odds toward path).
_FUNC_ODDSPATH_DAMAGING = 18.7      # strong functional evidence of damage -> pathogenic
_FUNC_ODDSPATH_NORMAL = 1 / 18.7    # normal function -> benign


def _is_damaging(feature: Feature) -> bool:
    """A functional readout indicating loss/abnormality (absent protein, low activity)."""
    v = feature.value
    if isinstance(v, bool):
        return v
    return str(v).lower() in {"absent", "abnormal", "damaging", "reduced", "low", "loss", "deficient"}


def _functional_variant_loglik(feature: Feature, variant_state: VariantState) -> float:
    """Variant side of a functional result: damaging -> supports PATH, normal -> BEN."""
    op = _FUNC_ODDSPATH_DAMAGING if _is_damaging(feature) else _FUNC_ODDSPATH_NORMAL
    # likelihood that a variant in this state would produce this functional result
    path_like = {VariantState.PATH: 1.0, VariantState.LP: 0.8, VariantState.VUS: 0.5,
                 VariantState.LB: 0.2, VariantState.BEN: 0.0}[variant_state]
    # interpolate the odds with how pathogenic the state is
    odds = op ** (2 * path_like - 1)        # path_like=1 -> op; 0 -> 1/op; 0.5 -> 1
    return math.log(max(odds, 1e-9))


def func_loglik(features: list[Feature], disease: Disease, variant_state: VariantState) -> float:
    """log P(E_func | D, V): disease-pattern LR + (for FUNCTIONAL) the variant PS3 effect."""
    ll = 0.0
    for f in features:
        if f.kind not in (FeatureKind.LAB, FeatureKind.FUNCTIONAL):
            continue
        entry = disease.feature_lr.get(f.id)
        if entry:                                   # disease side (pattern match)
            ll += math.log(max(feature_lr(float(entry[0]), f.observed), 1e-9))
        if f.kind == FeatureKind.FUNCTIONAL:        # variant side (PS3/BS3), counted here once
            ll += _functional_variant_loglik(f, variant_state)
    return ll
