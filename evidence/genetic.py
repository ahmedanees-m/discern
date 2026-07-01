"""Genetic evidence stream - variant-intrinsic likelihood (plan Phase 1.1).

Emits a likelihood over ``VariantState`` using **only the variant-intrinsic ACMG codes**
(PM2/BS1/BA1, PP3/BP4, PVS1, PM5/PS1, PM4/BP3, BP7) scored by the gene's VCEP spec.
PP4, PS3/BS3, PP1/PM3 are deliberately excluded here - they are owned by other factors
(disease coupling, functional, segregation/phasing), so each code enters the joint model
exactly once (the per-code circularity fix; Gate G3).

The likelihood is a Gaussian over the signed points: a variant truly in state V is
expected to accrue points near that band's center, so the observed points are a noisy
readout of V. This gives a smooth P(E_geno | V) for the joint factor graph.
"""
from __future__ import annotations

import math

from core.dx_schemas import VariantState
from rules.vcep.loader import VcepSpec

# Expected points center per variant state (Tavtigian bands: P>=10, LP 6-9, VUS 0-5,
# LB -1..-6, B<=-7).
_CENTERS = {
    VariantState.PATH: 12.0,
    VariantState.LP: 7.5,
    VariantState.VUS: 2.5,
    VariantState.LB: -3.5,
    VariantState.BEN: -9.0,
}


def variant_intrinsic_points(applied_codes: list[str], spec: VcepSpec) -> float:
    """Signed points from the variant-intrinsic codes only, scored by the VCEP spec."""
    intrinsic = spec.variant_intrinsic_codes(applied_codes)
    return sum(spec.strength_for(c).value for c in intrinsic)


def variant_intrinsic_likelihood(applied_codes: list[str], spec: VcepSpec,
                                 sigma: float = 4.0) -> dict[VariantState, float]:
    """P(E_geno | V) as a normalized distribution over VariantState (variant-intrinsic only)."""
    pts = variant_intrinsic_points(applied_codes, spec)
    raw = {s: math.exp(-((pts - c) ** 2) / (2 * sigma ** 2)) for s, c in _CENTERS.items()}
    z = sum(raw.values()) or 1.0
    return {s: v / z for s, v in raw.items()}
