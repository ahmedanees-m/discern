"""Ancestry-stratified resolution-gap evaluation (plan Section A.5, Step 6.2).

Measures the gap in resolution rate across genetic-similarity groups, with bootstrap
CIs, and compares the equity-aware config vs the non-equity-aware one. Target: the
equity-aware routing NARROWS the gap.
"""
from __future__ import annotations

from dataclasses import dataclass

from equity.dashboard import cohort_resolution_gap
from eval.stats import bootstrap_ci


@dataclass
class EquityComparison:
    gap_equity_off: float
    gap_equity_on: float
    narrowed: bool
    rates_off: dict
    rates_on: dict


def resolution_gap(resolved_by_group: dict[str, tuple[int, int]]) -> dict:
    return cohort_resolution_gap(resolved_by_group)


def compare_equity_configs(off: dict[str, tuple[int, int]],
                           on: dict[str, tuple[int, int]]) -> EquityComparison:
    g_off = cohort_resolution_gap(off)
    g_on = cohort_resolution_gap(on)
    return EquityComparison(
        gap_equity_off=g_off["resolution_gap"], gap_equity_on=g_on["resolution_gap"],
        narrowed=g_on["resolution_gap"] < g_off["resolution_gap"],
        rates_off=g_off["rates"], rates_on=g_on["rates"],
    )


def group_rate_ci(resolved_flags: list[int]) -> tuple[float, float, float]:
    """Mean resolution rate + bootstrap CI for one group's per-case 0/1 outcomes."""
    return bootstrap_ci([float(x) for x in resolved_flags])
