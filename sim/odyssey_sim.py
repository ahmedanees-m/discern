"""Parametric odyssey simulator + baseline policies.

Generates synthetic VUS cases with a *known* ground-truth label and a set of
attainable, costed actions (assay sens/spec/cost from the bleeding catalog). A
policy picks actions sequentially; each action returns a result drawn from its
sens/spec given the true label; the rule engine updates the points; the case
resolves when the posterior crosses a reporting threshold (P/LP or B/LB) or the
budget is exhausted.

Policies compared:
  * voi        - value-of-information ranking (delta_utility per cost), keeps testing
  * voi_stop   - VOI with graceful degradation (stop when no action improves)
  * greedy     - best sensitivity per dollar (cost-aware; a strong baseline)
  * random     - uniform random action
  * fixed      - a fixed clinical default order

Finding (homogeneous functional-assay cohort): VOI clearly beats the naive
random baseline (cheaper + more accurate) and uses no more tests than greedy;
voi_stop achieves the lowest cost by not testing the unresolvable. VOI does NOT
dominate a cost-aware greedy on dollars here - an expected, reported
result when one cheap accurate assay suffices; VOI's belief-state / conflict /
mechanism-applicability advantages are expected to matter on heterogeneous and
retrospective-journey cohorts. Pure-Python (seeded) so it runs on the laptop / CI.
"""
from __future__ import annotations

import random
from dataclasses import dataclass

from core.schemas import Action
from engine.voi import value_density
from rules.posterior import posterior

# Reporting thresholds (posterior).
_HI, _LO = 0.90, 0.10


@dataclass
class SimAction:
    action: Action
    k: float = 4.0          # points if positive (PS3-class)


@dataclass
class SimCase:
    start_points: float
    is_pathogenic: bool
    actions: list[SimAction]


@dataclass
class SimResult:
    resolved: bool
    correct: bool
    cost_usd: float
    days: int
    n_tests: int


def _draw_result(rng: random.Random, a: Action, is_path: bool) -> bool:
    """Did the assay return a 'positive' (damaging) result, given the truth?"""
    return rng.random() < (a.sensitivity if is_path else (1 - a.specificity))


def _pick(policy: str, points: float, acts: list[SimAction], rng: random.Random,
          prior_p: float = 0.5, stop_when_futile: bool = False) -> SimAction | None:
    if policy == "random":
        return rng.choice(acts)
    if policy == "fixed":
        return acts[0]                                   # fixed default order
    if policy == "greedy":
        return max(acts, key=lambda sa: sa.action.sensitivity / max(sa.action.cost_usd, 1))
    # voi: rank by decision-utility per cost, using each action's TRUE point yield and
    # the current belief state (points) - this state-awareness is VOI's edge over greedy.
    # `stop_when_futile` enables the "graceful degradation" behavior (don't keep testing
    # when nothing improves the expected actionability); off for a like-for-like baseline.
    scored = [(sa, value_density(points, sa.action, objective="delta_utility",
                                 k_pos=sa.k, k_neg=sa.k, prior_p=prior_p)) for sa in acts]
    best = max(scored, key=lambda t: t[1][0])
    if stop_when_futile and best[1][1]["delta_utility"] <= 0:
        return None
    return best[0]


def run_case(case: SimCase, policy: str, rng: random.Random, budget_usd: float = 5000,
             prior_p: float = 0.5) -> SimResult:
    # "voi_stop" = VOI with graceful degradation (stop when no action improves);
    # "voi" = VOI that keeps testing, for a like-for-like comparison with greedy.
    base_policy = "voi" if policy == "voi_stop" else policy
    stop_when_futile = policy == "voi_stop"
    points = case.start_points
    cost = days = n = 0
    remaining = list(case.actions)
    while remaining and cost < budget_usd:
        sa = _pick(base_policy, points, remaining, rng, prior_p, stop_when_futile)
        if sa is None:                # VOI graceful degradation -> stop
            break
        remaining.remove(sa)
        cost += sa.action.cost_usd
        days += sa.action.turnaround_days
        n += 1
        positive = _draw_result(rng, sa.action, case.is_pathogenic)
        points += sa.k if positive else -sa.k
        p = posterior(points, prior_p)
        if p >= _HI or p <= _LO:
            called_path = p >= _HI
            return SimResult(True, called_path == case.is_pathogenic, cost, days, n)
    return SimResult(False, False, cost, days, n)


def make_actions() -> list[SimAction]:
    """A representative action set (functional + segregation), literature-seeded priors."""
    return [
        # cheap, high-sensitivity but LOW-yield (PP, +1) in-silico screen: greedy's
        # "sensitivity per dollar" over-prefers it, but it rarely resolves a case alone.
        SimAction(Action("cheap_screen", ["PP3"], 60, 2, 0.95, 0.85, "insilico"), k=1.0),
        SimAction(Action("factor_activity", ["PS3", "BS3"], 250, 5, 0.92, 0.95, "functional")),
        SimAction(Action("flow_cytometry", ["PS3"], 450, 7, 0.90, 0.92, "functional")),
        SimAction(Action("rna_seq_splice", ["PS3", "BS3"], 900, 28, 0.78, 0.92, "rna")),
        SimAction(Action("thrombin_gen", ["PS3"], 400, 10, 0.78, 0.88, "functional")),
        SimAction(Action("family_seg", ["PP1", "BS4"], 800, 21, 0.65, 0.95, "segregation"), k=1.0),
        SimAction(Action("trio_testing", ["PS2", "PM3"], 1200, 21, 0.72, 0.99, "segregation"), k=2.0),
    ]


def make_cohort(n: int, rng: random.Random) -> list[SimCase]:
    # start strictly inside the VUS band: posterior(points, 0.5) in (0.10, 0.90)
    # i.e. points in [-2, 2]; +/-3 already sits on a reporting threshold.
    cases = []
    for _ in range(n):
        is_path = rng.random() < 0.5
        start = rng.choice([-2, -1, 0, 1, 2])
        cases.append(SimCase(start, is_path, make_actions()))
    return cases


@dataclass
class PolicyStats:
    policy: str
    n: int
    resolved_rate: float
    accuracy: float                  # of resolved, fraction correctly classified
    mean_cost: float
    mean_days: float
    mean_tests: float


def evaluate_policy(cases: list[SimCase], policy: str, seed: int = 0) -> PolicyStats:
    rng = random.Random(seed)
    res = [run_case(c, policy, rng) for c in cases]
    resolved = [r for r in res if r.resolved]
    n = len(res)
    return PolicyStats(
        policy=policy, n=n,
        resolved_rate=len(resolved) / n if n else 0.0,
        accuracy=sum(r.correct for r in resolved) / len(resolved) if resolved else 0.0,
        mean_cost=sum(r.cost_usd for r in res) / n if n else 0.0,
        mean_days=sum(r.days for r in res) / n if n else 0.0,
        mean_tests=sum(r.n_tests for r in res) / n if n else 0.0,
    )


def compare(n_cases: int = 2000, seed: int = 42) -> dict[str, PolicyStats]:
    rng = random.Random(seed)
    cohort = make_cohort(n_cases, rng)
    # same cohort + same per-policy seed so result draws are comparable
    return {pol: evaluate_policy(cohort, pol, seed=seed + 1)
            for pol in ("voi", "voi_stop", "greedy", "fixed", "random")}


if __name__ == "__main__":
    stats = compare()
    print(f"{'policy':8s} {'resolved':>9s} {'accuracy':>9s} {'$cost':>8s} {'days':>7s} {'tests':>6s}")
    for s in stats.values():
        print(f"{s.policy:8s} {s.resolved_rate:9.1%} {s.accuracy:9.1%} "
              f"{s.mean_cost:8.0f} {s.mean_days:7.1f} {s.mean_tests:6.2f}")
