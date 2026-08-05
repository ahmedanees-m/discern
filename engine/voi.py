"""Value-of-information scoring - the decision core.

From the current points -> posterior `p`. For an action whose positive result adds
`+k_pos` points and negative adds `-k_neg`, with sensitivity/specificity, we compute:
  * the predictive probability of a positive result `q`,
  * expected information gain `EIG = H(p) - E[H(p')]` (Bernoulli entropy), and
  * decision utility `E[deltaU]` = expected gain in probability of reaching an
    *actionable* (reportable) classification.
Actions are ranked by value density = objective / cost, cost = alpha-$ + beta-weeks + gamma-burden.
Myopic by default; shallow lookahead (depth 2-3) when top actions tie.

The math is calibrated to the literature: functional-assay strengths derive from
OddsPath calibration (Brnich 2020); the posterior bridge is Tavtigian 2018.
"""
from __future__ import annotations

import math

from core.schemas import Action
from rules.posterior import posterior


def _H(p: float) -> float:
    """Binary (Bernoulli) entropy in bits."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log2(p) + (1 - p) * math.log2(1 - p))


def _U(p: float, hi: float = 0.90, lo: float = 0.10) -> float:
    """Actionability utility: 1.0 once a reporting threshold is reached, else ramps."""
    if p >= hi or p <= lo:
        return 1.0
    return max((p - 0.5) / (hi - 0.5), (0.5 - p) / (0.5 - lo), 0.0)


_CODE_K = {"PVS": 8.0, "PS": 4.0, "PM": 2.0, "PP": 1.0,
           "BVS": 8.0, "BS": 4.0, "BM": 2.0, "BP": 1.0}


def action_k(a: Action) -> tuple[float, float]:
    """(k_pos, k_neg) an action yields, split by direction.

    A positive result adds the strongest *pathogenic* code's points; a negative result
    subtracts the strongest *benign* code's points. These can differ - e.g. segregation
    yields PP1 (+1) on co-segregation but BS4 (-4) on non-segregation. Taking a single
    max over both would overstate the weaker direction.
    """
    pos, neg = [], []
    for c in a.yields_codes:
        prefix = c[:3] if c[:3] in _CODE_K else c[:2]
        k = _CODE_K.get(prefix, 1.0)
        (pos if prefix.startswith("P") else neg).append(k)
    k_pos = max(pos) if pos else (max(neg) if neg else 1.0)
    k_neg = max(neg) if neg else (max(pos) if pos else 1.0)
    return k_pos, k_neg


def score_action(points: float, a: Action, prior_p: float = 0.10,
                 k_pos: float | None = None, k_neg: float | None = None) -> dict:
    """Score one action at the current points. k defaults to the action's yielded strength."""
    if k_pos is None or k_neg is None:
        ak_pos, ak_neg = action_k(a)
        k_pos = ak_pos if k_pos is None else k_pos
        k_neg = ak_neg if k_neg is None else k_neg
    p = posterior(points, prior_p)
    q = a.sensitivity * p + (1 - a.specificity) * (1 - p)        # P(positive result)
    p_pos = posterior(points + k_pos, prior_p)
    p_neg = posterior(points - k_neg, prior_p)
    eig = _H(p) - (q * _H(p_pos) + (1 - q) * _H(p_neg))
    dU = q * _U(p_pos) + (1 - q) * _U(p_neg) - _U(p)
    return {"eig": eig, "delta_utility": dU, "q": q,
            "p": p, "p_pos": p_pos, "p_neg": p_neg}


def value_density(points: float, a: Action, alpha: float = 1.0, beta: float = 50.0,
                  gamma: float = 0.0, objective: str = "delta_utility", **kw) -> tuple[float, dict]:
    """Return (value-per-cost, score-detail). Cost scalarizes $ + time + burden.

    Non-improving actions (objective <= 0) return the raw (negative) objective and are
    NOT divided by cost - otherwise a larger cost would make a useless action look
    'least bad' and get top-ranked. This keeps every positive-utility action ranked
    above every non-improving one, and ranks non-improving ones by raw utility.
    """
    s = score_action(points, a, **kw)
    val = s[objective]
    cost = alpha * a.cost_usd + beta * a.turnaround_days / 7.0 + gamma * a.burden
    if val <= 0 or cost <= 0:
        return val, s
    return val / cost, s


def rank_actions(points: float, actions: list[Action], **kw) -> list[tuple[Action, float, dict]]:
    """Rank actions by value density (desc). kw -> value_density (alpha/beta/objective/...)."""
    scored = [(a, *value_density(points, a, **kw)) for a in actions]
    return sorted(scored, key=lambda t: t[1], reverse=True)


def lookahead(points: float, actions: list[Action], depth: int = 2,
              prior_p: float = 0.10, k: float = 4.0) -> tuple[float, list[Action]]:
    """Best expected actionability over short action-sequences (exhaustive, shallow).

    Handles the two cases myopic VOI misses: a cheap test that *enables* a better
    one, and two cheap tests that *jointly* cross threshold.
    """
    if depth == 0 or not actions:
        return _U(posterior(points, prior_p)), []
    best_val, best_seq = -math.inf, []
    for a in actions:
        s = score_action(points, a, prior_p=prior_p)
        rest = [x for x in actions if x is not a]
        v_pos, seq_pos = lookahead(points + k, rest, depth - 1, prior_p, k)
        v_neg, _ = lookahead(points - k, rest, depth - 1, prior_p, k)
        val = s["q"] * v_pos + (1 - s["q"]) * v_neg
        if val > best_val:
            best_val, best_seq = val, [a, *seq_pos]
    return best_val, best_seq


def outcome_variance(points: float, a: Action, prior_p: float = 0.10,
                     k_pos: float = 4.0, k_neg: float = 4.0) -> float:
    """Variance of the post-action actionability utility (for risk aversion, Section 3.6)."""
    s = score_action(points, a, prior_p, k_pos, k_neg)
    u_pos, u_neg = _U(s["p_pos"]), _U(s["p_neg"])
    mean = s["q"] * u_pos + (1 - s["q"]) * u_neg
    return s["q"] * (u_pos - mean) ** 2 + (1 - s["q"]) * (u_neg - mean) ** 2


def risk_gate(ranked: list[tuple], tau: float = 0.0) -> list[tuple]:
    """Drop actions whose expected deltaU is below tau (clinician-tunable, Section 3.6).

    ``ranked`` items are (action, value_density, score_detail).
    """
    return [t for t in ranked if t[2].get("delta_utility", 0.0) >= tau]


def pareto_frontier(points: float, actions: list[Action], prior_p: float = 0.10) -> list[dict]:
    """Cost-EIG Pareto frontier (Section 3.5): cheapest action for each achievable EIG level.

    Lets the clinician choose 'cheapest path to LP' vs 'fastest' rather than committing
    to a single scalarization.
    """
    scored = [(a, a.cost_usd, score_action(points, a, prior_p=prior_p)["eig"]) for a in actions]
    scored.sort(key=lambda t: (t[1], -t[2]))      # cost asc, eig desc
    frontier, best_eig = [], -math.inf
    for a, cost, eig in scored:
        if eig > best_eig + 1e-12:
            frontier.append({"action": a.name, "cost_usd": cost, "turnaround_days": a.turnaround_days,
                             "eig_bits": round(eig, 4)})
            best_eig = eig
    return frontier
