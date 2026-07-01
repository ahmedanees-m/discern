"""Verifiable learning loop - Bayesian prior updating (plan Step 5.2).

Only the action cost/yield/attainability *priors* update from realized outcomes - the
verdict logic is NEVER opaquely retrained. Each update is a Beta-Bernoulli posterior
on an assay's sensitivity/specificity, attributable to specific cases, so improvement
is auditable.
"""
from __future__ import annotations

from dataclasses import dataclass


def update_sensitivity(prior_alpha: float, prior_beta: float, resolved: bool) -> tuple[float, float]:
    """Beta-Bernoulli update of an assay's sensitivity prior from one realized outcome."""
    return (prior_alpha + (1 if resolved else 0), prior_beta + (0 if resolved else 1))


@dataclass
class BetaPrior:
    alpha: float = 1.0
    beta: float = 1.0

    @property
    def mean(self) -> float:
        return self.alpha / (self.alpha + self.beta)

    @property
    def n(self) -> float:
        return self.alpha + self.beta - 2          # observations beyond the uniform prior

    def update(self, success: bool) -> BetaPrior:
        a, b = update_sensitivity(self.alpha, self.beta, success)
        return BetaPrior(a, b)
