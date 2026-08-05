"""Points -> posterior probability of pathogenicity.

Tavtigian et al. 2018 Bayesian bridge. Needed by the value-of-information core,
which reasons over probabilities, not raw points.

Checks (prior_p = 0.10):
    C = 10 -> ~0.99   (Pathogenic threshold)
    C =  6 -> ~0.90   (Likely pathogenic threshold)
    C =  0 -> 0.10    (prior)
    C = -7 -> <0.001  (Benign)
"""
from __future__ import annotations


def posterior(points: float, prior_p: float = 0.10) -> float:
    """Posterior probability of pathogenicity given combined ACMG points.

    OddsPath = 350 ** (points / 8); the 350 and the /8 come from calibrating
    the four strength levels (PP/PM/PS/PVS = 1/2/4/8 points) to the ACMG
    odds ladder (Tavtigian 2018).
    """
    odds_path = 350.0 ** (points / 8.0)
    return (odds_path * prior_p) / ((odds_path - 1) * prior_p + 1)
