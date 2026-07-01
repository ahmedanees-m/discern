"""Programmatic manuscript figures (plan Phase 8) - 300 dpi, reproducible.

Generates the headline figures from the deterministic engine + simulator so the paper's
figures are regenerable from code. Run: ``python -m figures.generate_all``.
Requires matplotlib (in the API/dev env); writes PNGs into figures/out/.
"""
from __future__ import annotations

import os


def fig_policy_comparison(outdir: str) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from sim.odyssey_sim import compare
    stats = compare(n_cases=2000, seed=42)
    policies = list(stats)
    costs = [stats[p].mean_cost for p in policies]
    resolved = [stats[p].resolved_rate * 100 for p in policies]

    fig, ax1 = plt.subplots(figsize=(7, 4))
    ax1.bar(policies, costs, color="#4C72B0", alpha=0.8)
    ax1.set_ylabel("mean $ cost / case")
    ax2 = ax1.twinx()
    ax2.plot(policies, resolved, "o-", color="#C44E52")
    ax2.set_ylabel("resolved (%)")
    ax1.set_title("Policy comparison (odyssey simulator)")
    path = os.path.join(outdir, "fig_policy_comparison.png")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def fig_posterior_bridge(outdir: str) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from rules.posterior import posterior
    pts = list(range(-10, 13))
    ys = [posterior(p) for p in pts]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(pts, ys, "-o", color="#55A868")
    for thr in (10, 6, -6, -7):
        ax.axvline(thr, ls="--", color="grey", lw=0.8)
    ax.set_xlabel("combined ACMG points")
    ax.set_ylabel("posterior P(pathogenic)")
    ax.set_title("Tavtigian points -> posterior bridge")
    path = os.path.join(outdir, "fig_posterior_bridge.png")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def main() -> None:
    outdir = os.path.join(os.path.dirname(__file__), "out")
    os.makedirs(outdir, exist_ok=True)
    for fn in (fig_policy_comparison, fig_posterior_bridge):
        print("wrote", fn(outdir))


if __name__ == "__main__":
    main()
