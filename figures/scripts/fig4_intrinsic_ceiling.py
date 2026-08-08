"""Figure 4: the ceiling intrinsic evidence imposes on missense classification.

Rendered from the committed benchmark JSON, on the shared journal
style in style.py.
"""
from __future__ import annotations

import json
import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from figures.scripts import style as st

# scripts/ -> figures/ -> repository root
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BENCH = os.path.join(ROOT, "bench")
EVAL = os.path.join(ROOT, "eval")


def _load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def build(outdir):
    st.apply()
    """The lead argument."""
    from bench.track1b_erepo_headtohead import load_rows
    m = _load(os.path.join(BENCH, "phase_r_variant_metrics.json"))
    ca = m["ceiling_attribution"]
    ic = m["added_value_over_ranking_score"]["intrinsic_only_ceiling"]
    mis = [r for r in load_rows() if r["is_missense"] and r["label"] is not None]
    pts = [r["discern_points"] for r in mis]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(st.FULL_W, 3.0))

    lo, hi = int(min(pts)) - 1, 8
    ax1.hist(pts, bins=np.arange(lo, hi + 1) - 0.5, color=st.BLUE, alpha=0.85, edgecolor="white")
    ax1.axvline(6, color=st.VERMILION, lw=1.4, ls="--")
    ax1.text(6.15, ax1.get_ylim()[1] * 0.92, "Likely Pathogenic\nthreshold (6 points)",
             color=st.VERMILION, fontsize=6.5, va="top")
    ax1.axvspan(6, hi, color=st.VERMILION, alpha=0.06)
    ax1.text((6 + hi) / 2, ax1.get_ylim()[1] * 0.45, f"zero of {ca['n']}", color=st.VERMILION,
             fontsize=8, ha="center", fontweight="bold")
    ax1.annotate("observed\nmaximum {:.0f}".format(ic["max_discern_points_on_missense"]),
                 xy=(3, 10), xytext=(4.3, ax1.get_ylim()[1] * 0.28), fontsize=6.5, ha="center",
                 arrowprops=dict(arrowstyle="->", color=st.GREY, lw=0.8,
                                 connectionstyle="arc3,rad=0.25"))
    ax1.set_xlabel("total ACMG points, intrinsic evidence only")
    ax1.set_ylabel("missense variants")
    ax1.set_title("Nothing reaches a pathogenic band", loc="left", fontweight="bold")
    st.panel_label(ax1, "A", dx=-0.135, dy=1.07)

    steps = [("as\nscored", ca["reach_lp_as_scored"], st.GREY, None),
             ("+ no-input\nintrinsic", ca["reach_lp_if_intrinsic_codes_were_available"], st.ORANGE,
              "+ intrinsic criteria with no input here (PM1, PM5, PS1, PS4)"),
             ("+ routed by\npartition", ca["reach_lp_if_routed_codes_were_re_added"], st.GREEN,
              "+ criteria the partition routes away (PS3, PP4, PP1, PM3)"),
             ("both\nrestored", ca["reach_lp_with_both"], st.BLUE, None)]
    ax2.bar([s[0] for s in steps], [s[1] for s in steps], color=[s[2] for s in steps], alpha=0.85)
    for i, (_, v, _c, _d) in enumerate(steps):
        ax2.text(i, v + 7, str(v), ha="center", fontsize=7.5, fontweight="bold")
    ax2.legend(handles=[mpatches.Patch(color=s[2], alpha=0.85, label=s[3])
                        for s in steps if s[3]],
               loc="upper center", frameon=False, fontsize=5.5, bbox_to_anchor=(0.5, 0.78))
    n_path = m["eRepo_primary"]["n_path"]
    ax2.axhline(n_path, color=st.VERMILION, ls="--", lw=1.2)
    ax2.text(-0.45, n_path + 4, f"{n_path} truly pathogenic", color=st.VERMILION, fontsize=6.5,
             ha="left", va="bottom")
    ax2.set_ylabel("variants reaching Likely Pathogenic")
    ax2.set_ylim(0, n_path * 1.15)
    ax2.set_title("Neither cause explains the ceiling", loc="left", fontweight="bold")
    st.panel_label(ax2, "B", dx=-0.115, dy=1.07)
    ax2.tick_params(axis="x", labelsize=6.5)
    fig.tight_layout()
    return st.save(fig, outdir, "fig4_intrinsic_ceiling")


# --------------------------------------------------------------------------------------------
