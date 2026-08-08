"""Figure 3: per-criterion agreement with expert-panel applications.

Rendered from the committed benchmark JSON, on the shared journal
style in style.py.
"""
from __future__ import annotations

import json
import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

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
    """Agreement per criterion, ordered by kappa, with the three zero-application causes distinct.

    Panel A reads as a gradient from criteria a sequence-only pipeline can derive to those it
    cannot, which is the paper's thesis rather than a wall of low values. Panel B shows the applied
    counts on both sides, so a reader sees the direction and size of each disagreement without
    leaving the figure.
    """
    m = _load(os.path.join(BENCH, "track1b_erepo_metrics.json"))["per_code_kappa_vs_erepo"]
    routed = {"PS3", "BS3", "PP4", "PP1", "BS4", "PM3", "BP2", "PS2", "PM6"}
    not_implemented = {"BP7", "BS2"}

    def cause_colour(code):
        if code in routed:
            return st.GREEN
        return st.GREY if code in not_implemented else st.ORANGE

    scored = sorted(((c, v["kappa"]) for c, v in m.items() if v["kappa"] is not None),
                    key=lambda kv: kv[1])
    zero = sorted(((c, v) for c, v in m.items() if v["discern_applied"] == 0),
                  key=lambda kv: kv[1]["erepo_applied"])

    fig, axes = plt.subplots(1, 2, figsize=(st.FULL_W, 3.9), gridspec_kw={"width_ratios": [1, 1]})
    ax1, ax2 = axes

    # barh puts index 0 at the bottom, so the never-applied block is listed first to place it under
    # the kappa gradient. Read top to bottom the panel then runs from strongest agreement to none.
    labels = [c for c, _ in zero] + [""] + [c for c, _ in scored]
    values = [0] * len(zero) + [0] + [k for _, k in scored]
    colours = [cause_colour(c) for c, _ in zero] + ["none"] + [st.BLUE] * len(scored)
    ax1.barh(range(len(labels)), values, color=colours, alpha=0.9)
    ax1.set_yticks(range(len(labels)))
    ax1.set_yticklabels(labels, fontsize=7)
    for i, (_, k) in enumerate(scored):
        ax1.text(k + 0.02, len(zero) + 1 + i, f"{k:.2f}", va="center", fontsize=6.5)
    for j, (c, _) in enumerate(zero):
        ax1.text(0.02, j, "never applied", va="center", fontsize=6.5,
                 color=cause_colour(c), style="italic")
    ax1.set_xlim(0, 1.12)
    ax1.set_xlabel("Cohen's kappa vs expert applications", fontsize=8)
    ax1.set_title("Criterion-level agreement", loc="left", fontweight="bold", fontsize=st.TITLE)
    st.panel_label(ax1, "A", dx=-0.115, dy=1.045)
    ax1.tick_params(axis="x", labelsize=7)

    y = range(len(labels))
    ax2.barh([i + 0.2 for i in y],
             [m[c]["erepo_applied"] if c else 0 for c in labels], height=0.38,
             color="#555555", alpha=0.85, label="expert panels")
    ax2.barh([i - 0.2 for i in y],
             [m[c]["discern_applied"] if c else 0 for c in labels], height=0.38,
             color=st.BLUE, alpha=0.85, label="DISCERN")
    ax2.set_yticks(list(y))
    ax2.set_yticklabels([""] * len(labels))
    ax2.set_xlabel("variants on which the criterion was applied", fontsize=8)
    ax2.set_title("How often each side applied it", loc="left", fontweight="bold",
                  fontsize=st.TITLE)
    st.panel_label(ax2, "B", dx=-0.055, dy=1.045)
    ax2.legend(loc="upper right", frameon=False, fontsize=7)
    ax2.set_xlim(0, max(v["erepo_applied"] for v in m.values()) * 1.18)
    ax2.tick_params(axis="x", labelsize=7)

    fig.legend(handles=[
        mpatches.Patch(color=st.GREEN, alpha=0.9, label="routed to another factor by the partition"),
        mpatches.Patch(color=st.ORANGE, alpha=0.9, label="variant-intrinsic; no input available here"),
        mpatches.Patch(color=st.GREY, alpha=0.9, label="not implemented (outside missense scope)")],
        loc="lower center", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    return st.save(fig, outdir, "fig3_per_criterion_kappa")


# --------------------------------------------------------------------------------------------
