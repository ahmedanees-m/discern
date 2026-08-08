"""The graphical abstract required at submission (920 x 300 px, landscape).

Sized at exactly the requested pixel dimensions (9.2 x 3.0 in at 100 dpi) and kept to a single
idea, since it is displayed small and beneath the abstract. This is the one item that does not
take the 170 mm page width from style.py: the journal specifies its size in pixels.
"""
from __future__ import annotations

import json
import os

import matplotlib.pyplot as plt

from figures.scripts import style as st

# scripts/ -> figures/ -> repository root
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BENCH = os.path.join(ROOT, "bench")
EVAL = os.path.join(ROOT, "eval")

WIDTH_PX, HEIGHT_PX, DPI = 920, 300, 100


def _load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def build(outdir):
    st.apply()
    m = _load(os.path.join(BENCH, "phase_r_variant_metrics.json"))
    ca = m["ceiling_attribution"]
    ic = m["added_value_over_ranking_score"]["intrinsic_only_ceiling"]
    reached = ic["max_discern_points_on_missense"]
    need = ca["lp_threshold_points"]

    # Half a pixel of slack: 920/100 is 9.199999... in binary floating point, and the canvas
    # truncates, so the plain division saves a 919 px image against a spec that says 920.
    fig = plt.figure(figsize=((WIDTH_PX + 0.5) / DPI, (HEIGHT_PX + 0.5) / DPI), dpi=DPI)
    fig.patch.set_facecolor("white")

    ax = fig.add_axes([0.055, 0.36, 0.355, 0.30])
    ax.barh([0], [reached], color=st.BLUE, height=0.62, zorder=2)
    ax.axvline(need, color=st.ORANGE, lw=2.0, zorder=3)
    ax.text(need + 0.2, 0.0, f"{need:.0f} points needed for\nLikely Pathogenic",
            va="center", ha="left", fontsize=8.5, color=st.ORANGE, fontweight="bold",
            linespacing=1.4, zorder=4)
    ax.text(reached - 0.18, 0.0, f"{reached:.0f}", va="center", ha="right", fontsize=11,
            color="white", fontweight="bold", zorder=4)
    ax.set_xlim(0, 9.4)
    ax.set_ylim(-0.55, 0.55)
    ax.set_yticks([])
    ax.set_xticks(range(0, 10, 2))
    ax.tick_params(labelsize=8)
    # Left-aligned so the caption cannot run off the left edge of the canvas, which centring did.
    ax.set_xlabel("ACMG points reachable from sequence evidence alone", fontsize=8.5, loc="left")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)

    # A hairline keeps the claim and its evidence visually separate at display size.
    fig.add_artist(plt.Line2D([0.455, 0.455], [0.14, 0.86], color=st.LIGHT, lw=1.0))

    x = 0.487
    fig.text(x, 0.90, "Intrinsic sequence evidence cannot classify\na missense variant as pathogenic",
             fontsize=12.5, fontweight="bold", va="top", linespacing=1.35)
    fig.text(x, 0.575, f"{ca['reach_lp_as_scored']} of {ca['n']} expert-classified missense variants "
                       f"reach a\npathogenic band under a partitioned ACMG framework.",
             fontsize=9.5, va="top", linespacing=1.4)
    fig.text(x, 0.285, "Restoring every criterion the partition withholds or cannot\n"
                       f"annotate recovers only {ca['reach_lp_with_both']} of "
                       f"{m['eRepo_primary']['n_path']}.",
             fontsize=9.5, va="top", color="#555555", linespacing=1.4)

    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "graphical_abstract.png")
    tif = os.path.join(outdir, "graphical_abstract.tif")
    # The shared style saves every figure on a tight bounding box, which is right for a manuscript
    # figure and wrong here: the journal asks for exactly 920 x 300 px, so the canvas is saved whole.
    with plt.rc_context({"savefig.bbox": "standard", "savefig.pad_inches": 0.0}):
        fig.savefig(out, dpi=DPI, facecolor="white")
        fig.savefig(tif, dpi=DPI, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    st.flatten_tiff(tif, DPI)
    plt.close(fig)
    return [out, tif]


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ROOT)
    print(build(os.path.join(ROOT, "figures", "out")))
