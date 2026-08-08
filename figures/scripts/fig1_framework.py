"""Figure 1: the framework and its evidence partition.

Panel A is a schematic; panel B is generated from the live partition map, so it cannot drift from
the code. A sits above B and spans the full width, because the three evidence streams, the joint
posterior and the four outputs need horizontal room; side by side, the output labels collided with
panel B's category axis.
"""
from __future__ import annotations

import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from figures.scripts import style as st

# scripts/ -> figures/ -> repository root
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _box(ax, x, y, w, h, text, face, edge, fontsize=st.SMALL, weight="normal"):
    ax.add_patch(mpatches.FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.012",
        facecolor=face, edgecolor=edge, linewidth=0.9, alpha=0.30 if face != "white" else 1.0,
        zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize,
            fontweight=weight, zorder=3, linespacing=1.35)


def build(outdir):
    from rules.vcep.partition import FACTOR_OF

    st.apply()
    fig = plt.figure(figsize=(st.FULL_W, 5.1))
    # A on top, B beneath: two rows, A given the height its boxes need.
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.22,
                          left=0.085, right=0.985, top=0.945, bottom=0.075)

    # ---------------------------------------------------------------- Panel A
    ax = fig.add_subplot(gs[0])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Coupled disease-variant model", loc="left", fontweight="bold", pad=6)
    st.panel_label(ax, "A", dx=-0.055, dy=1.14)

    inputs = [("Phenotype\nHPO terms and\npertinent negatives", 0.70, st.SKY),
              ("Variant-intrinsic genetics\nfrequency, in-silico,\nnull consequence", 0.385, st.BLUE),
              ("Laboratory and\nfunctional results", 0.07, st.GREEN)]
    bw, bh = 0.235, 0.235
    for label, y, colour in inputs:
        _box(ax, 0.0, y, bw, bh, label, colour, colour)
        ax.annotate("", xy=(0.365, 0.505), xytext=(bw + 0.004, y + bh / 2),
                    arrowprops=dict(arrowstyle="-|>", color=st.GREY, lw=0.9,
                                    shrinkA=0, shrinkB=2,
                                    connectionstyle="arc3,rad=0.10"), zorder=1)

    _box(ax, 0.365, 0.335, 0.215, 0.34, "P(D, V | E)\n\njoint posterior\nover disease\nand variant",
         st.ORANGE, st.ORANGE, fontsize=st.SMALL, weight="bold")

    outs = ["Ranked differential\ndiagnosis", "Variant classification\nwith criterion trail",
            "Treatment-safety\nhard stop", "Next observation by\nexpected information gain"]
    ow, oh = 0.275, 0.185
    for i, label in enumerate(outs):
        y = 0.775 - i * 0.245
        _box(ax, 0.695, y, ow, oh, label, "white", st.GREY, fontsize=st.TINY)
        ax.annotate("", xy=(0.692, y + oh / 2), xytext=(0.584, 0.505),
                    arrowprops=dict(arrowstyle="-|>", color=st.GREY, lw=0.9,
                                    shrinkA=2, shrinkB=0,
                                    connectionstyle="arc3,rad=0.10"), zorder=1)

    # ---------------------------------------------------------------- Panel B
    ax2 = fig.add_subplot(gs[1])
    ax2.set_title("Evidence partition: every criterion enters exactly one factor",
                  loc="left", fontweight="bold", pad=6)
    st.panel_label(ax2, "B", dx=-0.055, dy=1.14)

    colours = {"variant_intrinsic": st.BLUE, "functional": st.GREEN, "disease_pp4": st.ORANGE,
               "segregation": st.PURPLE, "phasing": st.VERMILION, "denovo": st.GREY}
    pretty = {"variant_intrinsic": "variant-intrinsic", "functional": "functional",
              "disease_pp4": "disease (PP4)", "segregation": "segregation",
              "phasing": "phasing", "denovo": "de novo"}
    factors = {}
    for code, fac in FACTOR_OF.items():
        factors.setdefault(fac, []).append(code)

    order = ["variant_intrinsic", "functional", "disease_pp4", "segregation", "phasing", "denovo"]
    order = [f for f in order if f in factors]
    ys = list(range(len(order)))[::-1]
    for y, fac in zip(ys, order, strict=True):
        codes = sorted(factors[fac])
        ax2.barh(y, len(codes), color=colours[fac], alpha=0.85, height=0.62, zorder=2)
        # The variant-intrinsic factor owns 19 criteria. On one line they ran off the axis and
        # inside the bar they overflowed it vertically, so every list sits to the right of its bar,
        # wrapped, with the x limit extended to leave room.
        ax2.text(len(codes) + 0.5, y, st.wrap(", ".join(codes), 38), va="center", ha="left",
                 fontsize=st.TINY, linespacing=1.4, zorder=3)
    ax2.set_yticks(ys)
    ax2.set_yticklabels([pretty[f] for f in order])
    ax2.set_xlabel("number of ACMG/AMP criteria owned")
    ax2.set_xlim(0, 34)
    ax2.set_ylim(-0.7, len(order) - 0.3)

    return st.save(fig, outdir, "fig1_architecture_and_partition")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ROOT)
    print(build(os.path.join(ROOT, "figures", "out")))
