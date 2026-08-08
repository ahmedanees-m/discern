"""Figure S2: sensitivity of the diagnosis arm to the strength of the gene term.

Rendered from the committed benchmark JSON, on the shared journal style in style.py.

The y axis starts at 70, not 0: every series lies between 80 and 100, and a zero-based axis
compressed the whole result into the top fifth of the panel. The break is stated on the axis.
"""
from __future__ import annotations

import json
import os

import matplotlib.lines as mlines
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
    m = _load(os.path.join(BENCH, "phase_r_gene_term_sensitivity.json"))
    grid = sorted(m["grid"], key=lambda r: r["likelihood_ratio"])
    lr = [r["likelihood_ratio"] for r in grid]

    fig, ax = plt.subplots(figsize=(st.FULL_W, 2.85))
    fig.subplots_adjust(left=0.085, right=0.995, top=0.86, bottom=0.175)

    # Bands first, so the series sit on top of them.
    locked = [r for r in grid if not r["deciding_assay_can_overturn_gene"]]
    for r in locked:
        ax.axvspan(r["likelihood_ratio"] * 0.82, r["likelihood_ratio"] * 1.22,
                   color=st.VERMILION, alpha=0.07, lw=0, zorder=0)

    series = [("curated Top-1", "curated_top1", st.BLUE, "o", "-"),
              ("hard-stop sensitivity", "hardstop_sensitivity", st.GREEN, "s", "--"),
              ("hard-stop specificity", "hardstop_specificity", st.PURPLE, "^", ":")]
    for label, key, colour, marker, ls in series:
        ax.plot(lr, [r[key] * 100 for r in grid], marker=marker, ls=ls, color=colour,
                lw=1.2, ms=4.0, label=label, zorder=3)

    committed = next(r for r in grid if r["committed"])
    ax.axvline(committed["likelihood_ratio"], color=st.VERMILION, lw=1.0, ls="--", zorder=2)
    # Left of the line: the lower right belongs to the legend, and the region under the rising
    # Top-1 curve is the only empty space wide enough for three lines.
    ax.text(committed["likelihood_ratio"] * 0.92, 76.5,
            "committed value: the largest\nat which the deciding assay\nstill overturns the gene",
            color=st.VERMILION, fontsize=st.TINY, va="center", ha="right", linespacing=1.5,
            zorder=4)

    ax.set_xscale("log")
    ax.set_xticks(lr)
    ax.set_xticklabels([str(v) for v in lr])
    ax.minorticks_off()
    ax.set_xlabel("P(G|D) likelihood ratio")
    ax.set_ylabel("percent (axis starts at 70)")
    ax.set_ylim(70, 104)
    ax.set_yticks([70, 80, 90, 100])
    ax.grid(axis="y", color=st.LIGHT, lw=0.5, zorder=0)
    ax.set_axisbelow(True)

    handles = [mlines.Line2D([], [], color=c, marker=mk, ls=ls, lw=1.2, ms=4.0, label=lab)
               for lab, _k, c, mk, ls in series]
    if locked:
        handles.append(mpatches.Patch(facecolor=st.VERMILION, alpha=0.07,
                                      label="gene term too strong to overturn"))
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=st.TINY,
              ncol=2, handlelength=2.0, columnspacing=1.4, borderpad=0.1)
    ax.set_title("Gene-term sensitivity: safety invariant, diagnosis flat above 4",
                 loc="left", fontweight="bold", pad=5)
    return st.save(fig, outdir, "figS1_gene_term_sensitivity")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ROOT)
    print(build(os.path.join(ROOT, "figures", "out")))
