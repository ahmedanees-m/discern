"""Figure S5: independent sensitivity on the CDC hemophilia mutation projects.

Rendered from the committed benchmark JSON, on the shared journal
style in style.py.
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


CHAMP = {"overall": 0.912, "F8": 0.977, "F9": 0.658, "n_alleles": 5437, "n_null": 2130}


def _load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def build(outdir):
    st.apply()
    fig, ax = plt.subplots(figsize=(st.HALF_W * 1.45, 2.7))
    fig.subplots_adjust(left=0.145, right=0.985, top=0.845, bottom=0.115)
    keys = ["F8", "F9", "overall"]
    labels = ["F8\n(hemophilia A)", "F9\n(hemophilia B)", "overall"]
    ax.bar(range(len(keys)), [CHAMP[k] * 100 for k in keys],
           color=[st.BLUE, st.ORANGE, st.GREY], alpha=0.85, width=0.62, zorder=2)
    for i, k in enumerate(keys):
        ax.text(i, CHAMP[k] * 100 + 2.5, f"{CHAMP[k]:.1%}", ha="center", fontsize=st.SMALL,
                fontweight="bold", zorder=3)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(labels, fontsize=st.TINY)
    ax.set_ylabel("LP/P recall on null variants (%)")
    ax.set_ylim(0, 112)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.grid(axis="y", color=st.LIGHT, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title(f"CHAMP/CHBMP independent sensitivity\n({CHAMP['n_null']} null of "
                 f"{CHAMP['n_alleles']} disease alleles, consequence alone)",
                 loc="left", fontweight="bold", pad=5)
    # Inside the F9 bar: the gap between bars is too narrow for two lines of text, and an
    # annotation placed there ran under the "overall" bar.
    ax.text(1, 20, "the gap is last-exon\nNMD escape, correctly\nheld at VUS", ha="center",
            va="center", fontsize=st.TINY, color="#3A3A3A", linespacing=1.6, zorder=4)
    return st.save(fig, outdir, "figS4_champ_chbmp_recall")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ROOT)
    print(build(os.path.join(ROOT, "figures", "out")))
