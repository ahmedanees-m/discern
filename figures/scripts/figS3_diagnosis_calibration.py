"""Figure S4: diagnosis-posterior calibration on the curated benchmark.

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


def _load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def build(outdir):
    st.apply()
    m = _load(os.path.join(BENCH, "track3_metrics.json"))["calibration_diagnosis"]
    fig, ax = plt.subplots(figsize=(st.HALF_W * 1.30, st.HALF_W * 1.18))
    fig.subplots_adjust(left=0.165, right=0.985, top=0.865, bottom=0.145)
    ax.plot([0, 1], [0, 1], ls=":", color=st.GREY, lw=0.8)
    ax.plot([m["mean_confidence"]], [m["accuracy"]], "o", color=st.BLUE, ms=8)
    ax.annotate(f"mean confidence {m['mean_confidence']:.3f}\naccuracy {m['accuracy']:.3f}\n"
                f"ECE {m['ece']:.3f} (n={m['n']})",
                xy=(m["mean_confidence"], m["accuracy"]), xytext=(0.33, 0.66), fontsize=st.TINY,
                arrowprops=dict(arrowstyle="->", color=st.GREY, lw=0.8))
    ax.axhspan(0.8, 1.0, color=st.GREEN, alpha=0.07)
    ax.text(0.035, 0.965, "zero confidently-wrong calls above 0.8\n"
            f"({m['confidently_wrong_at_0.8']} of {m['n']})", fontsize=st.TINY, va="top",
            color=st.GREEN, linespacing=1.4)
    ax.set_xlabel("mean predicted confidence")
    ax.set_ylabel("observed accuracy")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.05)
    ax.set_title(f"Diagnosis-posterior calibration\n(proof of concept, n={m['n']})",
                 loc="left", fontweight="bold", pad=5)
    return st.save(fig, outdir, "figS3_diagnosis_calibration")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ROOT)
    print(build(os.path.join(ROOT, "figures", "out")))
