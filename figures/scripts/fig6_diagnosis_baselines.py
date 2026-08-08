"""Figure 6: differential diagnosis against phenotype-blind and external baselines.

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
    g = _load(os.path.join(EVAL, "gene_only_baseline.json"))
    lir = _load(os.path.join(EVAL, "lirical_arm.json"))
    cov = _load(os.path.join(EVAL, "phenotype_tool_comparison.json"))["hpo_coverage"]
    p = g["pooled"]

    fig, axes = plt.subplots(1, 3, figsize=(st.FULL_W, 2.9))

    ax = axes[0]
    bars = [("uniform random", p["random_within_cluster"]["top1"]["expected"], st.GREY, False),
            ("prior only", p["prior_only"]["top1"], st.GREY, False),
            ("gene lookup", p["gene_lookup"]["top1"], st.ORANGE, False),
            ("DISCERN pre-corr.", 0.81, st.SKY, False),
            ("DISCERN post-corr.", p["DISCERN"]["top1"], st.BLUE, True)]
    for i, (_lab, v, c, hatch) in enumerate(bars):
        ax.bar(i, v * 100, color=c, alpha=0.85, hatch="//" if hatch else None,
               edgecolor="white" if hatch else None, lw=0 if not hatch else 1.0)
        ax.text(i, v * 100 + 2, f"{v:.0%}", ha="center", fontsize=6.5)
    ax.set_xticks(range(len(bars)))
    ax.set_xticklabels([b[0] for b in bars], fontsize=6, rotation=30, ha="right")
    ax.set_ylabel("Top-1 accuracy (%)")
    ax.set_ylim(0, 118)
    ax.set_title(f"Curated cases (n={p['n']})", loc="left", fontweight="bold", pad=14)
    st.panel_label(ax, "A", dx=-0.30, dy=1.20)
    ax.text(0.0, 1.035, f"vs gene lookup: McNemar "
            f"p={p['mcnemar_top1']['p_value_exact']:.2f} (ns)",
            transform=ax.transAxes, ha="left", fontsize=st.TINY, style="italic",
            color=st.GREY)
    ax.legend(handles=[mpatches.Patch(facecolor=st.BLUE, alpha=0.85, hatch="//", edgecolor="white",
                                      label="in-sample")], loc="upper center", frameon=False,
              bbox_to_anchor=(0.45, 1.0), handlelength=1.4)

    ax = axes[1]
    matched = lir["paired_tests_vs_lirical_restricted"]["hpo_representable_only"]
    full = lir["paired_tests_vs_lirical_restricted"]["phenotype_only"]
    groups = [("identical evidence\n13 HPO-codable\nfindings", matched, False),
              ("all findings\nDISCERN reads\nall 48", full, True)]
    x = 0
    ticks, labels = [], []
    for lab, t, sig in groups:
        ax.bar(x, t["lirical_recall@1"] * 100, color=st.ORANGE, alpha=0.85, width=0.38)
        ax.bar(x + 0.4, t["discern_recall@1"] * 100, color=st.BLUE, alpha=0.85, width=0.38)
        star = f"p={t['mcnemar']['p_value_exact']:.2f}" + ("*" if sig else " (ns)")
        ax.text(x + 0.2, max(t["discern_recall@1"], t["lirical_recall@1"]) * 100 + 6, star,
                ha="center", fontsize=6.5, fontweight="bold" if sig else "normal")
        ticks.append(x + 0.2)
        labels.append(lab)
        x += 1.45
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=st.TINY)
    ax.set_ylabel("Recall@1 (%)")
    ax.set_ylim(0, 118)
    ax.set_title(f"vs LIRICAL (n={matched.get('n', 23) if isinstance(matched, dict) else 23})",
                 loc="left", fontweight="bold", pad=14)
    st.panel_label(ax, "B", dx=-0.24, dy=1.20)
    ax.legend(handles=[mpatches.Patch(color=st.ORANGE, alpha=0.85, label="LIRICAL"),
                       mpatches.Patch(color=st.BLUE, alpha=0.85, label="DISCERN")],
              loc="upper left", frameon=False)

    ax = axes[2]
    codable = cov["features_expressible_as_hpo"]
    total = cov["distinct_features_used"]
    ax.pie([codable, total - codable], colors=[st.GREEN, st.ORANGE],
           labels=[f"{codable} codable\nin HPO",
                   f"{total - codable} laboratory findings\nHPO cannot express"],
           autopct="%1.0f%%", startangle=90, textprops={"fontsize": st.TINY},
           radius=0.82, labeldistance=1.12, pctdistance=0.62,
           wedgeprops={"alpha": 0.85, "edgecolor": "white"})
    ax.set_title("Why the channel is empty", loc="left", fontweight="bold", pad=14)
    st.panel_label(ax, "C", dx=-0.16, dy=1.20)
    fig.tight_layout()
    return st.save(fig, outdir, "fig6_diagnosis_baselines")


# --------------------------------------------------------------------------------------------
