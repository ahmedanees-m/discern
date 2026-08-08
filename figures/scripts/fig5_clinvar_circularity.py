"""Figure 5: ClinVar circularity in a live-database classifier.

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
    from sklearn.metrics import roc_curve

    from bench.track1b_erepo_headtohead import TIMESPLIT_AFTER, load_rows
    st.apply()
    m = _load(os.path.join(BENCH, "track1b_erepo_metrics.json"))
    c1 = _load(os.path.join(BENCH, "track1_metrics.json"))["clinvar_circularity"]
    rows = load_rows()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(st.FULL_W, 2.9))
    for rs, style, tag in ((rows, "-", "eRepo-primary"),
                           ([r for r in rows if (r["approval_date"] or "") > TIMESPLIT_AFTER],
                            "--", "time-split")):
        mis = [r for r in rs if r["is_missense"] and r["label"] is not None]
        y = [r["label"] for r in mis]
        for name, key, c in (("GeneBe", "genebe_score", st.VERMILION), ("DISCERN", "discern_points", st.BLUE)):
            pairs = [(r[key], yy) for r, yy in zip(mis, y, strict=False) if r[key] is not None]
            fpr, tpr, _ = roc_curve([b for _, b in pairs], [float(a) for a, _ in pairs])
            surf = "eRepo_primary" if tag == "eRepo-primary" else "time_split"
            key2 = "GeneBe_exhibit" if name == "GeneBe" else "DISCERN"
            auc = m[surf]["discrimination"][key2]["auroc"]
            ax1.plot(fpr, tpr, style, color=c, lw=1.3, label=f"{name} {tag}: {auc:.3f}")
    ax1.plot([0, 1], [0, 1], ls=":", color=st.GREY, lw=0.8)
    ax1.set_xlabel("false positive rate")
    ax1.set_ylabel("true positive rate")
    ax1.set_title("A live-ClinVar tool is perfectly separable", loc="left", fontweight="bold")
    st.panel_label(ax1, "A", dx=-0.115, dy=1.07)
    ax1.legend(loc="lower right", frameon=False)

    match, total = c1["genebe_class_matches_clinvar_direction"].split("/")
    ax2.bar(["GeneBe\n(consumes ClinVar)", "DISCERN / REVEL /\nAlphaMissense"],
            [100 * int(match) / int(total), 0], color=[st.VERMILION, st.BLUE], alpha=0.85)
    ax2.text(0, 100 * int(match) / int(total) + 2, f"{match}/{total} calls match\nthe ClinVar direction",
             ha="center", fontsize=6.5)
    ax2.text(1, 4, "apply no\nClinVar-derived code", ha="center", fontsize=6.5)
    ax2.set_ylabel("% of resolved calls matching ClinVar")
    ax2.set_ylim(0, 118)
    ax2.set_title("Reproduction, not prediction", loc="left", fontweight="bold")
    st.panel_label(ax2, "B", dx=-0.135, dy=1.07)
    fig.tight_layout()
    return st.save(fig, outdir, "fig5_clinvar_circularity")


# --------------------------------------------------------------------------------------------
