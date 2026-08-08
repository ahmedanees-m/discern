"""Figure 2: discrimination and calibration on the expert-panel surface.

Rendered from the committed benchmark JSON, on the shared journal
style in style.py.
"""
from __future__ import annotations

import json
import os

import matplotlib.pyplot as plt
import numpy as np

from core.stats import fmt
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

    from bench.phase_r_variant import oof_isotonic, shared_folds
    from bench.track1b_erepo_headtohead import TIMESPLIT_AFTER, load_rows

    st.apply()
    m = _load(os.path.join(BENCH, "phase_r_variant_metrics.json"))
    rows = load_rows()

    def surface(rs):
        mis = [r for r in rs if r["is_missense"] and r["label"] is not None]
        return mis, [r["label"] for r in mis]

    fig, axes = plt.subplots(1, 3, figsize=(st.FULL_W, 2.5))

    # Panel A - ROC
    mis, y = surface(rows)
    ax = axes[0]
    auroc_lines = []
    for name, key, c in (("DISCERN", "discern_points", st.BLUE), ("REVEL", "revel", st.ORANGE),
                         ("AlphaMissense", "alphamissense", st.GREEN)):
        pairs = [(r[key], yy) for r, yy in zip(mis, y, strict=False) if r[key] is not None]
        fpr, tpr, _ = roc_curve([b for _, b in pairs], [float(a) for a, _ in pairs])
        d = m["eRepo_primary"]["discrimination"][name]
        ax.plot(fpr, tpr, color=c, lw=1.3, label=name)
        auroc_lines.append((c, f"{name}  {fmt(d['auroc'])} "
                               f"[{fmt(d['auroc_ci95'][0])}-{fmt(d['auroc_ci95'][1])}]"))
    ax.plot([0, 1], [0, 1], ls=":", color=st.GREY, lw=0.8)
    dv = m["eRepo_primary"]["discern_vs_revel"]["delong"]
    ax.set_title(f"eRepo missense (n={m['eRepo_primary']['n_missense_PB']})", loc="left",
                 fontweight="bold", pad=6)
    st.panel_label(ax, "A", dx=-0.30, dy=1.16)
    ax.set_xlabel("false positive rate")
    ax.set_ylabel("true positive rate")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.04)
    # The block starts high enough that the last line clears the x axis: at the previous spacing
    # the DeLong line sat on the spine.
    for i, (c, line) in enumerate(auroc_lines):
        ax.text(0.50, 0.345 - i * 0.082, line, transform=ax.transAxes, fontsize=st.MICRO,
                color=c, ha="left", va="center")
    ax.text(0.50, 0.345 - len(auroc_lines) * 0.082,
            f"DeLong p={dv['p_value']:.2f} (ns)", transform=ax.transAxes,
            fontsize=st.MICRO, style="italic", color=st.INK, ha="left", va="center")

    # Panels B and C - reliability, both surfaces
    for ax, surf, letter, title in (
            (axes[1], rows, "B", "Reliability, eRepo-primary"),
            (axes[2], [r for r in rows if (r["approval_date"] or "") > TIMESPLIT_AFTER],
             "C", "Reliability, time-split")):
        ece_lines = []
        mis, y = surface(surf)
        cal = m["eRepo_primary" if surf is rows else "time_split"]["calibration"]
        ax.plot([0, 1], [0, 1], ls=":", color=st.GREY, lw=0.8)
        for name, key, c, lab in (("DISCERN", "discern_points", st.BLUE, "DISCERN_isotonic_oof"),
                                  ("REVEL", "revel", st.ORANGE, "REVEL_isotonic_oof"),
                                  ("AlphaMissense", "alphamissense", st.GREEN,
                                   "AlphaMissense_isotonic_oof")):
            keep = [i for i, r in enumerate(mis) if r[key] is not None]
            s = np.array([float(mis[i][key]) for i in keep])
            yy = np.array([y[i] for i in keep], int)
            p = oof_isotonic(s, yy, shared_folds(yy))
            bins = np.linspace(0, 1, 11)
            xs, ys = [], []
            for i in range(10):
                sel = (p >= bins[i]) & (p < bins[i + 1] if i < 9 else p <= 1.0)
                if sel.sum() >= 5:
                    xs.append(p[sel].mean())
                    ys.append(yy[sel].mean())
            e = cal[lab]
            ax.plot(xs, ys, "o-", color=c, ms=3, lw=1.1, label=name)
            ece_lines.append((c, f"{fmt(e['ece'])} [{fmt(e['ece_ci95'][0])}-"
                                 f"{fmt(e['ece_ci95'][1])}]"))
        # Which curve is which goes in a short key at the top left; the numbers go alone at the
        # bottom right. Split this way neither block is wide enough to reach the curve, which a
        # combined "name + value + interval" line was.
        ax.legend(loc="upper left", frameon=False, fontsize=st.MICRO, handlelength=1.2,
                  borderpad=0.1, labelspacing=0.25)
        ax.text(0.97, 0.300, "ECE (95% CI)", transform=ax.transAxes, fontsize=st.MICRO,
                ha="right", va="center", fontweight="bold")
        for i, (c, line) in enumerate(ece_lines):
            ax.text(0.97, 0.222 - i * 0.062, line, transform=ax.transAxes, fontsize=st.MICRO,
                    color=c, ha="right", va="center")
        ax.set_title(title, loc="left", fontweight="bold", pad=6)
        st.panel_label(ax, letter, dx=-0.26, dy=1.16)
        ax.set_xlabel("predicted probability")
        ax.set_ylabel("observed frequency")
    fig.tight_layout()
    return st.save(fig, outdir, "fig2_discrimination_and_calibration")


# --------------------------------------------------------------------------------------------
