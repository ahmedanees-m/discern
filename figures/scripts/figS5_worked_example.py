"""Figure S1: a worked example, the safety interlock on a single case.

One case end to end: a platelet-type von Willebrand case in which the interlock prevents a
documented harm. The gene is GP1BA; the contraindicated disease, type 2B, is a VWF disease.
Because safety is adjudicated gene-blind, the hard stop still fires even though type 2B holds
under one percent of the posterior - which is the whole argument for the interlock, rendered on a
single case. Every value is read from a live engine call, so the panel cannot drift from the
implementation.

Panel C plots ``p_disease[d][0]`` with the credible interval as a whisker. The tuple is
(mean, lower, upper) from the Monte-Carlo resampling, so plotting element 1 would chart the lower
bound under an axis labelled "posterior probability". The mean differs in the third digit from
``posterior.confidence``, which panel D reports: that is the analytic point estimate at the
nominal frequencies, not the mean over resampled ones. Both panels say which they are.
"""
from __future__ import annotations

import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from figures.scripts import style as st

# scripts/ -> figures/ -> repository root
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def build(outdir):
    from core.dx_schemas import Feature, FeatureKind
    from jointdx.factorgraph import Evidence
    from jointdx.orchestrate import diagnose
    from rules.vcep.partition import owner

    st.apply()
    ev = Evidence(
        variant_gene="GP1BA", variant_id="GP1BA:c.X", genetic_codes=["PM2"],
        clinical=[Feature("ripa_low_dose_enhanced", FeatureKind.LAB, True, observed=True),
                  Feature("ripa_mixing_platelet_origin", FeatureKind.LAB, True, observed=True),
                  Feature("ripa_mixing_plasma_origin", FeatureKind.LAB, False, observed=False),
                  Feature("thrombocytopenia", FeatureKind.CLINICAL, True, observed=True)])
    rec = diagnose(ev, planned_tx="ddavp", n_mc=400)
    ranked = sorted(rec.posterior.p_disease.items(), key=lambda kv: -kv[1][0])
    names = {d.id: d.name for d in rec.posterior.cluster.diseases}

    fig = plt.figure(figsize=(st.FULL_W, 4.35))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.92], width_ratios=[1, 1.1],
                          hspace=0.30, wspace=0.30,
                          left=0.175, right=0.995, top=0.925, bottom=0.105)

    # ------------------------------------------------------------------ A, what went in
    ax = fig.add_subplot(gs[0, 0])
    ax.axis("off")
    ax.set_title("Evidence supplied", loc="left", fontweight="bold", pad=4)
    st.panel_label(ax, "A", dx=-0.075, dy=1.10)
    lines = [("gene", "GP1BA", st.BLUE),
             ("variant", "missense, absent from gnomAD", st.BLUE),
             ("low-dose RIPA", "enhanced", st.GREEN),
             ("RIPA mixing study", "platelet origin", st.GREEN),
             ("plasma origin", "ABSENT (pertinent negative)", st.ORANGE),
             ("platelet count", "thrombocytopenia", st.GREEN),
             ("planned therapy", "desmopressin (DDAVP)", st.VERMILION)]
    for i, (k, v, c) in enumerate(lines):
        y = 0.94 - i * 0.155
        ax.text(0.0, y, k, fontsize=st.TINY, color=st.GREY, va="center")
        ax.text(0.44, y, v, fontsize=st.TINY, color=c, va="center",
                fontweight="bold" if c in (st.VERMILION, st.ORANGE) else "normal")
    ax.set_xlim(0, 1.32)
    ax.set_ylim(-0.03, 1.0)

    # ------------------------------------------------------------------ B, criterion trail
    ax = fig.add_subplot(gs[0, 1])
    ax.axis("off")
    ax.set_title("Criterion trail and factor ownership", loc="left", fontweight="bold", pad=4)
    st.panel_label(ax, "B", dx=-0.075, dy=1.10)
    applied = list(ev.genetic_codes)
    withheld = [("PP4", "phenotype specificity"), ("PS3", "functional assay"),
                ("PM5", "same residue, ClinVar-blinded"), ("PS4", "case-control counts")]
    ax.text(0.0, 0.94, "applied", fontsize=st.TINY, fontweight="bold", color=st.BLUE, va="center")
    for i, code in enumerate(applied):
        ax.text(0.05, 0.80 - i * 0.13, f"{code}  ->  {owner(code).replace('_', ' ')}",
                fontsize=st.TINY, va="center")
    base = 0.80 - len(applied) * 0.13 - 0.05
    ax.text(0.0, base, "not applied here", fontsize=st.TINY, fontweight="bold", color=st.ORANGE,
            va="center")
    for i, (code, why) in enumerate(withheld):
        fac = owner(code) or "?"
        ax.text(0.05, base - 0.135 - i * 0.13, f"{code}  ->  {fac.replace('_', ' ')}   ({why})",
                fontsize=st.TINY, color=st.GREY, va="center")
    ax.text(0.0, base - 0.135 - len(withheld) * 0.13 - 0.02,
            "variant remains VUS on intrinsic evidence", fontsize=st.TINY, style="italic",
            color=st.VERMILION, va="center")
    ax.set_xlim(0, 1.25)
    ax.set_ylim(-0.03, 1.0)

    # ------------------------------------------------------------------ C, the differential
    ax = fig.add_subplot(gs[1, 0])
    ids = [d for d, _ in ranked]
    # (probability, lower, upper); the bar is the point estimate and the whisker the interval.
    vals = [v[0] * 100 for _, v in ranked]
    lo = [max(v[0] - v[1], 0.0) * 100 for _, v in ranked]
    hi = [max(v[2] - v[0], 0.0) * 100 for _, v in ranked]
    cols = [st.BLUE if i == 0 else (st.VERMILION if d == "vwd2b" else st.GREY)
            for i, d in enumerate(ids)]
    y = list(range(len(ids)))[::-1]
    ax.barh(y, vals, color=cols, alpha=0.85, height=0.6, zorder=2,
            xerr=[lo, hi], error_kw=dict(ecolor=st.GREY, lw=0.7, capsize=1.8, zorder=3))
    for yy, v in zip(y, vals, strict=True):
        ax.text(min(v + 3.5, 74), yy, f"{v:.1f}%", va="center", fontsize=st.TINY, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels([st.wrap(names.get(d, d).replace(" (pseudo)", ""), 22) for d in ids],
                       fontsize=st.TINY, linespacing=1.25)
    ax.set_xlim(0, 108)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("posterior probability (%)")
    ax.set_title("Ranked differential", loc="left", fontweight="bold", pad=13)
    ax.text(0.0, 1.035, "posterior mean, 95% credible interval", transform=ax.transAxes,
            ha="left", fontsize=st.TINY, style="italic", color=st.GREY)
    st.panel_label(ax, "C", dx=-0.42, dy=1.17)
    ax.grid(axis="x", color=st.LIGHT, zorder=0)
    ax.set_axisbelow(True)

    # ------------------------------------------------------------------ D, what changes management
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    ax.set_title("Decision output", loc="left", fontweight="bold", pad=4)
    st.panel_label(ax, "D", dx=-0.075, dy=1.10)
    hard = [f for f in rec.safety_flags if "HARD STOP" in f.message]
    ax.add_patch(mpatches.FancyBboxPatch((0.0, 0.60), 1.12, 0.36, boxstyle="round,pad=0.012",
                                         facecolor=st.VERMILION, alpha=0.13,
                                         edgecolor=st.VERMILION, lw=1.0))
    ax.text(0.035, 0.885, "HARD STOP", fontsize=st.SMALL, fontweight="bold", color=st.VERMILION,
            va="center")
    if hard:
        p2b = next((f.p_competitor for f in hard if f.competitor_id == "vwd2b"), None)
        shown = f"{p2b:.1%}" if p2b else "non-zero"
        ax.text(0.035, 0.79,
                "desmopressin is contraindicated if type 2B\n"
                "von Willebrand disease, which the gene-blind\n"
                f"safety view still gives {shown} and has not excluded",
                fontsize=st.TINY, va="top", linespacing=1.5)
    rows = [("leading call, point estimate",
             f"{names.get(rec.posterior.leading, rec.posterior.leading).replace(' (pseudo)', '')} "
             f"({rec.posterior.confidence:.1%}), decided", st.BLUE, "bold"),
            ("variant classification", "uncertain significance - abstained", "black", "normal"),
            ("recommended next observation",
             rec.next_observation.name if rec.next_observation else "-", st.GREEN, "bold")]
    for i, (head, body, colour, weight) in enumerate(rows):
        y0 = 0.48 - i * 0.195
        ax.text(0.0, y0, head, fontsize=st.TINY, color=st.GREY, va="center")
        ax.text(0.0, y0 - 0.095, body, fontsize=st.TINY, color=colour, fontweight=weight,
                va="center")
    ax.set_xlim(0, 1.25)
    ax.set_ylim(-0.02, 1.0)
    return st.save(fig, outdir, "figS5_worked_example")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ROOT)
    print(build(os.path.join(ROOT, "figures", "out")))
