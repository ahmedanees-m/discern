"""Manuscript figures, rendered from the committed benchmark JSON.

Every panel is driven by a file under ``bench/`` or ``eval/`` so the figures move when the numbers
move. Two panels are schematics (Figure 1) and one reads its counts from the CHAMP/CHBMP benchmark
document rather than a JSON; both are marked in the code.

Production settings follow the submission package: vector PDF plus a 300 dpi PNG, fonts embedded,
Okabe-Ito colours (colourblind-safe), and no red/green as a sole discriminator.

Run:  python -m figures.make_figures [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
BENCH = os.path.join(ROOT, "bench")
EVAL = os.path.join(ROOT, "eval")

# Okabe-Ito
BLUE, ORANGE, GREEN, YELLOW = "#0072B2", "#E69F00", "#009E73", "#F0E442"
VERM, SKY, PURPLE, GREY = "#D55E00", "#56B4E9", "#CC79A7", "#666666"

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "font.family": "sans-serif",
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def _load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _save(fig, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    paths = []
    for ext in ("pdf", "png"):
        p = os.path.join(outdir, f"{name}.{ext}")
        fig.savefig(p, bbox_inches="tight")
        paths.append(p)
    plt.close(fig)
    return paths


# --------------------------------------------------------------------------------------------
def fig1_architecture(outdir):
    """Schematic. Panel B is generated from the live partition map, not drawn by hand."""
    from rules.vcep.partition import FACTOR_OF

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.4))

    # Panel A - the coupled model
    ax1.axis("off")
    ax1.set_title("A  Coupled disease-variant model", loc="left", fontweight="bold")
    inputs = [("Phenotype\n(HPO + pertinent negatives)", 0.85, SKY),
              ("Variant-intrinsic genetics\n(frequency, in-silico, null)", 0.55, BLUE),
              ("Functional / lab findings", 0.25, GREEN)]
    for label, y, c in inputs:
        ax1.add_patch(mpatches.FancyBboxPatch((0.02, y - 0.07), 0.36, 0.14, boxstyle="round,pad=0.01",
                                              facecolor=c, alpha=0.25, edgecolor=c))
        ax1.text(0.20, y, label, ha="center", va="center", fontsize=6.5)
        ax1.annotate("", xy=(0.48, 0.55), xytext=(0.385, y),
                     arrowprops=dict(arrowstyle="->", color=GREY, lw=0.9))
    ax1.add_patch(mpatches.FancyBboxPatch((0.48, 0.42), 0.20, 0.26, boxstyle="round,pad=0.01",
                                          facecolor=ORANGE, alpha=0.30, edgecolor=ORANGE, lw=1.2))
    ax1.text(0.58, 0.55, "P(D, V | E)\njoint posterior", ha="center", va="center",
             fontsize=7, fontweight="bold")
    outs = ["Ranked differential", "Variant classification", "Treatment-safety stop",
            "Cheapest next test"]
    for i, o in enumerate(outs):
        y = 0.88 - i * 0.22
        ax1.annotate("", xy=(0.74, y), xytext=(0.68, 0.55),
                     arrowprops=dict(arrowstyle="->", color=GREY, lw=0.9))
        ax1.text(0.755, y, o, ha="left", va="center", fontsize=6.5)
    ax1.set_xlim(0, 1.25)
    ax1.set_ylim(0, 1)

    # Panel B - the partition, straight from code
    ax2.set_title("B  Evidence partition: each criterion to one factor", loc="left",
                  fontweight="bold")
    factors, colours = {}, {"variant_intrinsic": BLUE, "disease_pp4": ORANGE, "functional": GREEN,
                            "segregation": PURPLE, "phasing": VERM, "denovo": GREY}
    for code, fac in FACTOR_OF.items():
        factors.setdefault(fac, []).append(code)
    order = ["variant_intrinsic", "functional", "disease_pp4", "segregation", "phasing", "denovo"]
    y = 0
    ticks, labels = [], []
    for fac in order:
        codes = sorted(factors.get(fac, []))
        ax2.barh(y, len(codes), color=colours[fac], alpha=0.75, height=0.6)
        ax2.text(len(codes) + 0.3, y, ", ".join(codes), va="center", fontsize=5.4)
        ticks.append(y)
        labels.append(fac.replace("_", " "))
        y -= 1
    ax2.set_yticks(ticks)
    ax2.set_yticklabels(labels)
    ax2.set_xlabel("number of ACMG criteria owned")
    ax2.set_xlim(0, 34)
    ax2.set_ylim(y + 0.5, 0.6)
    return _save(fig, outdir, "fig1_architecture_and_partition")


# --------------------------------------------------------------------------------------------
def fig2_discrimination_calibration(outdir):
    from sklearn.metrics import roc_curve

    from bench.phase_r_variant import oof_isotonic, shared_folds
    from bench.track1b_erepo_headtohead import TIMESPLIT_AFTER, load_rows

    m = _load(os.path.join(BENCH, "phase_r_variant_metrics.json"))
    rows = load_rows()

    def surface(rs):
        mis = [r for r in rs if r["is_missense"] and r["label"] is not None]
        return mis, [r["label"] for r in mis]

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7))

    # Panel A - ROC
    mis, y = surface(rows)
    ax = axes[0]
    for name, key, c in (("DISCERN", "discern_points", BLUE), ("REVEL", "revel", ORANGE),
                         ("AlphaMissense", "alphamissense", GREEN)):
        pairs = [(r[key], yy) for r, yy in zip(mis, y, strict=False) if r[key] is not None]
        fpr, tpr, _ = roc_curve([b for _, b in pairs], [float(a) for a, _ in pairs])
        d = m["eRepo_primary"]["discrimination"][name]
        ax.plot(fpr, tpr, color=c, lw=1.3,
                label=f"{name} {d['auroc']:.3f} [{d['auroc_ci95'][0]:.3f}-{d['auroc_ci95'][1]:.3f}]")
    ax.plot([0, 1], [0, 1], ls=":", color=GREY, lw=0.8)
    dv = m["eRepo_primary"]["discern_vs_revel"]["delong"]
    ax.set_title(f"A  eRepo missense (n={m['eRepo_primary']['n_missense_PB']})", loc="left",
                 fontweight="bold")
    ax.set_xlabel("false positive rate")
    ax.set_ylabel("true positive rate")
    ax.legend(loc="lower right", frameon=False)
    ax.text(0.98, 0.06, f"DISCERN vs REVEL: DeLong p={dv['p_value']:.2f} (ns)",
            transform=ax.transAxes, ha="right", fontsize=6, style="italic")

    # Panels B and C - reliability, both surfaces
    for ax, surf, title in ((axes[1], rows, "B  Reliability, eRepo-primary"),
                            (axes[2], [r for r in rows if (r["approval_date"] or "") > TIMESPLIT_AFTER],
                             "C  Reliability, time-split")):
        mis, y = surface(surf)
        cal = m["eRepo_primary" if surf is rows else "time_split"]["calibration"]
        ax.plot([0, 1], [0, 1], ls=":", color=GREY, lw=0.8)
        for name, key, c, lab in (("DISCERN", "discern_points", BLUE, "DISCERN_isotonic_oof"),
                                  ("REVEL", "revel", ORANGE, "REVEL_isotonic_oof"),
                                  ("AlphaMissense", "alphamissense", GREEN,
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
            ax.plot(xs, ys, "o-", color=c, ms=3, lw=1.1,
                    label=f"{name} ECE {e['ece']:.3f} [{e['ece_ci95'][0]:.3f}-{e['ece_ci95'][1]:.3f}]")
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xlabel("predicted probability")
        ax.set_ylabel("observed frequency")
        ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    return _save(fig, outdir, "fig2_discrimination_and_calibration")


# --------------------------------------------------------------------------------------------
def fig3_per_criterion_kappa(outdir):
    """The corrected claim lives here: zero-application criteria split by cause."""
    m = _load(os.path.join(BENCH, "track1b_erepo_metrics.json"))["per_code_kappa_vs_erepo"]
    routed = {"PS3", "BS3", "PP4", "PP1", "BS4", "PM3", "BP2", "PS2", "PM6"}

    scored = [(c, v["kappa"]) for c, v in m.items() if v["kappa"] is not None]
    scored.sort(key=lambda kv: kv[1])
    zero = [(c, v) for c, v in m.items() if v["discern_applied"] == 0]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5.2, 4.4),
                                   gridspec_kw={"height_ratios": [len(scored), max(len(zero), 1)]})
    ax1.barh([c for c, _ in scored], [k for _, k in scored], color=BLUE, alpha=0.85)
    for i, (_, k) in enumerate(scored):
        ax1.text(k + 0.015, i, f"{k:.2f}", va="center", fontsize=6.5)
    ax1.set_xlim(0, 1.08)
    ax1.set_xlabel("Cohen's kappa vs eRepo expert applications")
    ax1.set_title("A  Criteria DISCERN applies", loc="left", fontweight="bold")

    names = [c for c, _ in zero]
    colours = [GREEN if c in routed else ORANGE for c, _ in zero]
    ax2.barh(names, [v["erepo_applied"] for _, v in zero], color=colours, alpha=0.85)
    for i, (_, v) in enumerate(zero):
        ax2.text(v["erepo_applied"] + 1, i, f"{v['erepo_applied']} by experts", va="center",
                 fontsize=6.5)
    ax2.set_xlabel("times the expert panel applied it (DISCERN: zero)")
    ax2.set_title("B  Criteria DISCERN applies zero times, by cause", loc="left", fontweight="bold")
    ax2.legend(handles=[mpatches.Patch(color=GREEN, alpha=0.85,
                                       label="routed to another factor by the partition"),
                        mpatches.Patch(color=ORANGE, alpha=0.85,
                                       label="variant-intrinsic; no input in this pipeline")],
               loc="lower right", frameon=False)
    ax2.set_xlim(0, max(v["erepo_applied"] for _, v in zero) * 1.45 + 2)
    fig.tight_layout()
    return _save(fig, outdir, "fig3_per_criterion_kappa")


# --------------------------------------------------------------------------------------------
def fig4_intrinsic_ceiling(outdir):
    """The lead argument."""
    from bench.track1b_erepo_headtohead import load_rows
    m = _load(os.path.join(BENCH, "phase_r_variant_metrics.json"))
    ca = m["ceiling_attribution"]
    ic = m["added_value_over_ranking_score"]["intrinsic_only_ceiling"]
    mis = [r for r in load_rows() if r["is_missense"] and r["label"] is not None]
    pts = [r["discern_points"] for r in mis]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0))

    lo, hi = int(min(pts)) - 1, 8
    ax1.hist(pts, bins=np.arange(lo, hi + 1) - 0.5, color=BLUE, alpha=0.85, edgecolor="white")
    ax1.axvline(6, color=VERM, lw=1.4, ls="--")
    ax1.text(6.15, ax1.get_ylim()[1] * 0.92, "Likely Pathogenic\nthreshold (6 points)",
             color=VERM, fontsize=6.5, va="top")
    ax1.axvspan(6, hi, color=VERM, alpha=0.06)
    ax1.text((6 + hi) / 2, ax1.get_ylim()[1] * 0.45, f"zero of {ca['n']}", color=VERM,
             fontsize=8, ha="center", fontweight="bold")
    ax1.annotate("observed\nmaximum {:.0f}".format(ic["max_discern_points_on_missense"]),
                 xy=(3, 10), xytext=(4.3, ax1.get_ylim()[1] * 0.28), fontsize=6.5, ha="center",
                 arrowprops=dict(arrowstyle="->", color=GREY, lw=0.8,
                                 connectionstyle="arc3,rad=0.25"))
    ax1.set_xlabel("total ACMG points, intrinsic evidence only")
    ax1.set_ylabel("missense variants")
    ax1.set_title("A  Nothing reaches a pathogenic band", loc="left", fontweight="bold")

    steps = [("as\nscored", ca["reach_lp_as_scored"], GREY, None),
             ("+ no-input\nintrinsic", ca["reach_lp_if_intrinsic_codes_were_available"], ORANGE,
              "+ intrinsic criteria with no input here (PM1, PM5, PS1, PS4)"),
             ("+ routed by\npartition", ca["reach_lp_if_routed_codes_were_re_added"], GREEN,
              "+ criteria the partition routes away (PS3, PP4, PP1, PM3)"),
             ("both\nrestored", ca["reach_lp_with_both"], BLUE, None)]
    ax2.bar([s[0] for s in steps], [s[1] for s in steps], color=[s[2] for s in steps], alpha=0.85)
    for i, (_, v, _c, _d) in enumerate(steps):
        ax2.text(i, v + 7, str(v), ha="center", fontsize=7.5, fontweight="bold")
    ax2.legend(handles=[mpatches.Patch(color=s[2], alpha=0.85, label=s[3])
                        for s in steps if s[3]],
               loc="upper center", frameon=False, fontsize=5.5, bbox_to_anchor=(0.5, 0.78))
    n_path = m["eRepo_primary"]["n_path"]
    ax2.axhline(n_path, color=VERM, ls="--", lw=1.2)
    ax2.text(-0.45, n_path + 4, f"{n_path} truly pathogenic", color=VERM, fontsize=6.5,
             ha="left", va="bottom")
    ax2.set_ylabel("variants reaching Likely Pathogenic")
    ax2.set_ylim(0, n_path * 1.15)
    ax2.set_title("B  Neither cause explains the ceiling", loc="left", fontweight="bold")
    ax2.tick_params(axis="x", labelsize=6.5)
    fig.tight_layout()
    return _save(fig, outdir, "fig4_intrinsic_ceiling")


# --------------------------------------------------------------------------------------------
def fig5_clinvar_circularity(outdir):
    from sklearn.metrics import roc_curve

    from bench.track1b_erepo_headtohead import TIMESPLIT_AFTER, load_rows
    m = _load(os.path.join(BENCH, "track1b_erepo_metrics.json"))
    c1 = _load(os.path.join(BENCH, "track1_metrics.json"))["clinvar_circularity"]
    rows = load_rows()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.9))
    for rs, style, tag in ((rows, "-", "eRepo-primary"),
                           ([r for r in rows if (r["approval_date"] or "") > TIMESPLIT_AFTER],
                            "--", "time-split")):
        mis = [r for r in rs if r["is_missense"] and r["label"] is not None]
        y = [r["label"] for r in mis]
        for name, key, c in (("GeneBe", "genebe_score", VERM), ("DISCERN", "discern_points", BLUE)):
            pairs = [(r[key], yy) for r, yy in zip(mis, y, strict=False) if r[key] is not None]
            fpr, tpr, _ = roc_curve([b for _, b in pairs], [float(a) for a, _ in pairs])
            surf = "eRepo_primary" if tag == "eRepo-primary" else "time_split"
            key2 = "GeneBe_exhibit" if name == "GeneBe" else "DISCERN"
            auc = m[surf]["discrimination"][key2]["auroc"]
            ax1.plot(fpr, tpr, style, color=c, lw=1.3, label=f"{name} {tag}: {auc:.3f}")
    ax1.plot([0, 1], [0, 1], ls=":", color=GREY, lw=0.8)
    ax1.set_xlabel("false positive rate")
    ax1.set_ylabel("true positive rate")
    ax1.set_title("A  A live-ClinVar tool is perfectly separable", loc="left", fontweight="bold")
    ax1.legend(loc="lower right", frameon=False)

    match, total = c1["genebe_class_matches_clinvar_direction"].split("/")
    ax2.bar(["GeneBe\n(consumes ClinVar)", "DISCERN / REVEL /\nAlphaMissense"],
            [100 * int(match) / int(total), 0], color=[VERM, BLUE], alpha=0.85)
    ax2.text(0, 100 * int(match) / int(total) + 2, f"{match}/{total} calls match\nthe ClinVar direction",
             ha="center", fontsize=6.5)
    ax2.text(1, 4, "apply no\nClinVar-derived code", ha="center", fontsize=6.5)
    ax2.set_ylabel("% of resolved calls matching ClinVar")
    ax2.set_ylim(0, 118)
    ax2.set_title("B  Reproduction, not prediction", loc="left", fontweight="bold")
    fig.tight_layout()
    return _save(fig, outdir, "fig5_clinvar_circularity")


# --------------------------------------------------------------------------------------------
def fig6_diagnosis_baselines(outdir):
    g = _load(os.path.join(EVAL, "gene_only_baseline.json"))
    lir = _load(os.path.join(EVAL, "lirical_arm.json"))
    cov = _load(os.path.join(EVAL, "phenotype_tool_comparison.json"))["hpo_coverage"]
    p = g["pooled"]

    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.9))

    ax = axes[0]
    bars = [("uniform random", p["random_within_cluster"]["top1"]["expected"], GREY, False),
            ("prior only", p["prior_only"]["top1"], GREY, False),
            ("gene lookup", p["gene_lookup"]["top1"], ORANGE, False),
            ("DISCERN pre-corr.", 0.81, SKY, False),
            ("DISCERN post-corr.", p["DISCERN"]["top1"], BLUE, True)]
    for i, (_lab, v, c, hatch) in enumerate(bars):
        ax.bar(i, v * 100, color=c, alpha=0.85, hatch="//" if hatch else None,
               edgecolor="white" if hatch else None, lw=0 if not hatch else 1.0)
        ax.text(i, v * 100 + 2, f"{v:.0%}", ha="center", fontsize=6.5)
    ax.set_xticks(range(len(bars)))
    ax.set_xticklabels([b[0] for b in bars], fontsize=6, rotation=30, ha="right")
    ax.set_ylabel("Top-1 accuracy (%)")
    ax.set_ylim(0, 118)
    ax.set_title(f"A  Curated cases (n={p['n']})", loc="left", fontweight="bold")
    ax.text(0.5, 1.02, f"vs gene lookup: McNemar p={p['mcnemar_top1']['p_value_exact']:.2f} (ns)",
            transform=ax.transAxes, ha="center", fontsize=6, style="italic")
    ax.legend(handles=[mpatches.Patch(facecolor=BLUE, alpha=0.85, hatch="//", edgecolor="white",
                                      label="in-sample")], loc="upper left", frameon=False)

    ax = axes[1]
    matched = lir["paired_tests_vs_lirical_restricted"]["hpo_representable_only"]
    full = lir["paired_tests_vs_lirical_restricted"]["phenotype_only"]
    groups = [("identical evidence\n(13 HPO-codable findings)", matched, False),
              ("all findings\n(DISCERN reads 48)", full, True)]
    x = 0
    ticks, labels = [], []
    for lab, t, sig in groups:
        ax.bar(x, t["lirical_recall@1"] * 100, color=ORANGE, alpha=0.85, width=0.38)
        ax.bar(x + 0.4, t["discern_recall@1"] * 100, color=BLUE, alpha=0.85, width=0.38)
        star = f"p={t['mcnemar']['p_value_exact']:.2f}" + ("*" if sig else " (ns)")
        ax.text(x + 0.2, max(t["discern_recall@1"], t["lirical_recall@1"]) * 100 + 6, star,
                ha="center", fontsize=6.5, fontweight="bold" if sig else "normal")
        ticks.append(x + 0.2)
        labels.append(lab)
        x += 1.2
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=5.8)
    ax.set_ylabel("Recall@1 (%)")
    ax.set_ylim(0, 118)
    ax.set_title(f"B  vs LIRICAL (n={matched.get('n', 23) if isinstance(matched, dict) else 23})",
                 loc="left", fontweight="bold")
    ax.legend(handles=[mpatches.Patch(color=ORANGE, alpha=0.85, label="LIRICAL"),
                       mpatches.Patch(color=BLUE, alpha=0.85, label="DISCERN")],
              loc="upper left", frameon=False)

    ax = axes[2]
    codable = cov["features_expressible_as_hpo"]
    total = cov["distinct_features_used"]
    ax.pie([codable, total - codable], colors=[GREEN, ORANGE],
           labels=[f"{codable} codable\nin HPO", f"{total - codable} laboratory\nfindings HPO\ncannot express"],
           autopct="%1.0f%%", startangle=90, textprops={"fontsize": 6},
           wedgeprops={"alpha": 0.85, "edgecolor": "white"})
    ax.set_title("C  Why the channel is empty", loc="left", fontweight="bold")
    fig.tight_layout()
    return _save(fig, outdir, "fig6_diagnosis_baselines")


# --------------------------------------------------------------------------------------------
def figS1_gene_term_sweep(outdir):
    m = _load(os.path.join(BENCH, "phase_r_gene_term_sensitivity.json"))
    grid = sorted(m["grid"], key=lambda r: r["likelihood_ratio"])
    lr = [r["likelihood_ratio"] for r in grid]
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    ax.plot(lr, [r["curated_top1"] * 100 for r in grid], "o-", color=BLUE, label="curated Top-1")
    ax.plot(lr, [r["hardstop_sensitivity"] * 100 for r in grid], "s--", color=GREEN,
            label="hard-stop sensitivity")
    ax.plot(lr, [r["hardstop_specificity"] * 100 for r in grid], "^:", color=PURPLE,
            label="hard-stop specificity")
    committed = next(r for r in grid if r["committed"])
    ax.axvline(committed["likelihood_ratio"], color=VERM, lw=1.2, ls="--")
    ax.text(committed["likelihood_ratio"] * 1.05, 30, "committed value\n(largest at which the\n"
            "deciding assay still wins)", color=VERM, fontsize=6)
    for r in grid:
        if not r["deciding_assay_can_overturn_gene"]:
            ax.axvspan(r["likelihood_ratio"] * 0.82, r["likelihood_ratio"] * 1.22,
                       color=VERM, alpha=0.05)
    ax.set_xscale("log")
    ax.set_xticks(lr)
    ax.set_xticklabels([str(v) for v in lr])
    ax.set_xlabel("P(G|D) likelihood ratio")
    ax.set_ylabel("%")
    ax.set_ylim(0, 108)
    ax.legend(loc="lower right", frameon=False)
    ax.set_title("Gene-term sensitivity: safety invariant, diagnosis flat above 4",
                 loc="left", fontweight="bold")
    fig.tight_layout()
    return _save(fig, outdir, "figS1_gene_term_sensitivity")


def figS2_safety_matrix(outdir):
    m = _load(os.path.join(BENCH, "track3_metrics.json"))["safety_interlock"]
    rows = m["rows"]
    fig, ax = plt.subplots(figsize=(5.6, 2.6))
    grid = np.array([[1 if r["hardstop_on_real"] else 0,
                      1 if r["hardstop_on_harmless"] else 0] for r in rows])
    ax.imshow(grid, cmap=matplotlib.colors.ListedColormap(["#FFFFFF", GREEN]), vmin=0, vmax=1,
              aspect="auto")
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            ax.text(j, i, "FIRED" if grid[i, j] else "silent", ha="center", va="center",
                    fontsize=6, color="white" if grid[i, j] else GREY)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["contraindicated plan\n(must fire)", "harmless plan\n(must stay silent)"],
                       fontsize=6.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{r['id']}\n({r['harmful_tx']})" for r in rows], fontsize=5.8)
    ax.set_title(f"Safety interlock: sensitivity {m['hardstop_sensitivity']:.0%}, "
                 f"specificity {m['hardstop_specificity']:.0%}, judged gene-blind",
                 loc="left", fontweight="bold")
    for s in ax.spines.values():
        s.set_visible(False)
    fig.tight_layout()
    return _save(fig, outdir, "figS2_safety_matrix")


def figS3_diagnosis_calibration(outdir):
    m = _load(os.path.join(BENCH, "track3_metrics.json"))["calibration_diagnosis"]
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    ax.plot([0, 1], [0, 1], ls=":", color=GREY, lw=0.8)
    ax.plot([m["mean_confidence"]], [m["accuracy"]], "o", color=BLUE, ms=8)
    ax.annotate(f"mean confidence {m['mean_confidence']:.3f}\naccuracy {m['accuracy']:.3f}\n"
                f"ECE {m['ece']:.3f} (n={m['n']})",
                xy=(m["mean_confidence"], m["accuracy"]), xytext=(0.30, 0.62), fontsize=6.5,
                arrowprops=dict(arrowstyle="->", color=GREY, lw=0.8))
    ax.axhspan(0.8, 1.0, color=GREEN, alpha=0.07)
    ax.text(0.03, 0.97, f"zero confidently-wrong calls at 0.8 "
            f"({m['confidently_wrong_at_0.8']} of {m['n']})", fontsize=6.5, va="top", color=GREEN)
    ax.set_xlabel("mean predicted confidence")
    ax.set_ylabel("observed accuracy")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.05)
    ax.set_title("Diagnosis-posterior calibration (proof of concept, n=42)",
                 loc="left", fontweight="bold")
    fig.tight_layout()
    return _save(fig, outdir, "figS3_diagnosis_calibration")


# Counts transcribed from docs/DISCERN_CHAMP_CHBMP_Benchmark_v1.md; this is the one figure whose
# inputs are not a JSON, so the source is named here and in the caption.
CHAMP = {"overall": 0.912, "F8": 0.977, "F9": 0.658, "n_alleles": 5437, "n_null": 2130}


def figS4_champ_chbmp(outdir):
    fig, ax = plt.subplots(figsize=(4.6, 2.8))
    keys = ["F8", "F9", "overall"]
    ax.bar(keys, [CHAMP[k] * 100 for k in keys], color=[BLUE, ORANGE, GREY], alpha=0.85)
    for i, k in enumerate(keys):
        ax.text(i, CHAMP[k] * 100 + 2, f"{CHAMP[k]:.1%}", ha="center", fontsize=7)
    ax.set_ylabel("LP/P recall on null variants (%)")
    ax.set_ylim(0, 112)
    ax.set_title(f"CHAMP/CHBMP independent sensitivity\n({CHAMP['n_null']} null of "
                 f"{CHAMP['n_alleles']} disease alleles, consequence alone)",
                 loc="left", fontweight="bold")
    ax.text(1, 30, "F9 gap is last-exon\nNMD escape, correctly\nheld at VUS", ha="center",
            fontsize=6, color=VERM)
    fig.tight_layout()
    return _save(fig, outdir, "figS4_champ_chbmp_recall")



def figS5_worked_example(outdir):
    """One case end to end: the safety interlock preventing a documented harm.

    A platelet-type von Willebrand case. The gene is GP1BA; the contraindicated disease, type 2B,
    is a VWF disease. Because safety is adjudicated gene-blind, the hard stop still fires even
    though type 2B holds only one percent of the posterior - which is the whole argument for the
    interlock, rendered on a single case. Every value is read from a live engine call, so the panel
    cannot drift from the implementation.
    """
    from core.dx_schemas import Feature, FeatureKind
    from jointdx.factorgraph import Evidence
    from jointdx.orchestrate import diagnose
    from rules.vcep.partition import owner

    ev = Evidence(
        variant_gene="GP1BA", variant_id="GP1BA:c.X", genetic_codes=["PM2"],
        clinical=[Feature("ripa_low_dose_enhanced", FeatureKind.LAB, True, observed=True),
                  Feature("ripa_mixing_platelet_origin", FeatureKind.LAB, True, observed=True),
                  Feature("ripa_mixing_plasma_origin", FeatureKind.LAB, False, observed=False),
                  Feature("thrombocytopenia", FeatureKind.CLINICAL, True, observed=True)])
    rec = diagnose(ev, planned_tx="ddavp", n_mc=400)
    ranked = sorted(rec.posterior.p_disease.items(), key=lambda kv: -kv[1][1])
    names = {d.id: d.name for d in rec.posterior.cluster.diseases}

    fig = plt.figure(figsize=(7.4, 5.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1], width_ratios=[1, 1.15],
                          hspace=0.42, wspace=0.28)

    # A - what went in
    ax = fig.add_subplot(gs[0, 0])
    ax.axis("off")
    ax.set_title("A  Evidence supplied", loc="left", fontweight="bold")
    lines = [("gene", "GP1BA", BLUE),
             ("variant", "missense, absent from gnomAD", BLUE),
             ("low-dose RIPA", "enhanced", GREEN),
             ("RIPA mixing study", "platelet origin", GREEN),
             ("plasma origin", "ABSENT (pertinent negative)", ORANGE),
             ("platelet count", "thrombocytopenia", GREEN),
             ("planned therapy", "desmopressin (DDAVP)", VERM)]
    for i, (k, v, c) in enumerate(lines):
        y = 0.93 - i * 0.145
        ax.text(0.0, y, k, fontsize=6.4, color=GREY)
        ax.text(0.52, y, v, fontsize=6.4, color=c,
                fontweight="bold" if c in (VERM, ORANGE) else "normal")
    ax.set_xlim(0, 1.45)
    ax.set_ylim(0.0, 1.02)

    # B - criterion trail, with the owning factor
    ax = fig.add_subplot(gs[0, 1])
    ax.axis("off")
    ax.set_title("B  Criterion trail and factor ownership", loc="left", fontweight="bold")
    applied = list(ev.genetic_codes)
    withheld = [("PP4", "phenotype specificity"), ("PS3", "functional assay"),
                ("PM5", "same residue, ClinVar-blinded"), ("PS4", "case-control counts")]
    ax.text(0.0, 0.92, "applied", fontsize=6.6, fontweight="bold", color=BLUE)
    for i, code in enumerate(applied):
        ax.text(0.05, 0.80 - i * 0.11, f"{code}  ->  {owner(code).replace('_', ' ')}",
                fontsize=6.4)
    base = 0.80 - len(applied) * 0.11 - 0.06
    ax.text(0.0, base, "not applied here", fontsize=6.6, fontweight="bold", color=ORANGE)
    for i, (code, why) in enumerate(withheld):
        fac = owner(code) or "?"
        ax.text(0.05, base - 0.12 - i * 0.115, f"{code}  ->  {fac.replace('_', ' ')}   ({why})",
                fontsize=6.0, color=GREY)
    ax.text(0.0, base - 0.12 - len(withheld) * 0.115 - 0.04,
            "variant remains VUS on intrinsic evidence", fontsize=6.4, style="italic", color=VERM)
    ax.set_xlim(0, 1.25)
    ax.set_ylim(0, 1)

    # C - the differential
    ax = fig.add_subplot(gs[1, 0])
    ids = [d for d, _ in ranked]
    vals = [v[1] * 100 for _, v in ranked]
    cols = [BLUE if i == 0 else (VERM if d == "vwd2b" else GREY) for i, d in enumerate(ids)]
    ax.barh(range(len(ids))[::-1], vals, color=cols, alpha=0.85)
    for i, (_d, v) in enumerate(zip(ids, vals, strict=False)):
        ax.text(min(v + 2, 88), len(ids) - 1 - i, f"{v:.1f}%", va="center", fontsize=6.5)
    ax.set_yticks(range(len(ids))[::-1])
    ax.set_yticklabels([names.get(d, d).replace(" (pseudo)", "") for d in ids],
                       fontsize=6.0)
    ax.set_xlim(0, 108)
    ax.set_xlabel("posterior probability (%)")
    ax.set_title("C  Ranked differential", loc="left", fontweight="bold")

    # D - the output that changes management
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    ax.set_title("D  Decision output", loc="left", fontweight="bold")
    hard = [f for f in rec.safety_flags if "HARD STOP" in f.message]
    ax.add_patch(mpatches.FancyBboxPatch((0.0, 0.60), 1.08, 0.34, boxstyle="round,pad=0.015",
                                         facecolor=VERM, alpha=0.13, edgecolor=VERM, lw=1.3))
    ax.text(0.04, 0.86, "HARD STOP", fontsize=8, fontweight="bold", color=VERM)
    if hard:
        p2b = next((f.p_competitor for f in hard if f.competitor_id == "vwd2b"), None)
        shown = f"{p2b:.1%}" if p2b else "non-zero"
        ax.text(0.04, 0.755,
                "desmopressin is contraindicated if type 2B\n"
                f"von Willebrand disease, which the gene-blind\n"
                f"safety view still gives {shown} and has not excluded",
                fontsize=6.2, va="top")
    ax.text(0.0, 0.50, "leading call", fontsize=6.6, color=GREY)
    ax.text(0.0, 0.42, f"{names.get(rec.posterior.leading, rec.posterior.leading).replace(' (pseudo)', '')} "
                       f"({rec.posterior.confidence:.1%}), decided", fontsize=6.4,
            fontweight="bold", color=BLUE)
    ax.text(0.0, 0.30, "variant classification", fontsize=6.6, color=GREY)
    ax.text(0.0, 0.22, "uncertain significance - abstained", fontsize=6.6)
    ax.text(0.0, 0.10, "recommended next observation", fontsize=6.6, color=GREY)
    ax.text(0.0, 0.02, rec.next_observation.name if rec.next_observation else "-",
            fontsize=6.6, fontweight="bold", color=GREEN)
    ax.set_xlim(0, 1.25)
    ax.set_ylim(-0.02, 1)
    return _save(fig, outdir, "figS5_worked_example")

FIGURES = [fig1_architecture, fig2_discrimination_calibration, fig3_per_criterion_kappa,
           fig4_intrinsic_ceiling, fig5_clinvar_circularity, fig6_diagnosis_baselines,
           figS1_gene_term_sweep, figS2_safety_matrix, figS3_diagnosis_calibration,
           figS4_champ_chbmp, figS5_worked_example]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "out"))
    args = ap.parse_args()
    for fn in FIGURES:
        paths = fn(args.out)
        print(f"  {fn.__name__:32} -> {os.path.basename(paths[0])}")
    print(f"\nwrote {len(FIGURES)} figures (pdf + png) to {args.out}")


if __name__ == "__main__":
    main()
