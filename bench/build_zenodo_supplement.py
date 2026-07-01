"""Assemble the DISCERN Paper 1 Zenodo supporting-data supplement (DISCERN_Zenodo_Deposit_Manifest).

Deposits DERIVED results + variant IDs + public labels + manifests + configs + figures. It never
copies a raw third-party database (eRepo/ClinVar/gnomAD/REVEL/AlphaMissense/CHAMP-CHBMP), never any
patient-level data, and never copyrighted article text (the curated benchmark ships as PMID + feature
codes + expected diagnosis only). Everything here regenerates from the committed engine + metrics.

Run:  python -m bench.build_zenodo_supplement [OUTDIR]
Default OUTDIR:  ../Zenodo_files  (a sibling folder next to the repository root)
"""

from __future__ import annotations

import csv
import glob
import json
import os
import shutil
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402
from sklearn.metrics import roc_curve  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = (
    os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.join(REPO, "..", "Zenodo_files")
)


def o(*a):
    return os.path.join(OUT, *a)


def r(*a):
    return os.path.join(REPO, *a)


def mk(*a):
    os.makedirs(o(*a), exist_ok=True)


TREE = [
    "code",
    "configs/vcep_specs",
    "configs/cluster_definitions",
    "benchmarks/variant_arm",
    "benchmarks/trustworthiness",
    "benchmarks/diagnosis_arm",
    "benchmarks/champ_chbmp",
    "benchmarks/coupling_poc",
    "data_manifests",
    "figures",
    "tables",
    "preregistration",
    "docs",
]


def build_tree():
    os.makedirs(OUT, exist_ok=True)
    for d in TREE:
        mk(d)


def copy_direct():
    for name in [
        "DISCERN_Benchmark_Results_v1.md",
        "DISCERN_Reproducibility_Checklist.md",
        "DISCERN_v3.1_Claims_Map.md",
        "DISCERN_Execution_Summary_v3.1.md",
    ]:
        shutil.copy(r("docs", name), o("docs", name))
    shutil.copy(
        r("docs", "DISCERN_OSF_PreRegistration_v1.md"),
        o("preregistration", "DISCERN_OSF_PreRegistration_v1.md"),
    )
    for f in glob.glob(r("rules", "vcep", "specs", "*.yaml")):
        shutil.copy(f, o("configs", "vcep_specs", os.path.basename(f)))
    for f in glob.glob(r("diseases", "clusters", "*.yaml")):
        shutil.copy(f, o("configs", "cluster_definitions", os.path.basename(f)))
    shutil.copy(r("LICENSE"), o("LICENSE-code.txt"))
    shutil.copy(r("CITATION.cff"), o("CITATION.cff"))
    shutil.copy(
        r("bench", "track1_variant_headtohead.csv"),
        o("benchmarks", "variant_arm", "track1_clinvar_blinded.csv"),
    )
    shutil.copy(
        r("bench", "track3_risk_coverage.csv"),
        o("benchmarks", "trustworthiness", "risk_coverage_curve.csv"),
    )
    shutil.copy(
        r("bench", "data", "erepo_bleeding.tsv"),
        o("benchmarks", "variant_arm", "erepo_variant_set.tsv"),
    )
    shutil.copy(
        r("eval", "champ_chbmp_ps1_revel.tsv"),
        o("benchmarks", "champ_chbmp", "champ_chbmp_ps1_revel.tsv"),
    )


def _jload(name):
    return json.load(open(r("bench", name), encoding="utf-8"))


def gen_variant_arm():
    from bench.track1_variant_headtohead import _disc_prob_oof
    from bench.track1b_erepo_headtohead import TIMESPLIT_AFTER, load_rows

    rows = load_rows()
    # per-variant eRepo head-to-head
    with open(
        o("benchmarks", "variant_arm", "track1b_erepo_headtohead.csv"),
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "gene",
                "assertion",
                "label_pathogenic",
                "is_missense",
                "approval_date",
                "in_time_split",
                "discern_points",
                "revel",
                "alphamissense",
                "genebe_acmg_score",
            ]
        )
        for x in rows:
            ts = (x["approval_date"] or "") > TIMESPLIT_AFTER
            w.writerow(
                [
                    x["gene"],
                    x["assertion"],
                    x["label"],
                    x["is_missense"],
                    x["approval_date"],
                    int(ts),
                    x["discern_points"],
                    x["revel"],
                    x["alphamissense"],
                    x["genebe_score"],
                ]
            )
    # per-ACMG-code kappa
    m1b = _jload("track1b_erepo_metrics.json")
    with open(
        o("benchmarks", "variant_arm", "per_acmg_code_kappa.csv"), "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.writer(fh)
        w.writerow(["acmg_code", "erepo_applied", "discern_applied", "cohen_kappa"])
        for code, v in m1b["per_code_kappa_vs_erepo"].items():
            w.writerow([code, v["erepo_applied"], v["discern_applied"], v["kappa"]])
    # genebe circularity exhibit
    m1 = _jload("track1_metrics.json")
    c1, c1b = m1["clinvar_circularity"], m1b["eRepo_primary"]["discrimination"]["GeneBe_exhibit"]
    with open(
        o("benchmarks", "variant_arm", "genebe_circularity_exhibit.csv"),
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "surface",
                "genebe_auroc",
                "class_matches_clinvar_direction",
                "frac_named_pp5_bp6",
                "note",
            ]
        )
        w.writerow(
            [
                f"ClinVar P/B (n={m1['n_PB']})",
                c1["genebe_auroc_overall"],
                c1["genebe_class_matches_clinvar_direction"],
                c1["frac_named_clinvar_code"],
                "reproduces ClinVar; not a fair surface",
            ]
        )
        w.writerow(
            [
                f"ClinVar PP5/BP6-name-blind (n={c1['genebe_auroc_named_clinvar_code_blind']['n']})",
                c1["genebe_auroc_named_clinvar_code_blind"]["auroc"],
                "",
                "",
                "still ~1.0 -> dependence is systemic, not just named codes",
            ]
        )
        w.writerow(
            [
                "eRepo expert-panel missense",
                c1b["auroc"] if c1b else "",
                "",
                "",
                "still ~1.0 on the independent expert surface (exhibit only)",
            ]
        )
    # calibration reliability bins (DISCERN isotonic OOF vs REVEL vs AlphaMissense) on eRepo missense P/B
    mis = [x for x in rows if x["is_missense"] and x["label"] is not None]
    y = [x["label"] for x in mis]
    disc = _disc_prob_oof([x["discern_points"] for x in mis], y)
    series = {
        "DISCERN_isotonic": disc,
        "REVEL": [x["revel"] for x in mis],
        "AlphaMissense": [x["alphamissense"] for x in mis],
    }
    with open(
        o("benchmarks", "variant_arm", "calibration_reliability.csv"),
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        w = csv.writer(fh)
        w.writerow(["tool", "bin_low", "bin_high", "n", "mean_predicted", "observed_fraction"])
        edges = np.linspace(0, 1, 11)
        for name, sc in series.items():
            pr = [(float(s), yy) for s, yy in zip(sc, y, strict=False) if s is not None]
            p = np.array([a for a, _ in pr])
            yy = np.array([b for _, b in pr])
            for i in range(10):
                hi = edges[i + 1] if i < 9 else 1.0001
                mask = (p >= edges[i]) & (p < hi)
                if mask.sum():
                    w.writerow(
                        [
                            name,
                            round(edges[i], 2),
                            round(edges[i + 1], 2),
                            int(mask.sum()),
                            round(float(p[mask].mean()), 4),
                            round(float(yy[mask].mean()), 4),
                        ]
                    )
    return rows, m1, m1b


def gen_trustworthiness():
    m3 = _jload("track3_metrics.json")
    m1b = _jload("track1b_erepo_metrics.json")
    # calibration metrics (variant + diagnosis)
    with open(
        o("benchmarks", "trustworthiness", "calibration_metrics.csv"),
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        w = csv.writer(fh)
        w.writerow(["arm", "tool", "surface", "ECE", "Brier", "note"])
        cv = m1b["eRepo_primary"]["calibration"]
        w.writerow(
            [
                "variant",
                "DISCERN_isotonic",
                "eRepo missense P/B",
                cv["DISCERN_isotonic_oof"]["ece"],
                cv["DISCERN_isotonic_oof"]["brier"],
                "calibrated",
            ]
        )
        w.writerow(
            [
                "variant",
                "REVEL",
                "eRepo missense P/B",
                cv["REVEL_raw"]["ece"],
                cv["REVEL_raw"]["brier"],
                "raw predictor score, not gene-calibrated",
            ]
        )
        w.writerow(
            [
                "variant",
                "AlphaMissense",
                "eRepo missense P/B",
                cv["AlphaMissense_raw"]["ece"],
                cv["AlphaMissense_raw"]["brier"],
                "raw predictor score, not gene-calibrated",
            ]
        )
        cd = m3["calibration_diagnosis"]
        w.writerow(
            [
                "diagnosis",
                "DISCERN",
                "42 curated cases",
                cd["ece"],
                "",
                f"confidently-wrong at 0.8 = {cd['confidently_wrong_at_0.8']}",
            ]
        )
    # safety scenarios
    with open(
        o("benchmarks", "trustworthiness", "safety_scenarios.csv"),
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "scenario_id",
                "contraindicated_treatment",
                "hardstop_on_real_treatment",
                "hardstop_on_harmless_control",
            ]
        )
        for row in m3["safety_interlock"]["rows"]:
            w.writerow(
                [
                    row["id"],
                    row["harmful_tx"],
                    int(row["hardstop_on_real"]),
                    int(row["hardstop_on_harmless"]),
                ]
            )


def gen_diagnosis_arm():
    cases = yaml.safe_load(open(r("eval", "cases", "curated_cases.yaml"), encoding="utf-8"))[
        "cases"
    ]
    with open(
        o("benchmarks", "diagnosis_arm", "curated_cases.tsv"), "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "case_id",
                "cluster",
                "gene",
                "expected_diagnosis",
                "source_pmid",
                "present_features",
                "absent_features",
            ]
        )
        for c in cases:
            fe = c.get("features", {})
            present = ";".join(k for k, v in fe.items() if v)
            absent = ";".join(k for k, v in fe.items() if v is False)
            w.writerow(
                [
                    c["id"],
                    c["cluster"],
                    c.get("gene", ""),
                    c["true_dx"],
                    c["source_pmid"],
                    present,
                    absent,
                ]
            )
    from eval.curated_case_benchmark import run as run_curated

    s = run_curated()
    with open(
        o("benchmarks", "diagnosis_arm", "curated_case_results.csv"),
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "case_id",
                "expected",
                "leading",
                "top1_correct",
                "top3_correct",
                "decided",
                "source_pmid",
            ]
        )
        for row in s["rows"]:
            w.writerow(
                [
                    row["id"],
                    row["true"],
                    row["lead"],
                    int(row["top1"]),
                    int(row["top3"]),
                    int(row["decided"]),
                    row["pmid"],
                ]
            )
    return s


def gen_champ_chbmp():
    # No raw CDC data is redistributed. Ship the derived recall summary (from the committed benchmark)
    # + the committed PS1/REVEL subset; the full per-variant set regenerates from the pinned CDC source.
    with open(
        o("benchmarks", "champ_chbmp", "champ_chbmp_recall.csv"), "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.writer(fh)
        w.writerow(["stratum", "n_null_variants", "lp_p_recall_by_consequence_alone", "note"])
        w.writerow(
            ["overall (F8+F9 null)", 2130, 0.912, "recall by consequence alone, no predictor"]
        )
        w.writerow(["F8 null", "", 0.977, ""])
        w.writerow(
            ["F9 null", "", 0.658, "gap explained by the large terminal exon / NMD-escape rules"]
        )
        w.writerow(
            ["missense arm (novel to LP)", 1818, 0.015, "PS1+PM2+PP3 only (motivates coupling)"]
        )


def gen_coupling_poc():
    recs = [
        json.loads(x)
        for x in open(r("eval", "data", "bleeding_subset.jsonl"), encoding="utf-8")
        if x.strip()
    ]
    with open(
        o("benchmarks", "coupling_poc", "phenopacket_poc_cases.tsv"),
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["case_id", "gene", "hgvs_c", "hgvs_p", "disease", "pmid", "n_hpo_present"])
        for c in recs:
            w.writerow(
                [
                    c.get("id"),
                    c.get("gene"),
                    c.get("hgvs_c"),
                    c.get("hgvs_p"),
                    c.get("disease"),
                    c.get("pmid"),
                    len(c.get("hpo_present", [])),
                ]
            )
    with open(
        o("benchmarks", "coupling_poc", "phenopacket_poc_results.csv"),
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        w = csv.writer(fh)
        w.writerow(["endpoint", "matched", "mismatched", "lift", "interpretation"])
        w.writerow(
            [
                "primary binary upgrade to LP/P",
                "0/2",
                "0/2",
                0.0,
                "no case crossed the LP threshold; underpowered (n=2)",
            ]
        )
        w.writerow(
            [
                "secondary continuous P(path+LP)",
                0.208,
                0.082,
                0.126,
                "directionally consistent; motivates the cohort study (Gate G13)",
            ]
        )


def gen_tables(m1, m1b, curated):
    with open(
        o("tables", "table1_variant_headtohead.csv"), "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.writer(fh)
        w.writerow(["tool", "surface", "AUROC_missense", "AUPRC", "sens_at_90spec", "calibratable"])
        d1 = m1["discrimination"]
        w.writerow(
            [
                "DISCERN",
                "ClinVar missense P/B",
                d1["DISCERN_points"]["auroc"],
                d1["DISCERN_points"]["auprc"],
                d1["DISCERN_points"]["sens@90spec"],
                "yes (ECE 0.017)",
            ]
        )
        w.writerow(
            [
                "REVEL",
                "ClinVar missense P/B",
                d1["REVEL"]["auroc"],
                d1["REVEL"]["auprc"],
                d1["REVEL"]["sens@90spec"],
                "no (raw score)",
            ]
        )
        w.writerow(
            [
                "AlphaMissense",
                "ClinVar missense P/B",
                d1["AlphaMissense"]["auroc"],
                d1["AlphaMissense"]["auprc"],
                d1["AlphaMissense"]["sens@90spec"],
                "no (raw score)",
            ]
        )
        w.writerow(
            [
                "InterVar (full-DB, prior)",
                "ClinVar missense",
                d1["InterVar_full_DB_prior"]["auroc"],
                "",
                "",
                "no (class only)",
            ]
        )
        w.writerow(
            [
                "GeneBe (exhibit)",
                "ClinVar missense (circular)",
                d1["GeneBe_acmg_score"]["auroc"],
                "",
                "",
                "no (class only)",
            ]
        )
        dp = m1b["eRepo_primary"]["discrimination"]
        dt = m1b["time_split"]["discrimination"]
        w.writerow(
            [
                "DISCERN",
                "eRepo-primary missense",
                dp["DISCERN"]["auroc"],
                dp["DISCERN"]["auprc"],
                dp["DISCERN"]["sens@90spec"],
                "yes (ECE 0.017)",
            ]
        )
        w.writerow(
            [
                "DISCERN",
                "eRepo time-split missense",
                dt["DISCERN"]["auroc"],
                dt["DISCERN"]["auprc"],
                dt["DISCERN"]["sens@90spec"],
                "yes (ECE 0.019)",
            ]
        )
    cv = m1b["eRepo_primary"]["calibration"]
    with open(
        o("tables", "table2_calibration_and_partition.csv"), "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "DISCERN", "REVEL", "AlphaMissense", "note"])
        w.writerow(
            [
                "variant ECE (eRepo)",
                cv["DISCERN_isotonic_oof"]["ece"],
                cv["REVEL_raw"]["ece"],
                cv["AlphaMissense_raw"]["ece"],
                "lower is better; only DISCERN is calibrated",
            ]
        )
        w.writerow(
            [
                "variant Brier (eRepo)",
                cv["DISCERN_isotonic_oof"]["brier"],
                cv["REVEL_raw"]["brier"],
                cv["AlphaMissense_raw"]["brier"],
                "",
            ]
        )
        k = m1b["per_code_kappa_vs_erepo"]
        w.writerow(["per-code kappa PVS1", k["PVS1"]["kappa"], "", "", "DISCERN vs eRepo experts"])
        w.writerow(["per-code kappa PP3", k["PP3"]["kappa"], "", "", ""])
        w.writerow(
            [
                "PS3/PS4/PM1 applied by DISCERN",
                k.get("PS3", {}).get("discern_applied", 0),
                "",
                "",
                "0 = the partition (evidence-stream codes not derived from sequence alone)",
            ]
        )
    with open(o("tables", "table3_diagnosis_curated.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "value"])
        w.writerow(["n_cases", curated["n"]])
        w.writerow(["top1_accuracy", curated["top1"]])
        w.writerow(["top3_accuracy", curated["top3"]])
        w.writerow(["abstention_rate", curated["abstention_rate"]])


def _savefig(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(o("figures", name + "." + ext), bbox_inches="tight", dpi=300)
    plt.close(fig)


def gen_figures(rows, m1b):
    mis = [x for x in rows if x["is_missense"] and x["label"] is not None]
    y = np.array([x["label"] for x in mis])
    from bench.track1_variant_headtohead import _disc_prob_oof

    # fig1 architecture (simple schematic)
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.axis("off")
    ax.text(0.5, 0.95, "DISCERN coupled model", ha="center", fontsize=13, weight="bold")
    for i, t in enumerate(
        ["Variant (intrinsic ACMG)", "Phenotype (HPO features)", "Functional / lab"]
    ):
        ax.add_patch(plt.Rectangle((0.03, 0.62 - i * 0.22), 0.30, 0.15, fc="#e8f1f1", ec="#0c6b6b"))
        ax.text(0.18, 0.695 - i * 0.22, t, ha="center", va="center", fontsize=9)
    ax.add_patch(plt.Rectangle((0.40, 0.35), 0.22, 0.30, fc="#fdfaf3", ec="#9a5a06"))
    ax.text(0.51, 0.50, "Joint\nP(D,V|E)\npartitioned", ha="center", va="center", fontsize=9)
    for i, t in enumerate(
        ["Differential", "VUS reclass", "Safety hard-stop", "Next-best test", "Abstention"]
    ):
        ax.add_patch(plt.Rectangle((0.70, 0.80 - i * 0.17), 0.28, 0.12, fc="#f4f6f8", ec="#5a6772"))
        ax.text(0.84, 0.86 - i * 0.17, t, ha="center", va="center", fontsize=8)
    _savefig(fig, "fig1_architecture")
    # fig2 ROC + calibration
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.2))
    for name, key, col in [
        ("DISCERN", "discern_points", "#0c6b6b"),
        ("REVEL", "revel", "#b07a2b"),
        ("AlphaMissense", "alphamissense", "#5a6772"),
    ]:
        pr = [(float(x[key]), yy) for x, yy in zip(mis, y, strict=False) if x[key] is not None]
        sc = np.array([p for p, _ in pr])
        yy = np.array([b for _, b in pr])
        fpr, tpr, _ = roc_curve(yy, sc)
        a1.plot(fpr, tpr, color=col, label=name)
    a1.plot([0, 1], [0, 1], "--", color="#ccc")
    a1.set_xlabel("FPR")
    a1.set_ylabel("TPR")
    a1.set_title("Variant discrimination (eRepo missense)")
    a1.legend(fontsize=8)
    disc = _disc_prob_oof([x["discern_points"] for x in mis], y.tolist())
    for name, sc, col in [
        ("DISCERN (calibrated)", disc, "#0c6b6b"),
        ("REVEL (raw)", [x["revel"] for x in mis], "#b07a2b"),
    ]:
        pr = [(float(s), yy) for s, yy in zip(sc, y, strict=False) if s is not None]
        p = np.array([a for a, _ in pr])
        yy = np.array([b for _, b in pr])
        edges = np.linspace(0, 1, 11)
        xs = []
        ob = []
        for i in range(10):
            m = (p >= edges[i]) & (p < (edges[i + 1] if i < 9 else 1.001))
            if m.sum():
                xs.append(p[m].mean())
                ob.append(yy[m].mean())
        a2.plot(xs, ob, "o-", color=col, label=name)
    a2.plot([0, 1], [0, 1], "--", color="#ccc")
    a2.set_xlabel("mean predicted")
    a2.set_ylabel("observed")
    a2.set_title("Reliability")
    a2.legend(fontsize=8)
    _savefig(fig, "fig2_variant_roc_and_calibration")
    # fig3 per-code kappa
    k = m1b["per_code_kappa_vs_erepo"]
    codes = [c for c in k if k[c]["kappa"] is not None]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(codes, [k[c]["kappa"] for c in codes], color="#0c6b6b")
    ax.set_ylabel("Cohen kappa vs eRepo")
    ax.set_title("Per-ACMG-code agreement (the partition)")
    ax.set_ylim(0, 1)
    plt.xticks(rotation=45)
    _savefig(fig, "fig3_per_code_kappa")
    # fig4 risk-coverage
    rc = list(
        csv.DictReader(
            open(o("benchmarks", "trustworthiness", "risk_coverage_curve.csv"), encoding="utf-8")
        )
    )
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(
        [float(x["coverage"]) for x in rc],
        [float(x["accuracy_on_retained"]) for x in rc],
        "-",
        color="#0c6b6b",
    )
    ax.set_xlabel("coverage")
    ax.set_ylabel("accuracy on retained")
    ax.set_title("Risk-coverage (abstention pays off)")
    _savefig(fig, "fig4_risk_coverage")
    # fig5 genebe circularity
    m1 = _jload("track1_metrics.json")
    d = m1["discrimination"]
    tools = ["DISCERN", "REVEL", "AlphaMissense", "GeneBe\n(circular)"]
    vals = [
        d["DISCERN_points"]["auroc"],
        d["REVEL"]["auroc"],
        d["AlphaMissense"]["auroc"],
        d["GeneBe_acmg_score"]["auroc"],
    ]
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.bar(tools, vals, color=["#0c6b6b", "#b07a2b", "#5a6772", "#a51f12"])
    ax.set_ylabel("AUROC (ClinVar missense P/B)")
    ax.set_ylim(0.8, 1.02)
    ax.set_title("GeneBe reaches 1.0 by reproducing ClinVar (not a fair win)")
    _savefig(fig, "fig5_genebe_circularity")


def main():
    build_tree()
    copy_direct()
    rows, m1, m1b = gen_variant_arm()
    gen_trustworthiness()
    curated = gen_diagnosis_arm()
    gen_champ_chbmp()
    gen_coupling_poc()
    gen_tables(m1, m1b, curated)
    gen_figures(rows, m1b)
    print("supplement assembled at", OUT)
    n = sum(len(files) for _, _, files in os.walk(OUT))
    print("total files:", n)


if __name__ == "__main__":
    main()
