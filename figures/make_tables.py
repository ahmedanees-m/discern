"""Manuscript tables, written as CSV from the committed benchmark JSON.

Three main tables and eight supplementary ones, each regenerated from the same files the figures
read, so a change to a benchmark propagates to both. Table 1 keeps its capability columns: once the
calibration-uniqueness claim was retired those columns are what carry the comparison.

Run:  python -m figures.make_tables [--out DIR]
"""
from __future__ import annotations

import argparse
import csv
import json
import os

import yaml

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
BENCH = os.path.join(ROOT, "bench")
EVAL = os.path.join(ROOT, "eval")


def _load(p):
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def _write(outdir, name, header, rows, note=None):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{name}.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
        if note:
            w.writerow([])
            w.writerow([f"NOTE: {note}"])
    return path


def _ci(v):
    return "" if not v else f"[{v[0]:.3f}-{v[1]:.3f}]"


# --------------------------------------------------------------------------------------------
def table1(outdir):
    m = _load(os.path.join(BENCH, "phase_r_variant_metrics.json"))
    # GeneBe is scored by the Track 1' harness, not the Phase R one (Phase R only re-analyses the
    # three tools that emit a continuous score DISCERN can be calibrated against), so its row is
    # read from there. Without this the table's own footnote referred to a row that did not exist.
    t1b = _load(os.path.join(BENCH, "track1b_erepo_metrics.json"))
    header = ["tool", "surface", "n", "AUROC", "AUROC_95CI", "ECE", "ECE_95CI", "Brier",
              "paired_vs_DISCERN_DeLong_p", "calibratable", "emits_ACMG_class", "criterion_trail",
              "abstains", "ClinVar_blinded"]
    rows = []
    for surf, label in (("eRepo_primary", "eRepo-primary"), ("time_split", "time-split")):
        s = m[surf]
        p = s["discern_vs_revel"]["delong"]["p_value"]
        for tool, calkey, cal_ok, cls, trail, absta in (
                ("DISCERN", "DISCERN_isotonic_oof", "yes (as delivered)", "yes", "yes", "yes"),
                ("REVEL (raw)", "REVEL_raw", "n/a", "no", "no", "no"),
                ("REVEL (calibrated)", "REVEL_isotonic_oof", "yes (post-hoc, needs labels)", "no", "no", "no"),
                ("AlphaMissense (raw)", "AlphaMissense_raw", "n/a", "no", "no", "no"),
                ("AlphaMissense (calibrated)", "AlphaMissense_isotonic_oof",
                 "yes (post-hoc, needs labels)", "no", "no", "no")):
            base = tool.split(" ")[0]
            d = s["discrimination"].get(base, {})
            c = s["calibration"].get(calkey, {})
            rows.append([tool, label, d.get("n", ""), d.get("auroc", ""), _ci(d.get("auroc_ci95")),
                         c.get("ece", ""), _ci(c.get("ece_ci95")), c.get("brier", ""),
                         f"{p:.3f}" if base == "REVEL" else "", cal_ok, cls, trail, absta, "yes"])
        g = t1b.get(surf, {}).get("discrimination", {}).get("GeneBe_exhibit") or {}
        if g:
            rows.append(["GeneBe (circularity exhibit)", label, g.get("n", ""), g.get("auroc", ""),
                         "", "", "", "", "", "no (class only)", "yes", "yes",
                         "no", "NO - consumes ClinVar"])
    rows.append(["InterVar (by citation)", "prior ClinVar set", "", 0.811, "", "", "", "", "",
                 "no (class only)", "yes", "yes", "no", "yes"])
    rows.append(["BIAS-2015 (by citation)", "published eRepo metrics", "", "", "", "", "", "", "",
                 "no (class only)", "yes", "yes", "no", "yes"])
    return _write(outdir, "table1_variant_head_to_head", header, rows,
                  "InterVar and BIAS-2015 rows are published metrics on different sets and are NOT "
                  "like-for-like with the AUROC column. GeneBe is shown as a circularity exhibit, "
                  "not a fair comparator.")


def table2(outdir):
    m = _load(os.path.join(BENCH, "phase_r_variant_metrics.json"))["added_value_over_ranking_score"]
    dq = m["decision_quality"]
    header = ["rule", "coverage", "accuracy_on_resolved", "pathogenic_recall", "benign_recall",
              "n_resolved", "calls_pathogenic", "calls_benign"]
    rows = []
    for label, key in (("DISCERN ACMG class (its own operating point)", "DISCERN_acmg_class"),
                       ("DISCERN calibrated probability, coverage-matched",
                        "DISCERN_calibrated_prob_matched_coverage"),
                       ("REVEL at ClinGen thresholds", "REVEL_clingen_thresholds"),
                       ("REVEL, coverage-matched to DISCERN", "REVEL_matched_coverage")):
        d = dq.get(key) or {}
        rows.append([label, d.get("coverage", ""), d.get("accuracy_on_resolved", ""),
                     d.get("pathogenic_recall", ""), d.get("benign_recall", ""),
                     d.get("n_resolved", ""), d.get("called_pathogenic", ""),
                     d.get("called_benign", "")])
    return _write(outdir, "table2_decision_quality_matched_coverage", header, rows,
                  "The risk-coverage claim is withdrawn (the diagnosis arm is at ceiling at n=42); "
                  "this table reports coverage-matched decision quality only. Every band DISCERN "
                  "resolves on this surface is benign-side, which Table S2 and Figure 4 explain.")


def table3(outdir):
    g = _load(os.path.join(EVAL, "gene_only_baseline.json"))
    lir = _load(os.path.join(EVAL, "lirical_arm.json"))
    p, st = g["pooled"], g["by_stratum"]
    header = ["method", "n", "Top-1", "Top-3", "Top-1_unique_gene_stratum",
              "Top-1_shared_gene_stratum", "Top-1_no_gene_stratum", "McNemar_p_vs_gene_lookup",
              "in_sample"]

    def s(stratum, key, sub="top1"):
        return st.get(stratum, {}).get(key, {}).get(sub, "")

    rows = [
        ["DISCERN post-correction", p["n"], p["DISCERN"]["top1"], p["DISCERN"]["top3"],
         s("unique_gene", "DISCERN"), s("shared_gene", "DISCERN"), s("no_gene", "DISCERN"),
         p["mcnemar_top1"]["p_value_exact"], "YES - these cases exposed the correction"],
        ["DISCERN pre-correction", p["n"], 0.81, 1.0, 0.78, 1.0, 1.0, "", "no"],
        ["gene lookup (phenotype-blind)", p["n"], p["gene_lookup"]["top1"], p["gene_lookup"]["top3"],
         s("unique_gene", "gene_lookup"), s("shared_gene", "gene_lookup"),
         s("no_gene", "gene_lookup"), "", "no"],
        ["prior only (gene-blind)", p["n"], p["prior_only"]["top1"], "",
         s("unique_gene", "prior_only"), s("shared_gene", "prior_only"),
         s("no_gene", "prior_only"), "", "no"],
        ["uniform random within cluster", p["n"],
         p["random_within_cluster"]["top1"]["expected"],
         p["random_within_cluster"]["top3"]["expected"], "", "", "", "", "no"],
    ]
    for label, key in (("DISCERN phenotype-only (no gene), vs LIRICAL subset",
                        "DISCERN_phenotype_only_no_gene"),
                       ("DISCERN HPO-representable only (no gene), vs LIRICAL subset",
                        "DISCERN_hpo_representable_only"),
                       ("LIRICAL restricted to the same cluster", "LIRICAL_restricted_to_cluster")):
        d = lir[key]
        rows.append([label, d["n"], d["recall@1"], d["recall@3"], "", "", "", "", "no"])
    sizes = g["strata_sizes"]
    return _write(outdir, "table3_diagnosis_vs_baselines", header, rows,
                  f"Strata partition the benchmark exactly: unique_gene {sizes.get('unique_gene')} "
                  f"+ shared_gene {sizes.get('shared_gene')} + no_gene {sizes.get('no_gene')} = "
                  f"{sum(sizes.values())}. The post-correction row is in-sample and must not be "
                  f"quoted against an external tool.")


# --------------------------------------------------------------------------------------------
def tableS1(outdir):
    m = _load(os.path.join(BENCH, "track1b_erepo_metrics.json"))["per_code_kappa_vs_erepo"]
    routed = {"PS3", "BS3", "PP4", "PP1", "BS4", "PM3", "BP2", "PS2", "PM6"}
    rows = []
    for code, v in m.items():
        if v["discern_applied"] == 0:
            cause = ("routed to another factor by the partition" if code in routed
                     else "variant-intrinsic; no input in this pipeline")
        else:
            cause = "applied"
        rows.append([code, v["erepo_applied"], v["discern_applied"],
                     "" if v["kappa"] is None else v["kappa"], cause])
    return _write(outdir, "tableS1_per_criterion_kappa",
                  ["ACMG_criterion", "erepo_applied", "discern_applied", "kappa",
                   "zero_application_cause"], rows)


def tableS2(outdir):
    ca = _load(os.path.join(BENCH, "phase_r_variant_metrics.json"))["ceiling_attribution"]
    rows = [["as scored (intrinsic evidence only)", ca["reach_lp_as_scored"], ca["n"]],
            ["+ intrinsic criteria with no input here",
             ca["reach_lp_if_intrinsic_codes_were_available"], ca["n"]],
            ["+ criteria the partition routes away",
             ca["reach_lp_if_routed_codes_were_re_added"], ca["n"]],
            ["both restored", ca["reach_lp_with_both"], ca["n"]]]
    rows.append([])
    rows.append(["code", "reason unapplied", "kind of limit"])
    for code, why in ca["why_each_intrinsic_code_has_no_input"].items():
        kind = ("engine scope gap" if code == "PM1"
                else "data availability" if code == "PS4" else "protocol choice (ClinVar-blinded)")
        rows.append([code, why, kind])
    return _write(outdir, "tableS2_ceiling_attribution",
                  ["restoration_condition", "variants_reaching_LP", "of_n"], rows, ca["caveat"])


def tableS3(outdir):
    rows = []
    cdir = os.path.join(ROOT, "diseases", "clusters")
    for fn in sorted(os.listdir(cdir)):
        if not fn.endswith(".yaml"):
            continue
        with open(os.path.join(cdir, fn), encoding="utf-8") as fh:
            y = yaml.safe_load(fh)
        for d in y.get("diseases", []):
            for feat, lr in (d.get("feature_lr") or {}).items():
                val, n, pmid = (list(lr) + [None, None, None])[:3]
                rows.append([y["id"], y["name"], d["id"], d["name"], ",".join(d.get("genes", [])),
                             feat, val, n, pmid])
    return _write(outdir, "tableS3_cluster_likelihood_ratios",
                  ["cluster_id", "cluster_name", "disease_id", "disease_name", "genes",
                   "feature", "P(feature|disease)", "sample_size", "source_pmid"], rows,
                  "Every likelihood ratio carries a source PMID; a CI guard fails the build if one "
                  "is missing.")


def tableS4(outdir):
    s = _load(os.path.join(BENCH, "track3_metrics.json"))["safety_interlock"]
    rows = [[r["id"], r["harmful_tx"], "fire", "fired" if r["hardstop_on_real"] else "SILENT",
             "stay silent", "fired" if r["hardstop_on_harmless"] else "silent"] for r in s["rows"]]
    return _write(outdir, "tableS4_safety_scenarios",
                  ["scenario", "planned_management", "expected_on_contraindicated",
                   "observed_on_contraindicated", "expected_on_harmless", "observed_on_harmless"],
                  rows, f"sensitivity {s['hardstop_sensitivity']:.0%}, specificity "
                        f"{s['hardstop_specificity']:.0%} over {s['n_scenarios']} scenarios; "
                        f"adjudicated gene-blind against the engine's own leading call.")


def tableS5(outdir):
    from eval.phenotype_tool_comparison import case_hpo, feature_to_hpo
    cases = yaml.safe_load(open(os.path.join(EVAL, "cases", "curated_cases.yaml"),
                                encoding="utf-8"))["cases"]
    g = {r["id"]: r for r in _load(os.path.join(EVAL, "gene_only_baseline.json"))["per_case"]}
    f2h = feature_to_hpo()
    rows = []
    for c in cases:
        r = g.get(c["id"], {})
        rows.append([c["id"], c.get("source_pmid", ""), c.get("gene", ""), c["cluster"],
                     c["true_dx"], ";".join(case_hpo(c, f2h)), r.get("stratum", ""),
                     r.get("discern_lead", ""), r.get("discern_top1", ""),
                     r.get("gene_lookup_lead", ""), r.get("gene_lookup_top1", "")])
    return _write(outdir, "tableS5_curated_cases",
                  ["case_id", "source_pmid", "gene", "cluster", "expected_diagnosis",
                   "extracted_hpo_terms", "stratum", "discern_lead", "discern_top1",
                   "gene_lookup_lead", "gene_lookup_top1"], rows,
                  "PMIDs and extracted HPO terms only. No article text is reproduced (Gate G7).")


def tableS6(outdir):
    from figures.make_figures import CHAMP
    rows = [["F8", CHAMP["F8"], "NMD-triggering nulls resolved by PVS1"],
            ["F9", CHAMP["F9"], "gap is the large terminal exon: last-exon PTCs escape NMD and are "
                                "correctly held at VUS"],
            ["overall", CHAMP["overall"], f"{CHAMP['n_null']} null variants of "
                                          f"{CHAMP['n_alleles']} disease alleles"]]
    return _write(outdir, "tableS6_champ_chbmp_recall",
                  ["gene", "LP_P_recall_on_null", "note"], rows,
                  "Recall by consequence alone, with no frequency or predictor input. Source: "
                  "docs/DISCERN_CHAMP_CHBMP_Benchmark_v1.md")


def tableS7(outdir):
    """Data source manifest - the sources that are referenced, not redistributed."""
    rows = [
        ["ClinGen Evidence Repository (eRepo)", "erepo_classifications.tab",
         "https://erepo.clinicalgenome.org/evrepo/", "2026 snapshot", "expert-panel truth surface"],
        ["ClinVar", "clinvar_20260503.vcf.gz + variant_summary.txt.gz",
         "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/", "20260503", "supporting surface; blinded in the primary protocol"],
        ["gnomAD", "gnomad.v4.0.constraint_metrics.tsv",
         "https://gnomad.broadinstitute.org/downloads", "v4.0", "population frequency"],
        ["Human Phenotype Ontology", "hp.json, phenotype.hpoa, genes_to_phenotype.txt",
         "https://hpo.jax.org/app/data/ontology", "2026 release", "HPO terms and gene-disease annotation"],
        ["REVEL", "via ANNOVAR dbnsfp / GeneBe annotation", "https://sites.google.com/site/revelgenomics/",
         "as bundled", "PP3/BP4 missense signal"],
        ["AlphaMissense", "AlphaMissense_hg38.tsv.gz",
         "https://zenodo.org/records/8360242", "2023", "alternative missense predictor"],
        ["CHAMP / CHBMP", "CHAMP-Variant-List-2022.xlsx, CHBMP-Variant-List-2022.xlsx",
         "https://www.cdc.gov/hemophilia/php/mutation-project/", "2022", "independent F8/F9 truth"],
        ["GA4GH Phenopacket Store", "bleeding subset", "https://github.com/monarch-initiative/phenopacket-store",
         "as cited", "coupling proof-of-concept corpus"],
        ["LIRICAL", "lirical-cli-2.4.1 + its data bundle",
         "https://github.com/TheJacksonLaboratory/LIRICAL", "2.4.1", "external phenotype-driven comparator"],
    ]
    return _write(outdir, "tableS7_data_source_manifest",
                  ["source", "file(s) used", "url", "version_or_snapshot", "role"], rows,
                  "None of these is redistributed in this deposit. Checksums for the files actually "
                  "used are in MANIFEST.md; re-pull at the pinned versions to reproduce.")


def tableS8(outdir):
    m = _load(os.path.join(EVAL, "lirical_arm.json"))
    rows = []
    for key in ("LIRICAL_genome_wide", "LIRICAL_restricted_to_cluster",
                "DISCERN_hpo_representable_only", "DISCERN_phenotype_only_no_gene",
                "DISCERN_pre_gene_term_fix", "DISCERN_post_fix_IN_SAMPLE"):
        d = m[key]
        rows.append([key, d["n"], d["recall@1"], d["recall@3"], d["recall@5"], d["mrr"]])
    rows.append([])
    rows.append(["paired test", "arm", "delta", "CI95_low", "CI95_high", "McNemar_p"])
    for name, t in m["paired_tests_vs_lirical_restricted"].items():
        rows.append(["vs LIRICAL restricted", name, t["delta"], t["delta_ci95"][0],
                     t["delta_ci95"][1], t["mcnemar"]["p_value_exact"]])
    rows.append([])
    rows.append(["phenotype-only arm is invariant to the gene-term fix",
                 m["phenotype_only_arm_is_fix_invariant"]])
    return _write(outdir, "tableS8_lirical_arms",
                  ["arm", "n", "recall@1", "recall@3", "recall@5", "MRR"], rows,
                  "LIRICAL 2.4.1, phenotype-only mode, observed and negated HPO terms supplied, run "
                  "in a container. Quote the hpo_representable_only arm for inference claims and "
                  "the phenotype_only arm for the architectural claim.")


def tableS9(outdir):
    """Software and tool versions. Reproducibility reviewers check for this specifically."""
    import platform
    import subprocess

    def _v(mod):
        try:
            return __import__(mod).__version__
        except Exception:                                            # noqa: BLE001
            return "not installed in this environment"

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                            text=True).stdout.strip() or "unknown"
    rows = [
        ["DISCERN", commit, "this work", "https://github.com/ahmedanees-m/discern"],
        ["Python", platform.python_version(), "runtime", ""],
        ["numpy", _v("numpy"), "arithmetic and resampling", ""],
        ["scipy", _v("scipy"), "binomial test, normal quantiles", ""],
        ["scikit-learn", _v("sklearn"), "isotonic calibration, AUROC, Brier, kappa", ""],
        ["matplotlib", _v("matplotlib"), "figure rendering", ""],
        ["PyYAML", _v("yaml"), "cluster and specification loading", ""],
        ["ruff", "0.15.12 (pinned in CI)", "lint", ""],
        ["LIRICAL", "2.4.1", "external phenotype-driven comparator",
         "https://github.com/TheJacksonLaboratory/LIRICAL"],
        ["eclipse-temurin JRE", "17", "container runtime for LIRICAL", ""],
        ["REVEL", "as bundled in the GeneBe annotation cache", "PP3/BP4 missense signal", ""],
        ["AlphaMissense", "hg38 release, Zenodo record 8360242", "alternative missense predictor",
         "https://zenodo.org/records/8360242"],
        ["gnomAD", "v4.0 constraint metrics; v2.1.1/v4 allele frequencies", "population frequency",
         "https://gnomad.broadinstitute.org"],
        ["ClinGen eRepo", "2026 snapshot", "expert-panel truth surface",
         "https://erepo.clinicalgenome.org"],
        ["ClinVar", "20260503", "supporting surface, blinded in the primary protocol",
         "https://www.ncbi.nlm.nih.gov/clinvar"],
        ["Human Phenotype Ontology", "2026 release", "HPO terms and gene-disease annotation",
         "https://hpo.jax.org"],
        ["ANNOVAR", "as used for the prior InterVar comparison", "annotation for the H4 arm", ""],
        ["InterVar", "prior full-database run", "ACMG classifier comparator", ""],
        ["CHAMP / CHBMP", "2022 variant lists", "independent F8/F9 truth",
         "https://www.cdc.gov/hemophilia-data/champ-mutation-project/"],
    ]
    return _write(outdir, "tableS9_software_versions",
                  ["component", "version", "role in this work", "resource"], rows,
                  f"Generated on the machine that produced the deposit. DISCERN's own commit is the "
                  f"authoritative version: {commit}. Python package versions record the environment "
                  f"the committed results were produced in; continuous integration pins ruff and "
                  f"installs the remainder from pyproject.toml.")


TABLES = [table1, table2, table3, tableS1, tableS2, tableS3, tableS4, tableS5, tableS6, tableS7,
          tableS8, tableS9]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "tables"))
    args = ap.parse_args()
    for fn in TABLES:
        print(f"  {fn.__name__:12} -> {os.path.basename(fn(args.out))}")
    print(f"\nwrote {len(TABLES)} tables to {args.out}")


if __name__ == "__main__":
    main()
