"""Tests for the diagnosis-arm baselines and the HPO-coverage measurement.

These two harnesses produce results the manuscript leans on heavily: the phenotype-blind gene
lookup that the headline Top-1 must be quoted against, the stratification that decides which cases
are informative, and the finding that only a minority of discriminating features have any HPO
representation. Neither needs third-party data, so unlike the variant-arm harnesses they can and
should run in continuous integration.
"""
from __future__ import annotations

import numpy as np
import pytest

from diseases.ontology import all_clusters, cluster_for
from eval.curated_case_benchmark import load_cases
from eval.gene_only_baseline import (
    _mcnemar,
    _random_expectation,
    _stratum,
    _topk,
    gene_lookup_rank,
    prior_only_rank,
)
from eval.phenotype_tool_comparison import (
    case_hpo,
    coverage,
    discern_ranking,
    feature_to_hpo,
)


@pytest.fixture(scope="module")
def cases():
    return load_cases()


# ---------------------------------------------------------------- gene lookup
def test_gene_lookup_puts_the_matching_disease_first(cases):
    """The baseline must be a fair opponent: it should exploit the gene fully."""
    for c in cases:
        gene = (c.get("gene") or "").strip()
        if not gene:
            continue
        order = gene_lookup_rank(c)
        cluster = cluster_for(c["cluster"])
        matching = {d.id for d in cluster.diseases if gene in d.genes}
        if matching:
            assert order[0] in matching, c["id"]


def test_gene_lookup_ranks_every_disease_in_the_cluster_exactly_once(cases):
    for c in cases:
        order = gene_lookup_rank(c)
        ids = [d.id for d in cluster_for(c["cluster"]).diseases]
        assert sorted(order) == sorted(ids)


def test_gene_lookup_uses_no_phenotype(cases):
    """Stripping every finding must not change the baseline's ranking."""
    for c in cases[:12]:
        stripped = dict(c, features={})
        assert gene_lookup_rank(c) == gene_lookup_rank(stripped)


def test_prior_only_ignores_the_gene(cases):
    for c in cases[:12]:
        assert prior_only_rank(c) == prior_only_rank(dict(c, gene=""))


def test_prior_only_is_sorted_by_descending_prior(cases):
    order = prior_only_rank(cases[0])
    priors = {d.id: d.prior for d in cluster_for(cases[0]["cluster"]).diseases}
    assert [priors[i] for i in order] == sorted(priors.values(), reverse=True)


# ---------------------------------------------------------------- stratification
def test_strata_partition_the_benchmark_exactly(cases):
    """36 + 3 + 3 = 42; a reviewer will add these up (consistency check C1)."""
    counts = {}
    for c in cases:
        counts[_stratum(c)] = counts.get(_stratum(c), 0) + 1
    assert sum(counts.values()) == len(cases) == 42
    assert set(counts) <= {"no_gene", "shared_gene", "unique_gene", "gene_outside_cluster"}
    assert counts.get("gene_outside_cluster", 0) == 0


def test_stratum_labels_match_the_cluster_definition(cases):
    for c in cases:
        gene = (c.get("gene") or "").strip()
        s = _stratum(c)
        if not gene:
            assert s == "no_gene"
            continue
        n = sum(1 for d in cluster_for(c["cluster"]).diseases if gene in d.genes)
        assert s == ("shared_gene" if n > 1 else "unique_gene")


def test_shared_gene_stratum_is_the_informative_one(cases):
    """By construction these are the cases the gene alone cannot settle."""
    shared = [c for c in cases if _stratum(c) == "shared_gene"]
    for c in shared:
        gene = c["gene"].strip()
        matching = [d.id for d in cluster_for(c["cluster"]).diseases if gene in d.genes]
        assert len(matching) > 1


# ---------------------------------------------------------------- helpers
def test_topk_is_monotone_and_exact():
    order = ["a", "b", "c", "d"]
    assert _topk(order, "a", 1) == 1
    assert _topk(order, "c", 1) == 0
    assert _topk(order, "c", 3) == 1
    assert _topk(order, "z", 4) == 0


def test_random_floor_matches_the_analytic_expectation(cases):
    rng = np.random.default_rng(0)
    r1 = _random_expectation(cases, 1, rng, n_boot=200)
    sizes = [len(cluster_for(c["cluster"]).diseases) for c in cases]
    assert r1["expected"] == pytest.approx(float(np.mean([1 / s for s in sizes])), abs=1e-4)
    assert r1["sim_ci95"][0] <= r1["expected"] <= r1["sim_ci95"][1]


def test_random_floor_rises_with_k(cases):
    rng = np.random.default_rng(0)
    a = _random_expectation(cases, 1, rng, n_boot=100)["expected"]
    b = _random_expectation(cases, 3, rng, n_boot=100)["expected"]
    assert b > a


def test_mcnemar_wrapper_keeps_the_reported_key_names():
    m = _mcnemar([1, 1, 0], [0, 0, 1])
    assert set(m) == {"discern_only_correct", "baseline_only_correct", "discordant_pairs",
                      "p_value_exact"}
    assert m["discern_only_correct"] == 2
    assert m["baseline_only_correct"] == 1


# ---------------------------------------------------------------- HPO coverage
def test_crosswalk_inverts_without_losing_terms():
    f2h = feature_to_hpo()
    assert f2h
    for feature, terms in f2h.items():
        assert isinstance(feature, str) and terms
        assert all(t.startswith("HP:") for t in terms)


def test_hpo_coverage_reports_the_finding_the_manuscript_states():
    c = coverage()
    assert c["distinct_features_used"] == 48
    assert c["features_expressible_as_hpo"] == 13
    assert c["hpo_terms_per_case"]["zero"] == 19
    assert c["runnable_cases"] == 23
    assert c["runnable_cases"] + c["hpo_terms_per_case"]["zero"] == c["n_cases"] == 42


def test_unmappable_features_are_the_laboratory_assays():
    """The motivating negative: what HPO cannot express is the deciding lab work."""
    unmapped = set(coverage()["features_with_no_hpo_term"])
    for assay in ("flow_cd42_reduced", "loss_hmw_multimers", "normal_lta",
                  "abnormal_prothrombinase_assay", "aiib3_expression_absent"):
        assert assay in unmapped, assay


def test_case_hpo_only_returns_present_findings(cases):
    f2h = feature_to_hpo()
    for c in cases:
        absent_only = {k: False for k in (c.get("features") or {})}
        assert case_hpo(dict(c, features=absent_only), f2h) == []


# ---------------------------------------------------------------- ranking arms
def test_phenotype_only_arm_is_invariant_to_the_gene_term(cases):
    """The property that makes the LIRICAL comparison reportable: withholding the gene makes
    P(G|D) inert, so this arm is identical before and after the Phase R correction."""
    for c in cases:
        with_term = discern_ranking(c, gene_term=True, drop_gene=True)
        without = discern_ranking(c, gene_term=False, drop_gene=True)
        assert with_term == without, c["id"]


def test_hpo_restricted_arm_sees_only_codable_features(cases):
    """The evidence-matched arm must discard exactly the findings LIRICAL cannot receive."""
    f2h = feature_to_hpo()
    codable = set(f2h)
    for c in cases[:15]:
        stripped = dict(c, features={k: v for k, v in (c.get("features") or {}).items()
                                     if k in codable})
        assert discern_ranking(c, drop_gene=True, hpo_representable_only=True) == \
            discern_ranking(stripped, drop_gene=True)


def test_every_ranking_arm_covers_the_whole_cluster(cases):
    for c in cases[:10]:
        n = len(cluster_for(c["cluster"]).diseases)
        for kwargs in ({}, {"gene_term": False}, {"drop_gene": True},
                       {"drop_gene": True, "hpo_representable_only": True}):
            assert len(discern_ranking(c, **kwargs)) == n


def test_clusters_and_cases_match_the_reported_coverage(cases):
    cl = all_clusters()
    assert len(cl) == 10
    assert sum(len(x.diseases) for x in cl.values()) == 31
    assert len({g for x in cl.values() for d in x.diseases for g in d.genes}) == 29
    assert {c["cluster"] for c in cases} <= set(cl)
