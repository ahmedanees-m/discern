"""PHASE R - what the phenotype channel actually adds to the diagnosis arm.

The curated benchmark hands the engine the causal gene. For several clusters the gene names the
disease on its own, so a headline Top-1 over all 42 cases is confounded by the cases where there was
never a question to answer. This module measures the confound instead of assuming it away.

Three reference points, none of which sees a phenotype:
  gene lookup   - rank the cluster with the diseases that list the case gene first, ordered by their
                  prior. A static table, no likelihood ratios, no findings.
  prior only    - rank the cluster by prior and ignore the gene. Separates "the gene helps" from
                  "one disease dominates this cluster".
  random        - draw uniformly within the cluster. The floor.

and the split that decides whether the phenotype engine earns its place:
  shared-gene cases   - the gene maps to more than one disease inside its cluster, so something other
                        than the gene has to break the tie. This subset is the real result.
  unique-gene cases   - the gene maps to exactly one disease in the cluster. The gene is the answer.

Paired case-level comparison against DISCERN uses McNemar (exact) plus a bootstrap CI on the delta.

Run:  python -m eval.gene_only_baseline
"""
from __future__ import annotations

import json
import os

import numpy as np
from scipy.stats import binomtest

from diseases.ontology import cluster_for
from eval.curated_case_benchmark import load_cases
from eval.curated_case_benchmark import run as discern_run

HERE = os.path.dirname(__file__)
OUT_JSON = os.path.join(HERE, "gene_only_baseline.json")
N_BOOT = 1000
SEED = 0


def _ranked_by_prior(diseases):
    return [d.id for d in sorted(diseases, key=lambda d: (-d.prior, d.id))]


def gene_lookup_rank(case) -> list[str]:
    """Cluster ranking from the gene alone: gene-matching diseases first, each block by prior."""
    cluster = cluster_for(case["cluster"])
    gene = (case.get("gene") or "").strip()
    match = [d for d in cluster.diseases if gene and gene in d.genes]
    rest = [d for d in cluster.diseases if not (gene and gene in d.genes)]
    return _ranked_by_prior(match) + _ranked_by_prior(rest)


def prior_only_rank(case) -> list[str]:
    return _ranked_by_prior(cluster_for(case["cluster"]).diseases)


def _stratum(case) -> str:
    gene = (case.get("gene") or "").strip()
    if not gene:
        return "no_gene"
    cluster = cluster_for(case["cluster"])
    n = sum(1 for d in cluster.diseases if gene in d.genes)
    return "shared_gene" if n > 1 else ("unique_gene" if n == 1 else "gene_outside_cluster")


def _topk(order, true_dx, k):
    return int(true_dx in order[:k])


def _random_expectation(cases, k, rng, n_boot=N_BOOT):
    """Uniform-within-cluster floor: exact expectation plus a simulated spread."""
    sizes = [len(cluster_for(c["cluster"]).diseases) for c in cases]
    exact = float(np.mean([min(k, s) / s for s in sizes]))
    sims = []
    for _ in range(n_boot):
        hits = 0
        for s in sizes:
            hits += int(rng.integers(0, s) < min(k, s))
        sims.append(hits / len(sizes))
    return {"expected": round(exact, 4),
            "sim_ci95": [round(float(np.percentile(sims, 2.5)), 4),
                         round(float(np.percentile(sims, 97.5)), 4)]}


def _mcnemar(a_correct, b_correct):
    """Exact McNemar on paired case-level correctness (a = DISCERN, b = baseline)."""
    b_only = sum(1 for a, bb in zip(a_correct, b_correct, strict=True) if a and not bb)
    c_only = sum(1 for a, bb in zip(a_correct, b_correct, strict=True) if bb and not a)
    n = b_only + c_only
    p = float(binomtest(b_only, n, 0.5).pvalue) if n else 1.0
    return {"discern_only_correct": b_only, "baseline_only_correct": c_only,
            "discordant_pairs": n, "p_value_exact": round(p, 6)}


def _delta_ci(a_correct, b_correct, rng, n_boot=N_BOOT):
    a, b = np.asarray(a_correct, int), np.asarray(b_correct, int)
    idx = rng.integers(0, len(a), size=(n_boot, len(a)))
    d = a[idx].mean(axis=1) - b[idx].mean(axis=1)
    return {"delta": round(float(a.mean() - b.mean()), 4),
            "delta_ci95": [round(float(np.percentile(d, 2.5)), 4),
                           round(float(np.percentile(d, 97.5)), 4)]}


def _summary(cases, discern_rows, rng, label):
    idx = {r["id"]: r for r in discern_rows}
    disc1 = [int(idx[c["id"]]["top1"]) for c in cases]
    disc3 = [int(idx[c["id"]]["top3"]) for c in cases]
    gl1 = [_topk(gene_lookup_rank(c), c["true_dx"], 1) for c in cases]
    gl3 = [_topk(gene_lookup_rank(c), c["true_dx"], 3) for c in cases]
    pr1 = [_topk(prior_only_rank(c), c["true_dx"], 1) for c in cases]
    n = len(cases)
    return {
        "stratum": label, "n": n,
        "DISCERN": {"top1": round(float(np.mean(disc1)), 4), "top3": round(float(np.mean(disc3)), 4)},
        "gene_lookup": {"top1": round(float(np.mean(gl1)), 4), "top3": round(float(np.mean(gl3)), 4)},
        "prior_only": {"top1": round(float(np.mean(pr1)), 4)},
        "random_within_cluster": {"top1": _random_expectation(cases, 1, rng),
                                  "top3": _random_expectation(cases, 3, rng)},
        "discern_minus_gene_lookup_top1": _delta_ci(disc1, gl1, rng),
        "mcnemar_top1": _mcnemar(disc1, gl1),
    }


def run() -> dict:
    cases = load_cases()
    discern_rows = discern_run()["rows"]
    rng = np.random.default_rng(SEED)

    strata = {}
    for c in cases:
        strata.setdefault(_stratum(c), []).append(c)

    out = {
        "n_cases": len(cases),
        "strata_sizes": {k: len(v) for k, v in sorted(strata.items())},
        "pooled": _summary(cases, discern_rows, rng, "pooled"),
        "by_stratum": {k: _summary(v, discern_rows, rng, k) for k, v in sorted(strata.items())},
        "per_case": [
            {"id": c["id"], "cluster": c["cluster"], "gene": c.get("gene", ""),
             "stratum": _stratum(c), "true_dx": c["true_dx"],
             "discern_lead": next(r["lead"] for r in discern_rows if r["id"] == c["id"]),
             "gene_lookup_lead": gene_lookup_rank(c)[0],
             "discern_top1": int(next(r["top1"] for r in discern_rows if r["id"] == c["id"])),
             "gene_lookup_top1": _topk(gene_lookup_rank(c), c["true_dx"], 1)}
            for c in cases
        ],
    }
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    return out


def main():
    o = run()
    print("== PHASE R: does the phenotype channel beat knowing the gene? ==")
    print(f"n={o['n_cases']}  strata: {o['strata_sizes']}")
    for key in ["pooled", *sorted(o["by_stratum"])]:
        s = o["pooled"] if key == "pooled" else o["by_stratum"][key]
        rnd = s["random_within_cluster"]
        print(f"\n-- {key} (n={s['n']}) --")
        print(f"   DISCERN       Top-1={s['DISCERN']['top1']:.0%}  Top-3={s['DISCERN']['top3']:.0%}")
        print(f"   gene lookup   Top-1={s['gene_lookup']['top1']:.0%}  Top-3={s['gene_lookup']['top3']:.0%}")
        print(f"   prior only    Top-1={s['prior_only']['top1']:.0%}")
        print(f"   random floor  Top-1={rnd['top1']['expected']:.0%}  Top-3={rnd['top3']['expected']:.0%}")
        d = s["discern_minus_gene_lookup_top1"]
        m = s["mcnemar_top1"]
        print(f"   delta (DISCERN - gene lookup) = {d['delta']:+.0%}  95% CI [{d['delta_ci95'][0]:+.0%}, {d['delta_ci95'][1]:+.0%}]")
        print(f"   McNemar: DISCERN-only right={m['discern_only_correct']}  baseline-only right="
              f"{m['baseline_only_correct']}  discordant={m['discordant_pairs']}  p={m['p_value_exact']}")
    print(f"\nwrote {os.path.relpath(OUT_JSON)}")


if __name__ == "__main__":
    main()
