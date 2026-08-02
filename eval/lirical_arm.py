"""PHASE R - the LIRICAL head-to-head on the curated cases.

Runs the comparison the pre-submission plan asked for rather than declining it, and is explicit
about what the number can support. Two things make this comparison awkward, and both are reported
instead of smoothed over:

  1. Only 13 of the benchmark's 48 discriminating features have any HPO representation, so 19 of
     the 42 cases carry no phenotype term at all and LIRICAL has nothing to rank them on. It is run
     on the 23 cases it can accept, with DISCERN restricted to the identical subset.
  2. LIRICAL ranks diseases genome-wide by OMIM/ORPHA identifier; DISCERN ranks within a cluster of
     three to eight look-alikes. Scoring is therefore done at gene/disease-family resolution - a
     LIRICAL hit counts if its top-k contains any disease identifier annotated to the case's causal
     gene. That is the generous reading for LIRICAL: it cannot be penalised for OMIM granularity,
     but it also means the score measures gene prioritisation rather than within-cluster
     discrimination, which is the thing DISCERN exists to do.

The gene-to-disease crosswalk is derived from the committed HPO annotation file rather than
asserted, so every identifier traces to a public source.

Run:  python -m eval.lirical_arm inputs        # write the per-case HPO input TSV
      python -m eval.lirical_arm score FILE    # score LIRICAL's collected output
"""
from __future__ import annotations

import csv
import json
import os
import sys
from collections import defaultdict

from diseases.ontology import cluster_for
from eval.curated_case_benchmark import load_cases
from eval.phenotype_tool_comparison import case_hpo, discern_ranking, feature_to_hpo

HERE = os.path.dirname(__file__)
G2P = os.path.join(HERE, "..", "data", "raw", "hpo", "genes_to_phenotype.txt")
INPUT_TSV = os.path.join(HERE, "data", "lirical_cases.tsv")
OUT_JSON = os.path.join(HERE, "lirical_arm.json")


def negated_hpo(case, f2h) -> list[str]:
    """Pertinent negatives, so LIRICAL gets the same absent findings DISCERN uses."""
    terms: set[str] = set()
    for fid, present in (case.get("features") or {}).items():
        if not present and fid in f2h:
            terms |= set(f2h[fid])
    return sorted(terms)


def gene_to_disease_ids() -> dict[str, set[str]]:
    """Gene symbol -> the OMIM/ORPHA disease identifiers HPO annotates to it."""
    out: dict[str, set[str]] = defaultdict(set)
    with open(G2P, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("disease_id"):
                out[row["gene_symbol"]].add(row["disease_id"])
    return out


def disease_ids_for(disease_id: str, cluster_id: str, g2d) -> set[str]:
    cluster = cluster_for(cluster_id)
    d = next((x for x in cluster.diseases if x.id == disease_id), None)
    if d is None:
        return set()
    ids: set[str] = set()
    for g in d.genes:
        ids |= g2d.get(g, set())
    return ids


def write_inputs() -> int:
    f2h = feature_to_hpo()
    os.makedirs(os.path.dirname(INPUT_TSV), exist_ok=True)
    n = 0
    with open(INPUT_TSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["case_id", "observed_hpo", "negated_hpo"])
        for c in load_cases():
            obs, neg = case_hpo(c, f2h), negated_hpo(c, f2h)
            if not obs:
                continue                       # LIRICAL cannot rank a case with no phenotype term
            w.writerow([c["id"], ",".join(obs), ",".join(neg)])
            n += 1
    return n


def _metrics(hit_at):
    n = len(hit_at)
    if not n:
        return {"n": 0}
    def r(k):
        return round(sum(1 for h in hit_at if h is not None and h <= k) / n, 4)
    mrr = sum(1.0 / h for h in hit_at if h is not None) / n
    return {"n": n, "recall@1": r(1), "recall@3": r(3), "recall@5": r(5),
            "recall@10": r(10), "mrr": round(mrr, 4)}


def score(path: str) -> dict:
    """Score LIRICAL output: TSV of case_id <TAB> ranked disease identifiers, comma separated."""
    cases = {c["id"]: c for c in load_cases()}
    g2d = gene_to_disease_ids()

    lirical_hits, lirical_cluster_hits = [], []
    discern_hits, discern_prefix_hits, discern_nogene_hits, rows = [], [], [], []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("case_id"):
                continue
            cid, _, order = line.partition("\t")
            c = cases.get(cid)
            if c is None:
                continue
            ranked = [x.strip() for x in order.split(",") if x.strip()]
            truth_ids = disease_ids_for(c["true_dx"], c["cluster"], g2d)
            hit = next((i + 1 for i, d in enumerate(ranked) if d in truth_ids), None)
            lirical_hits.append(hit)

            # the like-for-like arm: collapse LIRICAL's genome-wide list onto the same cluster
            # DISCERN is choosing within, so both tools rank the same three to eight diseases
            members = {d.id: disease_ids_for(d.id, c["cluster"], g2d)
                       for d in cluster_for(c["cluster"]).diseases}
            seen, cluster_rank = set(), []
            for did in ranked:
                for mid, ids in members.items():
                    if did in ids and mid not in seen:
                        seen.add(mid)
                        cluster_rank.append(mid)
            c_hit = cluster_rank.index(c["true_dx"]) + 1 if c["true_dx"] in cluster_rank else None
            lirical_cluster_hits.append(c_hit)

            d_rank = discern_ranking(c)
            d_hit = d_rank.index(c["true_dx"]) + 1 if c["true_dx"] in d_rank else None
            discern_hits.append(d_hit)

            # the two arms that make this reportable (see the module docstring)
            pre = discern_ranking(c, gene_term=False)
            discern_prefix_hits.append(pre.index(c["true_dx"]) + 1 if c["true_dx"] in pre else None)
            ng = discern_ranking(c, drop_gene=True)
            discern_nogene_hits.append(ng.index(c["true_dx"]) + 1 if c["true_dx"] in ng else None)

            rows.append({"id": cid, "gene": c.get("gene", ""), "true_dx": c["true_dx"],
                         "lirical_rank": hit, "lirical_within_cluster_rank": c_hit,
                         "discern_rank": d_hit, "lirical_top1": ranked[0] if ranked else None,
                         "n_truth_ids": len(truth_ids)})

    return {
        "scoring_resolution": ("gene / disease-family: a LIRICAL hit counts if its top-k contains "
                               "any HPO-annotated disease identifier for the case's causal gene"),
        "headline_arm": "DISCERN_phenotype_only_no_gene",
        "why": ("These 42 cases are what exposed the missing P(G|D) term, so post-fix DISCERN has "
                "seen them and LIRICAL has not. Two arms escape that problem. The phenotype-only "
                "arm withholds the gene entirely, which makes P(G|D) inert - it is therefore "
                "identical before and after the fix, uncontaminated by construction, and it is also "
                "the only arm whose inputs match what LIRICAL receives. It is the headline. The "
                "pre-fix arm reproduces the engine as it stood before these cases informed it, and "
                "is reported as a second uncontaminated reference. The post-fix arm is labelled "
                "in-sample and must never be quoted as a head-to-head result."),
        "phenotype_only_arm_is_fix_invariant": (
            _metrics(discern_nogene_hits) == _metrics(
                [discern_ranking(cases[r["id"]], gene_term=False, drop_gene=True).index(
                    cases[r["id"]]["true_dx"]) + 1
                 if cases[r["id"]]["true_dx"] in discern_ranking(
                     cases[r["id"]], gene_term=False, drop_gene=True) else None
                 for r in rows])),
        "LIRICAL_genome_wide": _metrics(lirical_hits),
        "LIRICAL_restricted_to_cluster": _metrics(lirical_cluster_hits),
        "DISCERN_pre_gene_term_fix": _metrics(discern_prefix_hits),
        "DISCERN_phenotype_only_no_gene": _metrics(discern_nogene_hits),
        "DISCERN_post_fix_IN_SAMPLE": _metrics(discern_hits),
        "input_parity": {
            "LIRICAL_receives": "HPO terms only (observed and negated); no gene",
            "DISCERN_phenotype_only_receives": "the same findings, no gene - the matched-input arm",
            "DISCERN_pre_and_post_fix_receive": "the same findings plus the causal gene",
            "note": ("The gene is worth most of this benchmark: a phenotype-blind gene lookup "
                     "scores 93 percent on the full 42. Any arm that receives the gene is therefore "
                     "not comparable to LIRICAL, which does not. Read "
                     "DISCERN_phenotype_only_no_gene against LIRICAL_restricted_to_cluster."),
        },
        "caveat": ("The genome-wide row is not a like-for-like contest: LIRICAL ranks roughly 8,600 "
                   "diseases from HPO terms, while DISCERN ranks the three to eight members of one "
                   "cluster. The restricted row is the fair ranking surface. Even there the tools "
                   "optimise different objectives, and the reason LIRICAL does poorly is not that "
                   "it is a weak ranker but that the findings which separate these diseases are "
                   "laboratory assays the Human Phenotype Ontology does not encode. That is a fact "
                   "about the domain, not a win."),
        "rows": rows,
    }


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "inputs":
        n = write_inputs()
        print(f"wrote {n} runnable cases to {os.path.relpath(INPUT_TSV)}")
        return
    if len(sys.argv) > 2 and sys.argv[1] == "score":
        o = score(sys.argv[2])
        with open(OUT_JSON, "w", encoding="utf-8") as fh:
            json.dump(o, fh, indent=2)
        print("== PHASE R: LIRICAL (phenotype-only) vs DISCERN, identical case subset ==")
        print(f"   scoring at {o['scoring_resolution']}")
        for k in ("LIRICAL_genome_wide", "LIRICAL_restricted_to_cluster",
                  "DISCERN_phenotype_only_no_gene", "DISCERN_pre_gene_term_fix",
                  "DISCERN_post_fix_IN_SAMPLE"):
            m = o[k]
            print(f"   {k:32} n={m['n']}  R@1={m['recall@1']:.0%}  R@3={m['recall@3']:.0%}  "
                  f"R@5={m['recall@5']:.0%}  R@10={m['recall@10']:.0%}  MRR={m['mrr']:.3f}")
        print(f"\n   {o['caveat']}")
        print(f"   wrote {os.path.relpath(OUT_JSON)}")
        return
    print(__doc__)


if __name__ == "__main__":
    main()
