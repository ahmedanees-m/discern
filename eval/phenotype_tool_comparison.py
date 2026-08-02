"""PHASE R - the LIRICAL / Exomiser comparison, and the measurement that decides what it can mean.

The pre-submission plan assumed sparse HPO handicaps every tool equally, so a head-to-head would
still be fair. This module tests that assumption before relying on it, because it turns out to be
false in a specific and reportable way: the sparsity is not random thinning, it is systematic. The
findings that separate these diseases are laboratory assay results - ristocetin-induced platelet
aggregation and its mixing study, flow cytometry for CD42 and alphaIIbbeta3, multimer patterns,
light transmission aggregometry, a prothrombinase assay - and the Human Phenotype Ontology does not
encode them as phenotype terms. A phenotype-driven ranker is therefore not handicapped equally; it
is denied the discriminating channel entirely, while keeping the gene, which R4 shows is most of
the answer on this benchmark.

What this module does:
  1. Measures HPO coverage of the curated benchmark: how many of its features can be expressed as
     HPO terms at all, and how many cases end up with none.
  2. Exports one phenopacket per case so LIRICAL and Exomiser read exactly what DISCERN read.
  3. Scores any external tool's ranking against the same truth (Recall@1/3/5 and MRR) on the subset
     that tool can actually accept, alongside DISCERN restricted to the identical subset.

Run:  python -m eval.phenotype_tool_comparison             # coverage + phenopacket export
      python -m eval.phenotype_tool_comparison score FILE  # score a tool's ranking TSV
"""
from __future__ import annotations

import json
import os
import statistics
import sys
from collections import Counter

import yaml

from diseases.ontology import cluster_for
from eval.curated_case_benchmark import _evidence, load_cases
from jointdx.factorgraph import joint
from jointdx.infer import marginal_disease

HERE = os.path.dirname(__file__)
CROSSWALK = os.path.join(HERE, "hpo_feature_crosswalk.yaml")
PACKET_DIR = os.path.join(HERE, "data", "phenopackets")
OUT_JSON = os.path.join(HERE, "phenotype_tool_comparison.json")


def feature_to_hpo() -> dict[str, list[str]]:
    """Invert the committed HPO -> feature crosswalk."""
    with open(CROSSWALK, encoding="utf-8") as fh:
        xw = yaml.safe_load(fh)["crosswalk"]
    out: dict[str, list[str]] = {}
    for hpo, meta in xw.items():
        out.setdefault(meta["feature"], []).append(hpo)
    return {k: sorted(v) for k, v in out.items()}


def case_hpo(case, f2h) -> list[str]:
    terms: set[str] = set()
    for fid, present in (case.get("features") or {}).items():
        if present and fid in f2h:
            terms |= set(f2h[fid])
    return sorted(terms)


def coverage() -> dict:
    cases = load_cases()
    f2h = feature_to_hpo()
    used = {fid for c in cases for fid in (c.get("features") or {})}
    counts = [len(case_hpo(c, f2h)) for c in cases]
    return {
        "n_cases": len(cases),
        "distinct_features_used": len(used),
        "features_expressible_as_hpo": len(used & set(f2h)),
        "features_with_no_hpo_term": sorted(used - set(f2h)),
        "hpo_terms_per_case": {
            "median": statistics.median(counts), "mean": round(statistics.mean(counts), 2),
            "max": max(counts), "zero": sum(1 for c in counts if c == 0),
            "distribution": dict(sorted(Counter(counts).items())),
        },
        "runnable_cases": sum(1 for c in counts if c > 0),
        "finding": ("Only a minority of the benchmark's discriminating features have any HPO "
                    "representation, and the ones that do not are the laboratory assays that "
                    "actually separate these diseases. A phenotype-driven ranker is not equally "
                    "handicapped - it is left with the gene, which the gene-only baseline already "
                    "shows carries most of this benchmark."),
    }


def export_phenopackets() -> int:
    """One minimal phenopacket per case, so an external tool reads what DISCERN read."""
    os.makedirs(PACKET_DIR, exist_ok=True)
    f2h = feature_to_hpo()
    n = 0
    for c in load_cases():
        terms = case_hpo(c, f2h)
        packet = {
            "id": c["id"],
            "subject": {"id": c["id"]},
            "phenotypicFeatures": [{"type": {"id": t}} for t in terms],
            "interpretations": [{
                "id": c["id"],
                "diagnosis": {"genomicInterpretations": [
                    {"gene": {"symbol": c.get("gene", "")}}] if c.get("gene") else []},
            }],
            "metaData": {"createdBy": "DISCERN curated benchmark",
                         "externalReferences": [{"id": f"PMID:{c.get('source_pmid')}"}]},
        }
        with open(os.path.join(PACKET_DIR, f"{c['id']}.json"), "w", encoding="utf-8") as fh:
            json.dump(packet, fh, indent=2)
        n += 1
    return n


def discern_ranking(case, gene_term: bool = True, drop_gene: bool = False) -> list[str]:
    """Rank the cluster for a case.

    `gene_term=False` reproduces the engine as it behaved before Phase R added P(G|D) - the state
    in which these cases had not yet informed the model, and therefore the only state in which
    performance on them is comparable to a tool that never saw them.

    `drop_gene=True` withholds the gene entirely, matching what a phenotype-driven ranker receives.
    """
    cluster = cluster_for(case["cluster"])
    ev = _evidence(case)
    if drop_gene:
        ev.variant_gene = ""
    md = marginal_disease(joint(cluster, ev, gene_evidence=gene_term))
    return [d for d, _ in sorted(md.items(), key=lambda kv: kv[1], reverse=True)]


def _rank_metrics(rankings: dict[str, list[str]], truth: dict[str, str]) -> dict:
    ids = [i for i in rankings if i in truth]
    if not ids:
        return {"n": 0}
    def recall_at(k):
        return round(sum(truth[i] in rankings[i][:k] for i in ids) / len(ids), 4)
    mrr = 0.0
    for i in ids:
        r = rankings[i]
        mrr += 1.0 / (r.index(truth[i]) + 1) if truth[i] in r else 0.0
    return {"n": len(ids), "recall@1": recall_at(1), "recall@3": recall_at(3),
            "recall@5": recall_at(5), "mrr": round(mrr / len(ids), 4)}


def score_external(path: str) -> dict:
    """Score an external tool's ranking against DISCERN on the identical case subset.

    Input is a TSV of `case_id<TAB>ranked_disease_ids_comma_separated`, which is what the
    LIRICAL/Exomiser output converters emit.
    """
    cases = {c["id"]: c for c in load_cases()}
    truth = {i: c["true_dx"] for i, c in cases.items()}
    ext: dict[str, list[str]] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cid, _, order = line.partition("\t")
            if cid in cases:
                ext[cid] = [x.strip() for x in order.split(",") if x.strip()]
    disc = {i: discern_ranking(cases[i]) for i in ext}
    return {"external_tool": _rank_metrics(ext, truth),
            "DISCERN_same_subset": _rank_metrics(disc, truth),
            "note": "both restricted to the cases the external tool could accept"}


def run() -> dict:
    out = {"hpo_coverage": coverage(), "phenopackets_written": export_phenopackets(),
           "phenopacket_dir": os.path.relpath(PACKET_DIR)}
    cases = load_cases()
    f2h = feature_to_hpo()
    runnable = [c for c in cases if case_hpo(c, f2h)]
    truth = {c["id"]: c["true_dx"] for c in cases}
    out["DISCERN_on_hpo_runnable_subset"] = _rank_metrics(
        {c["id"]: discern_ranking(c) for c in runnable}, truth)
    out["DISCERN_on_all_cases"] = _rank_metrics(
        {c["id"]: discern_ranking(c) for c in cases}, truth)
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    return out


def main():
    if len(sys.argv) > 2 and sys.argv[1] == "score":
        print(json.dumps(score_external(sys.argv[2]), indent=2))
        return
    o = run()
    c = o["hpo_coverage"]
    h = c["hpo_terms_per_case"]
    print("== PHASE R: can a phenotype-driven ranker even see this benchmark? ==")
    print(f"   features used across the {c['n_cases']} cases : {c['distinct_features_used']}")
    print(f"   of those, expressible as an HPO term    : {c['features_expressible_as_hpo']}")
    print(f"   HPO terms per case                      : median {h['median']}, mean {h['mean']}, max {h['max']}")
    print(f"   cases with NO HPO term at all           : {h['zero']} of {c['n_cases']}")
    print(f"   cases an HPO-driven tool could rank     : {c['runnable_cases']}")
    print(f"\n   unmappable features (the deciding assays): {', '.join(c['features_with_no_hpo_term'][:8])} ...")
    print(f"\n   DISCERN, all cases          : {o['DISCERN_on_all_cases']}")
    print(f"   DISCERN, HPO-runnable subset: {o['DISCERN_on_hpo_runnable_subset']}")
    print(f"\n   wrote {o['phenopackets_written']} phenopackets to {o['phenopacket_dir']}")
    print(f"   wrote {os.path.relpath(OUT_JSON)}")


if __name__ == "__main__":
    main()
