"""Curated published-case diagnosis benchmark (v3.1 Track B4).

Runs DISCERN over the hand-curated, citation-only case set (`eval/cases/curated_cases.yaml`)
and reports Top-1 / Top-3 diagnosis accuracy + abstention. Small-n BY DESIGN (the cohorts
carry the headline; Gate B/C). No patient identifiers are used.

Run:  python3 -m eval.curated_case_benchmark
"""
from __future__ import annotations

import os

import yaml

from core.dx_schemas import Feature, FeatureKind
from diseases.ontology import cluster_for
from jointdx.abstain import decide
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import marginal_disease

_CASES = os.path.join(os.path.dirname(__file__), "cases", "curated_cases.yaml")


def load_cases(path: str = _CASES) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)["cases"]


def _evidence(case: dict) -> Evidence:
    # Codebase convention (see tests + evidence/phenotype_lr.py): `observed` is the present/absent
    # flag the phenotype LR reads. A feature set false is a PERTINENT NEGATIVE (observed=False),
    # NOT a feature scored as present. (Earlier this hardcoded observed=True, which inverted every
    # pertinent negative - masked by the prior 10-case set, exposed by the expanded set.)
    clin = [Feature(fid, FeatureKind.LAB, bool(present), observed=bool(present))
            for fid, present in case.get("features", {}).items()]
    return Evidence(variant_gene=case.get("gene", ""), clinical=clin)


def run(path: str = _CASES) -> dict:
    cases = load_cases(path)
    top1 = top3 = decided = 0
    rows = []
    for c in cases:
        cluster = cluster_for(c["cluster"])
        ev = _evidence(c)
        md = marginal_disease(joint(cluster, ev))
        ranked = sorted(md.items(), key=lambda kv: kv[1], reverse=True)
        order = [d for d, _ in ranked]
        true = c["true_dx"]
        is1 = order[:1] == [true]
        is3 = true in order[:3]
        d = decide(cluster, ev, n_mc=40)
        top1 += is1
        top3 += is3
        decided += d.decided
        rows.append({"id": c["id"], "true": true, "top1": is1, "top3": is3,
                     "lead": order[0], "decided": d.decided, "pmid": c.get("source_pmid")})
    n = len(cases)
    return {"n": n, "top1": round(top1 / n, 4), "top3": round(top3 / n, 4),
            "abstention_rate": round(1 - decided / n, 4), "rows": rows}


def main():
    s = run()
    print(f"Curated cases: n={s['n']}  Top-1={s['top1']:.0%}  Top-3={s['top3']:.0%}  "
          f"abstention={s['abstention_rate']:.0%}")
    for r in s["rows"]:
        mark = "OK" if r["top1"] else ("top3" if r["top3"] else "MISS")
        print(f"  [{mark:4}] {r['id']:22} true={r['true']:10} lead={r['lead']:10} PMID:{r['pmid']}")


if __name__ == "__main__":
    main()
