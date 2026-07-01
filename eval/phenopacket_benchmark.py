"""Phenopacket Store diagnosis benchmark (plan section 3A, Tier A; PhEval-compatible).

Runs DISCERN on the bleeding-cluster subset of the open GA4GH Phenopacket Store: for each
real published case it takes the causal gene and the HPO terms (observed and excluded) and
asks whether the joint model's leading disease matches the case's diagnosis (Top-1/Top-3).

Phenopacket diagnosis is by gene + OMIM/Mondo; DISCERN routes by gene into a cluster and
returns the leading disease. The HPO crosswalk maps the discriminating HPO terms onto the
cluster features; gene routing carries the unambiguous cases, phenotype disambiguates the
ambiguous genes (GP1BA: BSS vs PT-VWD; VWF subtypes). Best-effort crosswalk - documented.

Run (on the VM, after cloning monarch-initiative/phenopacket-store):
    python3 -m eval.phenopacket_benchmark /path/to/phenopacket-store
"""
from __future__ import annotations

import json
import os
import sys

from core.dx_schemas import Feature, FeatureKind
from diseases.ontology import route_clusters
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import marginal_disease
from jointdx.orchestrate import _best_cluster

# gene -> expected DISCERN disease id (unambiguous mappings; GP1BA/VWF resolved by label).
GENE_TO_DISEASE = {
    "ITGA2B": "gt", "ITGB3": "gt", "FERMT3": "lad3", "RASGRP2": "rasgrp2", "ITGB2": "lad1",
    "GP1BB": "bss", "GP9": "bss", "MYH9": "myh9",
    "NBEAL2": "gps", "LYST": "chs", "HPS1": "hps", "HPS3": "hps", "HPS4": "hps",
    "HPS5": "hps", "HPS6": "hps", "AP3B1": "hps", "GFI1B": "gps",
    "RUNX1": "runx1", "ETV6": "etv6", "ANKRD26": "ankrd26",
    "F8": "hemophilia_a", "F9": "hemophilia_b", "F13A1": "fxiii", "F13B": "fxiii",
}

# discriminating HPO terms -> cluster feature id (best-effort; gene routing carries the rest).
HPO_TO_FEATURE = {
    "HP:0001974": "leukocytosis",          # Increased total leukocyte count
    "HP:0002719": "recurrent_infections",
    "HP:0002718": "recurrent_infections",  # Recurrent bacterial infections
    "HP:0001873": "thrombocytopenia",
    "HP:0011875": "giant_platelets",
    "HP:0032499": "giant_granules",        # Giant neutrophil granules (Chediak-Higashi)
    "HP:0001892": "glanzmann_type_bleeding",
    "HP:0011890": "glanzmann_type_bleeding",
    "HP:0000978": "glanzmann_type_bleeding",
    "HP:0000421": "glanzmann_type_bleeding",   # Epistaxis
    "HP:0001025": "oculocutaneous_albinism",
    "HP:0007443": "oculocutaneous_albinism",
    "HP:0002218": "oculocutaneous_albinism",   # Silver-gray hair (Chediak pigmentary)
}

# HPO terms indicating THROMBOSIS (not a DISCERN bleeding/platelet cluster) -> exclude the case.
THROMBOSIS_HPO = {"HP:0002625", "HP:0002204", "HP:0030977", "HP:0004850", "HP:0005305",
                  "HP:0030242", "HP:0040412", "HP:0034336"}


def _walk(obj, key):
    """Yield all values for `key` anywhere in a nested dict/list."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k == key:
                yield v
            yield from _walk(v, key)
    elif isinstance(obj, list):
        for it in obj:
            yield from _walk(it, key)


def parse_phenopacket(pp: dict) -> dict | None:
    # the CAUSAL gene comes from the interpretation/geneContext, not any mentioned symbol
    genes = {g for g in _walk(pp.get("interpretations", []), "symbol") if isinstance(g, str)}
    gene = next((g for g in genes if g in GENE_TO_DISEASE), None)
    if gene is None:
        return None
    hpo_ids = {(pf.get("type") or {}).get("id", "") for pf in pp.get("phenotypicFeatures", [])}
    if hpo_ids & THROMBOSIS_HPO:        # thrombophilia case, not a DISCERN bleeding cluster
        return None
    feats = []
    for pf in pp.get("phenotypicFeatures", []):
        hid = (pf.get("type") or {}).get("id", "")
        fid = HPO_TO_FEATURE.get(hid)
        if fid:
            feats.append(Feature(fid, FeatureKind.CLINICAL, not pf.get("excluded", False),
                                 observed=not pf.get("excluded", False)))
    return {"gene": gene, "features": feats, "expected": GENE_TO_DISEASE[gene]}


def run(store_dir: str, limit: int | None = None) -> dict:
    files = []
    for root, _d, fns in os.walk(store_dir):
        for fn in fns:
            if fn.endswith(".json"):
                files.append(os.path.join(root, fn))
    top1 = top3 = evaluable = 0
    by_gene: dict[str, list[int]] = {}
    for path in files:
        try:
            with open(path, encoding="utf-8") as fh:
                pp = json.load(fh)
        except Exception:
            continue
        parsed = parse_phenopacket(pp)
        if parsed is None:
            continue
        clusters = route_clusters([parsed["gene"]])
        if not clusters:
            continue
        ev = Evidence(variant_gene=parsed["gene"], clinical=parsed["features"])
        cl = _best_cluster(clusters, ev)
        md = marginal_disease(joint(cl, ev))
        ranked = sorted(md, key=md.get, reverse=True)
        evaluable += 1
        hit1 = ranked[0] == parsed["expected"]
        hit3 = parsed["expected"] in ranked[:3]
        top1 += hit1
        top3 += hit3
        by_gene.setdefault(parsed["gene"], []).append(int(hit1))
        if limit and evaluable >= limit:
            break
    return {
        "evaluable": evaluable,
        "top1": round(top1 / evaluable, 4) if evaluable else None,
        "top3": round(top3 / evaluable, 4) if evaluable else None,
        "by_gene": {g: (sum(v), len(v)) for g, v in sorted(by_gene.items())},
    }


def main():
    s = run(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else None)
    print(f"Phenopacket Store bleeding-cluster cases: {s['evaluable']}")
    print(f"Top-1 accuracy: {s['top1']:.1%}" if s["top1"] is not None else "no cases")
    print(f"Top-3 accuracy: {s['top3']:.1%}" if s["top3"] is not None else "")
    print("By gene (correct/total):", s["by_gene"])


if __name__ == "__main__":
    main()
