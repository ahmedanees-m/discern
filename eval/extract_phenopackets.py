"""Extract the bleeding/platelet-gene subset from a local Phenopacket Store clone (Coupling PoC).

Parses GA4GH Phenopacket v2 JSON (camelCase) and writes one record per case to JSONL plus an
HPO-term histogram (to drive the crosswalk data-first). Public literature cases (PHI-free);
source = github.com/monarch-initiative/phenopacket-store (Danis et al., HGG Adv 2024, PMID 39394689).

Usage: python3 extract_phenopackets.py <phenopacket-store-root> <out.jsonl> <hpo_hist.tsv>
"""
from __future__ import annotations

import json
import pathlib
import sys
from collections import Counter

# ISTH TIER1 bleeding/platelet gene set (broad; cluster-modelled subset is selected downstream).
GENES = {
    "ITGA2B", "ITGB3", "FERMT3", "RASGRP2", "ITGB2", "GP1BA", "GP1BB", "GP9", "MYH9", "VWF",
    "F8", "F9", "F13A1", "F13B", "F2", "F5", "F7", "F10", "F11", "FGA", "FGB", "FGG",
    "RUNX1", "ANKRD26", "ETV6", "NBEAL2", "GFI1B", "PLAU", "VPS33B", "LYST",
    "HPS1", "HPS3", "HPS4", "HPS5", "HPS6", "AP3B1", "ANO6", "ACTN1", "TUBB1", "WAS",
    "GATA1", "FLI1", "FLNA", "STIM1", "ORAI1", "P2RY12", "GP6", "ITGA2", "MPL", "THPO", "SLFN14",
}


def _first(d, *keys, default=None):
    for k in keys:
        if isinstance(d, dict) and k in d:
            d = d[k]
        else:
            return default
    return d


def parse(pp: dict) -> dict | None:
    gene = hgvs_c = hgvs_g = hgvs_p = disease = None
    for interp in pp.get("interpretations", []):
        diag = interp.get("diagnosis", {})
        disease = disease or _first(diag, "disease", "label")
        for gi in diag.get("genomicInterpretations", []):
            vd = _first(gi, "variantInterpretation", "variationDescriptor") or {}
            gc = vd.get("geneContext", {})
            if gc.get("symbol"):
                gene = gc["symbol"]
            for ex in vd.get("expressions", []):
                syn, val = ex.get("syntax", ""), ex.get("value", "")
                if syn == "hgvs.c" and not hgvs_c:
                    hgvs_c = val
                elif syn == "hgvs.g" and not hgvs_g:
                    hgvs_g = val
                elif syn == "hgvs.p" and not hgvs_p:
                    hgvs_p = val
                elif syn == "hgvs" and ":p." in val and not hgvs_p:
                    hgvs_p = val
    if gene is None or gene not in GENES:
        return None
    present, excluded = [], []
    for pf in pp.get("phenotypicFeatures", []):
        t = pf.get("type", {})
        item = {"id": t.get("id", ""), "label": t.get("label", "")}
        (excluded if pf.get("excluded") else present).append(item)
    if not disease:
        for d in pp.get("diseases", []):
            disease = disease or _first(d, "term", "label")
    pmid = ""
    for ref in _first(pp, "metaData", "externalReferences", default=[]) or []:
        if str(ref.get("id", "")).upper().startswith("PMID"):
            pmid = ref["id"]
            break
    return {"id": pp.get("id", ""), "gene": gene, "hgvs_c": hgvs_c, "hgvs_g": hgvs_g,
            "hgvs_p": hgvs_p, "disease": disease, "pmid": pmid,
            "hpo_present": present, "hpo_excluded": excluded}


def main(root: str, out_jsonl: str, hist_tsv: str):
    rootp = pathlib.Path(root)
    hist = Counter()
    n = 0
    with open(out_jsonl, "w", encoding="utf-8") as out:
        for jf in rootp.rglob("*.json"):
            try:
                pp = json.loads(jf.read_text(encoding="utf-8", errors="replace"))
            except Exception:
                continue
            if "phenotypicFeatures" not in pp and "interpretations" not in pp:
                continue
            rec = parse(pp)
            if rec is None:
                continue
            out.write(json.dumps(rec) + "\n")
            n += 1
            for h in rec["hpo_present"]:
                hist[(h["id"], h["label"])] += 1
    with open(hist_tsv, "w", encoding="utf-8") as fh:
        fh.write("hpo_id\tlabel\tcount\n")
        for (hid, lab), c in hist.most_common():
            fh.write(f"{hid}\t{lab}\t{c}\n")
    print(f"bleeding-gene cases: {n}; distinct HPO terms: {len(hist)}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
