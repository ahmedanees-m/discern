"""Extract the eRepo bleeding-gene benchmark set for Phase A2 Part A (eRepo-primary surface).

Reads the committed ClinGen eRepo export (data/raw/erepo/erepo_classifications.tab; FDA-recognized
expert-panel classifications), restricts to the DISCERN cluster gene panel, pulls the GRCh38
coordinate straight out of the HGVS Expressions column (the NC_0000* g. notation - no external
converter needed), and records the expert assertion (label), the Applied Evidence Codes (for the
per-code kappa), and the Approval/Published dates (for the time-split panel). SNVs only (the missense
headline axis is all SNVs); indels are counted and reported, not silently dropped.

Output: bench/data/erepo_bleeding.tsv  (the variant list GeneBe then annotates).
Run:  python -m bench.erepo_extract
"""
from __future__ import annotations

import csv
import glob
import os
import re

import yaml

HERE = os.path.dirname(__file__)
TAB = os.path.join(HERE, "..", "data", "raw", "erepo", "erepo_classifications.tab")
OUT = os.path.join(HERE, "data", "erepo_bleeding.tsv")

# GRCh38 RefSeq chromosome accessions -> chrom
NC38 = {
    "NC_000001.11": "1", "NC_000002.12": "2", "NC_000003.12": "3", "NC_000004.12": "4",
    "NC_000005.10": "5", "NC_000006.12": "6", "NC_000007.14": "7", "NC_000008.11": "8",
    "NC_000009.12": "9", "NC_000010.11": "10", "NC_000011.10": "11", "NC_000012.12": "12",
    "NC_000013.11": "13", "NC_000014.9": "14", "NC_000015.10": "15", "NC_000016.10": "16",
    "NC_000017.11": "17", "NC_000018.10": "18", "NC_000019.10": "19", "NC_000020.11": "20",
    "NC_000021.9": "21", "NC_000022.11": "22", "NC_000023.11": "X", "NC_000024.10": "Y",
}
SNV_RE = re.compile(r"(NC_0000\d{2}\.\d+):g\.(\d+)([ACGT]+)>([ACGT]+)")

ASSERT = {"Pathogenic": "P", "Likely Pathogenic": "LP", "Uncertain Significance": "VUS",
          "Likely Benign": "LB", "Benign": "B"}


def cluster_genes() -> set:
    genes = set()
    for f in glob.glob(os.path.join(HERE, "..", "diseases", "clusters", "*.yaml")):
        d = yaml.safe_load(open(f, encoding="utf-8"))
        for dis in d.get("diseases", []):
            genes.update(dis.get("genes", []))
    return genes


def _grch38_snv(hgvs_expressions: str):
    for m in SNV_RE.finditer(hgvs_expressions or ""):
        acc, pos, ref, alt = m.groups()
        if acc in NC38 and len(ref) == 1 and len(alt) == 1:
            return NC38[acc], pos, ref, alt
    return None


def run():
    genes = cluster_genes()
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    rows = list(csv.DictReader(open(TAB, encoding="utf-8"), delimiter="\t"))
    sub = [r for r in rows if r.get("HGNC Gene Symbol") in genes]
    out_rows, n_snv, n_indel = [], 0, 0
    for r in sub:
        lab = ASSERT.get((r.get("Assertion") or "").strip())
        if lab is None:
            continue
        snv = _grch38_snv(r.get("HGVS Expressions", ""))
        if snv is None:
            n_indel += 1
            continue
        n_snv += 1
        chrom, pos, ref, alt = snv
        out_rows.append({
            "gene": r["HGNC Gene Symbol"], "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "assertion": lab, "approval_date": r.get("Approval Date", ""),
            "published_date": r.get("Published Date", ""),
            "codes_met": (r.get("Applied Evidence Codes (Met)") or "").replace("\t", " "),
            "clinvar_vid": r.get("ClinVar Variation Id", ""),
        })
    cols = ["gene", "chrom", "pos", "ref", "alt", "assertion", "approval_date", "published_date",
            "codes_met", "clinvar_vid"]
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)
    return {"erepo_total": len(rows), "cluster_gene_rows": len(sub), "snv_extracted": n_snv,
            "indel_or_other_skipped": n_indel, "genes": len(genes)}


def main():
    import collections
    s = run()
    print(f"eRepo total={s['erepo_total']}  cluster-gene rows={s['cluster_gene_rows']}  "
          f"SNVs extracted={s['snv_extracted']}  indels/other skipped={s['indel_or_other_skipped']}")
    rows = list(csv.DictReader(open(OUT, encoding="utf-8"), delimiter="\t"))
    print("assertion dist:", collections.Counter(r["assertion"] for r in rows).most_common())
    yrs = collections.Counter((r["approval_date"] or "")[:4] for r in rows)
    print("approval-year dist:", sorted(yrs.items()))
    print(f"wrote {os.path.relpath(OUT)}")


if __name__ == "__main__":
    main()
