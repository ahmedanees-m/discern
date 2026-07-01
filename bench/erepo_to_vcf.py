"""Build a GRCh38 VCF from the eRepo bleeding SNV set, for Nirvana/BIAS-2015 (Phase A2 Part A).

Reads bench/data/erepo_bleeding.tsv and writes a sorted, deduplicated, chr-prefixed VCFv4.2 that
Nirvana (3.18.1) annotates. SNVs only (already filtered upstream). Output: bench/data/erepo_bleeding.vcf
"""
from __future__ import annotations

import csv
import os

HERE = os.path.dirname(__file__)
SRC = os.path.join(HERE, "data", "erepo_bleeding.tsv")
OUT = os.path.join(HERE, "data", "erepo_bleeding.vcf")

_ORD = {str(c): c for c in range(1, 23)}
_ORD["X"], _ORD["Y"] = 23, 24


def _key(r):
    return (_ORD.get(r["chrom"], 99), int(r["pos"]), r["ref"], r["alt"])


def run():
    rows = list(csv.DictReader(open(SRC, encoding="utf-8"), delimiter="\t"))
    seen, uniq = set(), []
    for r in rows:
        k = (r["chrom"], r["pos"], r["ref"], r["alt"])
        if k in seen:
            continue
        seen.add(k)
        uniq.append(r)
    uniq.sort(key=_key)
    with open(OUT, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in uniq:
            chrom = r["chrom"] if r["chrom"].startswith("chr") else "chr" + r["chrom"]
            fh.write(f"{chrom}\t{r['pos']}\t.\t{r['ref']}\t{r['alt']}\t.\t.\t.\n")
    return len(uniq)


if __name__ == "__main__":
    n = run()
    print(f"wrote {os.path.relpath(OUT)} with {n} variants")
