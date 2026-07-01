"""Build the ANNOVAR/InterVar variant set from eRepo x ClinVar VCF (v3.1 full-H4 prep).

Extracts the GRCh38 records for the DISCERN bleeding-cluster eRepo variants from the ClinVar
GRCh38 VCF (ID column = ClinVar VariationID) and writes:
  - <out>.vcf       : the subset VCF (for convert2annovar / InterVar)
  - <out>.meta.tsv  : vid, chrom, pos, ref, alt, gene, clnsig (binary P/B label source)

Run on the VM:
    python3 -m eval.build_variant_set data/raw/erepo/erepo_classifications.tab \
        data/raw/clinvar/clinvar_20260503.vcf.gz /tmp/h4set
"""
from __future__ import annotations

import csv
import gzip
import re
import sys

CLUSTER_GENES = {
    "ITGA2B", "ITGB3", "FERMT3", "RASGRP2", "ITGB2", "VWF", "GP1BA", "GP1BB", "GP9",
    "MYH9", "ACTN1", "TUBB1", "FLNA", "NBEAL2", "GFI1B", "LYST", "HPS1", "HPS3", "HPS4",
    "HPS5", "HPS6", "AP3B1", "RUNX1", "ETV6", "ANKRD26", "F8", "F9", "F11", "F7", "F10",
    "F5", "F13A1", "F13B", "FGA", "FGB", "FGG", "F2", "ANO6", "PLAU", "VPS33B",
}
_CLNSIG = re.compile(r"CLNSIG=([^;]+)")
_GENEINFO = re.compile(r"GENEINFO=([^:;]+)")


def run(erepo, vcf, out):
    want = {}
    with open(erepo, encoding="utf-8", errors="replace") as fh:
        r = csv.reader(fh, delimiter="\t")
        h = next(r)
        gi, vi = h.index("HGNC Gene Symbol"), h.index("ClinVar Variation Id")
        for row in r:
            if len(row) <= max(gi, vi):
                continue
            g, vid = row[gi].strip(), row[vi].strip()
            if g in CLUSTER_GENES and vid:
                want[vid] = g

    n = 0
    opener = gzip.open if vcf.endswith(".gz") else open
    with opener(vcf, "rt", encoding="utf-8", errors="replace") as fh, \
            open(out + ".vcf", "w", encoding="utf-8") as ov, \
            open(out + ".meta.tsv", "w", encoding="utf-8") as om:
        om.write("vid\tchrom\tpos\tref\talt\tgene\tclnsig\n")
        for line in fh:
            if line.startswith("#"):
                ov.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            vid = f[2]
            if vid not in want:
                continue
            chrom, pos, ref, alt, info = f[0], f[1], f[3], f[4], f[7]
            if "," in alt:                       # skip multiallelic for a clean 1-1 join
                continue
            clnsig = (_CLNSIG.search(info).group(1) if _CLNSIG.search(info) else "")
            gene = want[vid]
            ov.write(line)
            om.write(f"{vid}\t{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{clnsig}\n")
            n += 1
    print(f"cluster variants wanted: {len(want)}; written to subset (single-allele): {n}")


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], sys.argv[3])
