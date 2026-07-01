"""Reconcile the FULL-DB multianno to InterVar's column contract for the clean full-database run.

Unlike the prior reduced-config `fixup_multianno.py` (which stubbed avsnp147/1000g2015aug_all/
AAChange.knownGene because those DBs were unavailable), this run has them all present and real, so
NOTHING is stubbed. The only reconciliation is the frequency-DB column name: gnomAD v2.1.1 EXOME is
used (the coding-appropriate reference for an all-coding-variant benchmark; ~7x more coding samples
than the genome callset) and its AF columns are mapped into the gnomad_genome slot InterVar reads.
CLNSIG->CLINSIG and the phyloP alias are cosmetic (InterVar does not use them for ACMG criteria).
Position-independent (InterVar finds columns by name). v3.1 clean full-DB InterVar run.
"""
import sys

RENAME = {
    "AF": "gnomAD_genome_ALL", "AF_afr": "gnomAD_genome_AFR", "AF_amr": "gnomAD_genome_AMR",
    "AF_eas": "gnomAD_genome_EAS", "AF_nfe": "gnomAD_genome_NFE", "AF_fin": "gnomAD_genome_FIN",
    "AF_asj": "gnomAD_genome_ASJ", "AF_oth": "gnomAD_genome_OTH", "CLNSIG": "CLINSIG",
    "phyloP100way_vertebrate": "phyloP46way_placental",
}

inp, out = sys.argv[1], sys.argv[2]
with open(inp, encoding="utf-8", errors="replace") as f, open(out, "w", encoding="utf-8") as o:
    hdr = [RENAME.get(c, c) for c in f.readline().rstrip("\n").split("\t")]
    o.write("\t".join(hdr) + "\n")
    n = 0
    for line in f:
        o.write(line)
        n += 1
print(f"fixed full-DB multianno: {n} rows, {len(hdr)} cols (no stubs)")
