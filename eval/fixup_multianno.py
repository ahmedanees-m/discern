"""Reconcile an ANNOVAR multianno (built from available DBs) to the column contract InterVar
expects: rename gnomAD-exome AF* -> gnomAD_genome_*, CLNSIG -> CLINSIG, phyloP100way ->
phyloP46way_placental; stub the DBs we could not fetch (avsnp147, 1000g2015aug_all,
AAChange.knownGene). Position-independent (InterVar finds columns by name). v3.1 InterVar run.
"""
import sys

RENAME = {
    "AF": "gnomAD_genome_ALL", "AF_afr": "gnomAD_genome_AFR", "AF_amr": "gnomAD_genome_AMR",
    "AF_eas": "gnomAD_genome_EAS", "AF_nfe": "gnomAD_genome_NFE", "AF_fin": "gnomAD_genome_FIN",
    "AF_asj": "gnomAD_genome_ASJ", "AF_oth": "gnomAD_genome_OTH", "CLNSIG": "CLINSIG",
    "phyloP100way_vertebrate": "phyloP46way_placental",
}
STUB = ["avsnp147", "1000g2015aug_all", "AAChange.knownGene"]

inp, out = sys.argv[1], sys.argv[2]
with open(inp, encoding="utf-8", errors="replace") as f, open(out, "w", encoding="utf-8") as o:
    hdr = f.readline().rstrip("\n").split("\t")
    hdr = [RENAME.get(c, c) for c in hdr]
    o.write("\t".join(hdr + STUB) + "\n")
    n = 0
    for line in f:
        cols = line.rstrip("\n").split("\t")
        o.write("\t".join(cols + ["."] * len(STUB)) + "\n")
        n += 1
print(f"fixed multianno: {n} rows, {len(hdr) + len(STUB)} cols")
