"""CHAMP/CHBMP independent sensitivity benchmark for the DISCERN variant engine (v3.1 Track A2 ext).

CHAMP (CDC Hemophilia A Mutation Project, F8) and CHBMP (CDC Hemophilia B Mutation Project, F9) are
public-use catalogs of variants *reported in patients with haemophilia* - i.e. an independent,
curated set of disease-causing F8/F9 alleles with ISTH-criteria severity + inhibitor history
(Payne 2013 PMID 23280990; Li 2013 PMID 24498619; 2022 lists, downloaded 2024-11-07 from
cdc.gov/hemophilia/mutation-project). Verified in DISCERN_v1_Dataset_Map_Verification_Report.md.

Truth label for every catalogued allele = PATHOGENIC. This is therefore a SENSITIVITY (recall)
benchmark of DISCERN's gene-specific (CFD-VCEP F8/F9) scoring, broken out by variant consequence:

  * NULL (nonsense / frameshift / canonical-splice +/-1,2): scored on CONSEQUENCE ALONE
    (PVS1 decision tree, Abou Tayoun 2018) - no predictor/frequency input needed. Expect LP/P.
  * MISSENSE / in-frame / synonymous / UTR / promoter: intrinsic-only CEILING. With the CFD-VCEP
    F8/F9 spec a missense tops out at PP3_Supporting+PM2_Supporting = 2 pts = VUS without routed
    PS3/PP1/PP4 - this is the designed limitation the disease-variant coupling addresses.

No PHI, no patient-level data - a public variant catalog. Run: `python -m eval.champ_chbmp_benchmark`.
"""
from __future__ import annotations

import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import openpyxl

from rules.point_engine import Classification
from rules.variant_scoring import Annotations, score_variant

# Location of the downloaded CDC files (outside the repo - not committed; public CDC data).
# Override with the CHAMP_CHBMP_DIR environment variable; defaults to a repo-relative data dir.
DATA_DIR = Path(os.environ.get("CHAMP_CHBMP_DIR", "data/champ_chbmp"))
FILES = {
    "F8": ("CHAMP-Variant-List-2022.xlsx", "CHAMP Variant List", 2351, 26),   # NP_000123 precursor, 26 exons
    "F9": ("CHBMP-Variant-List-2022.xlsx", "CHBMP Variant List", 461, 8),      # NP_000124 precursor, 8 exons
}

# CDC "Variant Type" -> (DISCERN consequence, is_null_LOF). Non-LOF map to a benign-ish consequence.
NULL_TYPES = {
    "nonsense": "nonsense",
    "frameshift": "frameshift",
    "splice site change": "canonical_splice",   # only +/-1,2 kept as PVS1-eligible (see below)
}
STRUCT_TYPES = {"large structural change", "small structural change"}  # CNV - reported separately


def _norm(s) -> str:
    return re.sub(r"\s+", " ", str(s or "").strip().lower())


def _canonical_splice(hgvs_cdna: str) -> bool:
    """True iff an intronic +/-1 or +/-2 position (canonical donor/acceptor dinucleotide)."""
    # e.g. c.6115+1G>A, c.5999-2A>G, c.1538-1G>C ; exclude +3/-5/deep-intronic/exonic
    for m in re.finditer(r"[+-](\d+)", hgvs_cdna or ""):
        if m.group(1) in ("1", "2"):
            return True
    return False


def _exon_int(exon) -> int | None:
    m = re.match(r"^\s*(\d+)", str(exon or ""))
    return int(m.group(1)) if m else None


def load(gene: str):
    fname, sheet, prot_len, last_exon = FILES[gene]
    wb = openpyxl.load_workbook(DATA_DIR / fname, read_only=True, data_only=True)
    ws = wb[sheet]
    it = ws.iter_rows(values_only=True)
    next(it)  # header
    recs = []
    for r in it:
        if r[0] in (None, ""):
            continue
        recs.append({
            "gene": gene, "cdna": str(r[0]).strip(), "protein": str(r[2 if gene == "F8" else 3] or ""),
            "vtype": _norm(r[4 if gene == "F8" else 5]), "exon": _exon_int(r[6 if gene == "F8" else 7]),
            "codon": r[7 if gene == "F8" else 8], "prot_len": prot_len, "last_exon": last_exon,
        })
    wb.close()
    return recs


def to_annotations(rec) -> Annotations | None:
    """Map a catalog record to DISCERN Annotations using CONSEQUENCE + NMD/codon geometry only
    (no gnomAD/REVEL): a deliberately conservative, predictor-free lower bound for the null subset."""
    vt = rec["vtype"]
    null_key = next((k for k in NULL_TYPES if k in vt), None)
    if null_key is None:
        return None
    cons = NULL_TYPES[null_key]
    if cons == "canonical_splice" and not _canonical_splice(rec["cdna"]):
        return None  # non-canonical splice -> needs SpliceAI, not PVS1
    exon, last = rec["exon"], rec["last_exon"]
    if cons == "canonical_splice":
        # CDC lists the intron for splice sites (exon parse -> None); a canonical +/-1,2 site
        # disrupts splicing -> aberrant transcript predicted to undergo NMD (Tayoun Very Strong),
        # except when it affects the terminal exon/intron (then NMD-escape).
        nmd = (exon is None or exon < last)
    else:
        # PTC: NMD unless in the last (3'-most) exon, which escapes NMD.
        nmd = (exon is not None and exon < last)
    m = re.match(r"^\s*(\d+)", str(rec["codon"] or ""))
    codon = int(m.group(1)) if m else None
    removes = (codon is not None and codon < 0.9 * rec["prot_len"])  # truncation removes >10% C-term
    return Annotations(consequence=cons, nmd_predicted=nmd, removes_gt10pct=removes)


def run():
    PATHOGENIC = {Classification.P, Classification.LP}
    overall = {"null_total": 0, "null_path": 0}
    by_type = defaultdict(lambda: Counter())
    miscalls = []
    type_dist = {}

    for gene in ("F8", "F9"):
        recs = load(gene)
        type_dist[gene] = Counter(r["vtype"] for r in recs)
        for rec in recs:
            ann = to_annotations(rec)
            if ann is None:
                continue  # non-null (missense/in-frame/synonymous/UTR/promoter/non-canonical splice)
            sv = score_variant(gene, rec["cdna"], ann)
            is_path = sv.classification in PATHOGENIC
            overall["null_total"] += 1
            overall["null_path"] += int(is_path)
            key = (gene, ann.consequence)
            by_type[key]["total"] += 1
            by_type[key]["path"] += int(is_path)
            by_type[key][sv.classification.name] += 1
            if not is_path:
                miscalls.append((gene, rec["cdna"], rec["vtype"], ann.consequence, sv.classification.name,
                                 sv.points, ",".join(sv.codes)))

    print("=" * 78)
    print("CHAMP/CHBMP independent sensitivity benchmark - DISCERN F8/F9 (CFD-VCEP)")
    print("Truth = pathogenic (catalogued disease alleles). NULL subset scored on consequence alone.")
    print("=" * 78)
    for gene in ("F8", "F9"):
        tot = sum(type_dist[gene].values())
        print(f"\n[{gene}] catalog n={tot}  variant-type distribution:")
        for t, n in type_dist[gene].most_common():
            print(f"    {n:5d}  {t}")

    print("\n--- NULL-subset recall (LP/P) by consequence ---")
    for (gene, cons), c in sorted(by_type.items()):
        rec = 100.0 * c["path"] / c["total"] if c["total"] else 0.0
        bands = " ".join(f"{k}={c[k]}" for k in ("P", "LP", "VUS", "LB", "B") if c.get(k))
        print(f"  {gene:3} {cons:16} n={c['total']:4d}  LP/P={c['path']:4d} ({rec:5.1f}%)   [{bands}]")

    t, p = overall["null_total"], overall["null_path"]
    print(f"\n  OVERALL NULL subset: {p}/{t} = {100.0*p/t:.1f}% classified Likely-Pathogenic/Pathogenic")
    print("  (PVS1-driven, no gnomAD/REVEL input; adding PM2 absent-from-gnomAD only raises this)")

    if miscalls:
        print(f"\n  Non-pathogenic NULL calls ({len(miscalls)}) - inspect (last-exon NMD-escape etc.):")
        for m in miscalls[:25]:
            print("    ", m)

    # Missense ceiling (analytic): show the intrinsic-only band for a representative missense.
    print("\n--- MISSENSE intrinsic-only ceiling (why coupling is needed) ---")
    for gene in ("F8", "F9"):
        n_mis = type_dist[gene].get("missense", 0)
        # Best intrinsic case: strong predictor + absent gnomAD (PP3_Supporting + PM2_Supporting).
        sv = score_variant(gene, "missense", Annotations(consequence="missense", af=0.0, revel=0.95))
        print(f"  {gene}: {n_mis} catalogued missense; best intrinsic codes={sv.codes} "
              f"pts={sv.points} -> {sv.classification.name} (needs routed PS3/PP1/PP4 to reach LP)")
    return overall


if __name__ == "__main__":
    sys.exit(0 if run() else 1)
