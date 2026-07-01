"""gnomAD per-variant frequency cross-check (v3.1 Track A1, last sub-part; Gate G10).

Tests whether DISCERN's GENE-SPECIFIC CSpec frequency thresholds (spec.af), applied to the
REAL per-variant gnomAD allele frequency, reproduce the frequency code (BA1/BS1/PM2) the VCEP
actually applied. The per-variant gnomAD AFs are the values the expert curators CITE VERBATIM
in the eRepo "Summary of interpretation" narratives (e.g. "reported at a MAF of 0.00635 in the
African subpopulation in gnomAD v2.1.1") - authoritative and already on the VM, so no multi-TB
gnomAD download is needed. Restricted to genes with a verified DISCERN spec.

Run:  python3 -m eval.gnomad_freq_check data/raw/erepo/erepo_classifications.tab
"""
from __future__ import annotations

import csv
import re
import sys
from collections import Counter

from rules.vcep.loader import get_spec

# genes with a verified gene-specific spec (frequency thresholds extracted + checked)
SPEC_GENES = {"ITGA2B", "ITGB3", "F8", "F9", "VWF", "GP1BA"}

# observed AF the curator cites (NOT the threshold). Avoid "cut-off/criteria/threshold of".
_OBS = re.compile(
    r"(?:reported at|observed at|present at|filtering allele frequency of|"
    r"grpmax filtering allele frequency of|MAF of|frequency of|FAF of)\s*"
    r"(?:a\s+)?(?:high\s+)?(?:of\s+)?(\d+\.\d+(?:[eE]-?\d+)?)", re.I)
_THRESH_CTX = re.compile(r"(cut[- ]?off|criteria|threshold)", re.I)
# "absent / not observed in gnomAD" -> observed AF effectively 0 (the PM2-rare case)
_ABSENT = re.compile(r"(absent|not (?:been )?(?:observed|present|detected|reported|found))"
                     r"[^.]{0,40}(gnomad|population|databas)", re.I)


def _applied_freq_code(codes: str):
    bases = {t.strip().split("_")[0] for t in codes.split(",")}   # PM2_Supporting -> PM2
    for c in ("BA1", "BS1", "PM2"):
        if c in bases:
            return c
    return None


def _observed_af(narr: str):
    """Largest curator-cited observed AF in the narrative (subpop/grpmax = what BA1/BS1 use)."""
    afs = []
    for m in _OBS.finditer(narr):
        # skip if this number is part of a 'threshold/cut-off of X' phrase
        ctx = narr[max(0, m.start() - 30):m.start()]
        if _THRESH_CTX.search(ctx):
            continue
        try:
            afs.append(float(m.group(1)))
        except ValueError:
            pass
    if afs:
        return max(afs)
    if _ABSENT.search(narr):           # "absent in gnomAD" -> AF ~ 0 (the rare/PM2 case)
        return 0.0
    return None


def _predicted_code(spec, af: float):
    if af >= spec.af_threshold("ba1"):
        return "BA1"
    if af >= spec.af_threshold("bs1"):
        return "BS1"
    if af < spec.af_threshold("pm2"):
        return "PM2"
    return None                                   # between bs1 and pm2 -> no frequency code


def run(path: str) -> dict:
    total = with_af = with_code = checked = concordant = 0
    confusion: Counter = Counter()
    with open(path, encoding="utf-8", errors="replace") as fh:
        r = csv.reader(fh, delimiter="\t")
        h = next(r)
        gi, ci, si = (h.index("HGNC Gene Symbol"), h.index("Applied Evidence Codes (Met)"),
                      h.index("Summary of interpretation"))
        for row in r:
            if len(row) <= max(gi, ci, si) or row[gi].strip() not in SPEC_GENES:
                continue
            total += 1
            spec = get_spec(row[gi].strip())
            applied = _applied_freq_code(row[ci])
            af = _observed_af(row[si])
            if af is not None:
                with_af += 1
            if applied:
                with_code += 1
            if af is None or applied is None:
                continue
            checked += 1
            pred = _predicted_code(spec, af)
            ok = pred == applied
            concordant += ok
            confusion[f"VCEP:{applied} -> DISCERN:{pred}"] += 1
    return {
        "spec_gene_variants": total,
        "variants_with_parsed_af": with_af,
        "variants_with_freq_code": with_code,
        "cross_checked": checked,
        "concordance": round(concordant / checked, 4) if checked else None,
        "confusion": confusion.most_common(),
    }


def main():
    s = run(sys.argv[1])
    print(f"Spec-gene variants: {s['spec_gene_variants']}  "
          f"(parsed gnomAD AF in {s['variants_with_parsed_af']}, freq code in {s['variants_with_freq_code']})")
    print(f"Cross-checked (AF + applied freq code both present): {s['cross_checked']}")
    print(f"DISCERN gene-specific threshold reproduces the VCEP freq code: "
          f"{s['concordance']:.1%}" if s['concordance'] is not None else "n/a")
    print("Confusion (VCEP applied -> DISCERN predicted):")
    for k, v in s["confusion"]:
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
