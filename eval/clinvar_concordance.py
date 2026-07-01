"""ClinVar concordance of DISCERN's intrinsic-only band (v3.1 Track A1; Gate G10 / H3).

Tests whether DISCERN's VARIANT-INTRINSIC-ONLY band (after routing PP4/PS3/PP1/PM3 etc. out
to the other factors) still agrees with the ClinVar community label. Joins each ERepo variant
to ClinVar by VariationID, stratifies by review-status stars, and reports exact + within-one-bin
agreement per star level. Note: 3-star (expert-panel) ClinVar entries are often the
VCEP's own submission, so the 2-star (multiple-submitter) stratum is the more INDEPENDENT test.

Run on the VM:
    python3 eval/clinvar_concordance.py data/raw/erepo/erepo_classifications.tab \
        data/raw/clinvar/variant_summary.txt.gz
"""
from __future__ import annotations

import csv
import gzip
import sys
from collections import Counter, defaultdict

from rules.acmg_codes import code_points
from rules.point_engine import BANDS, Classification
from rules.vcep.partition import owner

_ORDER = [Classification.B, Classification.LB, Classification.VUS, Classification.LP, Classification.P]


def _band(pts: float, ba1: bool) -> Classification:
    if ba1:
        return Classification.B
    for thr, cls in BANDS:
        if pts >= thr:
            return cls
    return Classification.B


def _intrinsic_band(codes: list[str]) -> Classification:
    pts = 0.0
    ba1 = False
    for c in codes:
        p, is_ba1 = code_points(c)
        if is_ba1:
            ba1 = True
        if owner(c) == "variant_intrinsic":
            pts += p
    return _band(pts, ba1)


def _clinvar_class(txt: str):
    t = txt.lower()
    if "conflict" in t:
        return None
    if "pathogenic/likely pathogenic" in t or "pathogenic, likely pathogenic" in t:
        return Classification.P
    if "benign/likely benign" in t or "benign, likely benign" in t:
        return Classification.B
    if "likely pathogenic" in t:
        return Classification.LP
    if "likely benign" in t:
        return Classification.LB
    if "pathogenic" in t:
        return Classification.P
    if "benign" in t:
        return Classification.B
    if "uncertain" in t:
        return Classification.VUS
    return None


def _stars(review: str) -> int:
    r = review.lower()
    if "practice guideline" in r:
        return 4
    if "reviewed by expert panel" in r:
        return 3
    if "multiple submitters" in r and "no conflict" in r:
        return 2
    if "conflicting" in r:
        return 1
    if "single submitter" in r:
        return 1
    return 0


def _col(header: list[str], *candidates: str) -> int:
    for c in candidates:
        if c in header:
            return header.index(c)
    raise KeyError(f"none of {candidates} in ClinVar header")


def run(erepo_path: str, clinvar_path: str) -> dict:
    # 1) ERepo -> {VariationID: intrinsic_band}
    want: dict[str, Classification] = {}
    with open(erepo_path, encoding="utf-8", errors="replace") as fh:
        r = csv.reader(fh, delimiter="\t")
        h = next(r)
        vi = h.index("ClinVar Variation Id")
        ci = h.index("Applied Evidence Codes (Met)")
        for row in r:
            if len(row) <= max(vi, ci):
                continue
            vid = row[vi].strip()
            codes = [c.strip() for c in row[ci].split(",") if c.strip()]
            if vid and codes:
                want[vid] = _intrinsic_band(codes)

    # 2) stream ClinVar variant_summary; capture class + stars for wanted VariationIDs
    seen: dict[str, tuple] = {}
    opener = gzip.open if clinvar_path.endswith(".gz") else open
    with opener(clinvar_path, "rt", encoding="utf-8", errors="replace") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        header = [c.lstrip("#") for c in header]
        vid_i = _col(header, "VariationID")
        cls_i = _col(header, "Germline classification", "ClinicalSignificance")
        rev_i = _col(header, "Germline review status", "ReviewStatus")
        for row in r:
            if len(row) <= max(vid_i, cls_i, rev_i):
                continue
            vid = row[vid_i].strip()
            if vid not in want or vid in seen:
                continue
            cv = _clinvar_class(row[cls_i])
            if cv is None:
                continue
            seen[vid] = (cv, _stars(row[rev_i]))

    # 3) agreement, stratified by stars
    buckets: dict[str, dict] = defaultdict(lambda: {"n": 0, "exact": 0, "within1": 0})
    cls_dist: Counter = Counter()
    for vid, dband in want.items():
        if vid not in seen:
            continue
        cv, stars = seen[vid]
        cls_dist[stars] += 1
        di, ci2 = _ORDER.index(dband), _ORDER.index(cv)
        for key in ("all", f"{stars}star", ">=2star" if stars >= 2 else "<2star"):
            b = buckets[key]
            b["n"] += 1
            b["exact"] += (di == ci2)
            b["within1"] += (abs(di - ci2) <= 1)

    out = {"erepo_with_varid": len(want), "matched_in_clinvar": len(seen),
           "stars_distribution": dict(sorted(cls_dist.items()))}
    for key, b in buckets.items():
        if b["n"]:
            out[key] = {"n": b["n"], "exact": round(b["exact"] / b["n"], 4),
                        "within1": round(b["within1"] / b["n"], 4)}
    return out


def main():
    s = run(sys.argv[1], sys.argv[2])
    print(f"ERepo variants with a ClinVar VariationID: {s['erepo_with_varid']}; "
          f"matched (non-conflicting) in ClinVar: {s['matched_in_clinvar']}")
    print(f"Stars distribution: {s['stars_distribution']}")
    for key in ("all", ">=2star", "<2star", "2star", "3star", "4star", "1star", "0star"):
        if key in s:
            b = s[key]
            print(f"  {key:9} n={b['n']:6}  exact={b['exact']:.1%}  within-1-bin={b['within1']:.1%}")


if __name__ == "__main__":
    main()
