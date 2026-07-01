"""Genome-wide per-code partition validation on the FULL ClinGen ERepo (v3.1 Track A1).

Generalizes eval/erepo_reconstruction.py by dropping the bleeding-cluster gene restriction:
the partition (rules/vcep/partition.owner) is gene-agnostic, so it runs over EVERY eRepo
variant. Reports the two A1 hypotheses:

  H1 (Gate G9) - partition vocabulary coverage genome-wide: recognized_codes / total_codes,
     with the ENUMERATED uncovered code vocabulary (the <=1% the partition does not map).
  H2 - inflation-prevented rate genome-wide with a Wilson 95% CI: among variants carrying a
     non-genetic owned code, how many a naive all-codes band would push above the
     intrinsic-only band. Compared against the 20.7% bleeding-panel result.

Run on the VM:  python3 eval/erepo_genomewide.py data/raw/erepo/erepo_classifications.tab
"""
from __future__ import annotations

import csv
import math
import sys
from collections import Counter

from rules.acmg_codes import code_points
from rules.point_engine import BANDS, Classification
from rules.vcep.partition import base_code, owner

NON_GENETIC = {"disease_pp4", "functional", "segregation", "phasing", "denovo"}
_ORDER = [Classification.B, Classification.LB, Classification.VUS, Classification.LP, Classification.P]


def _band(pts: float, ba1: bool) -> Classification:
    if ba1:
        return Classification.B
    for thr, cls in BANDS:
        if pts >= thr:
            return cls
    return Classification.B


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, centre - half), min(1.0, centre + half))


def run(path: str) -> dict:
    total = at_risk = inflated = 0
    total_codes = recognized = 0
    unknown: Counter = Counter()
    genes: Counter = Counter()
    with open(path, encoding="utf-8", errors="replace") as fh:
        r = csv.reader(fh, delimiter="\t")
        h = next(r)
        gi = h.index("HGNC Gene Symbol")
        ci = h.index("Applied Evidence Codes (Met)")
        for row in r:
            if len(row) <= max(gi, ci):
                continue
            codes = [c.strip() for c in row[ci].split(",") if c.strip()]
            if not codes:
                continue
            all_pts = intrinsic = 0.0
            ba1 = has_owned = False
            for c in codes:
                base = base_code(c)
                fac = owner(c)
                total_codes += 1
                if fac is None:
                    unknown[base] += 1
                else:
                    recognized += 1
                pts, is_ba1 = code_points(c)
                if is_ba1:
                    ba1 = True
                all_pts += pts
                if fac == "variant_intrinsic":
                    intrinsic += pts
                elif fac in NON_GENETIC:
                    has_owned = True
            total += 1
            genes[row[gi].strip()] += 1
            if has_owned:
                at_risk += 1
                if _ORDER.index(_band(all_pts, ba1)) > _ORDER.index(_band(intrinsic, ba1)):
                    inflated += 1
    cov = recognized / total_codes if total_codes else None
    inf_lo, inf_hi = wilson(inflated, total)
    return {
        "variants": total,
        "genes": len(genes),
        "total_codes": total_codes,
        "partition_coverage": round(cov, 5) if cov is not None else None,
        "uncovered_codes": unknown.most_common(),
        "n_uncovered": sum(unknown.values()),
        "variants_with_owned_codes": at_risk,
        "at_risk_fraction": round(at_risk / total, 4) if total else None,
        "inflation_prevented": inflated,
        "inflation_prevented_fraction": round(inflated / total, 4) if total else None,
        "inflation_prevented_ci95": [round(inf_lo, 4), round(inf_hi, 4)],
    }


def main():
    s = run(sys.argv[1])
    print(f"Genome-wide ERepo: {s['variants']} variants across {s['genes']} genes, {s['total_codes']} codes")
    print(f"H1 partition coverage: {s['partition_coverage']:.4%}  "
          f"({s['n_uncovered']} uncovered codes: {s['uncovered_codes'] or 'none'})")
    print(f"H2 at-risk (non-genetic owned codes): {s['variants_with_owned_codes']} ({s['at_risk_fraction']:.1%})")
    ci = s["inflation_prevented_ci95"]
    print(f"H2 inflation prevented: {s['inflation_prevented']} "
          f"({s['inflation_prevented_fraction']:.1%}, 95% CI {ci[0]:.1%}-{ci[1]:.1%}) "
          f"[bleeding-panel reference: 20.7%]")


if __name__ == "__main__":
    main()
