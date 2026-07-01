"""ACMG combining-rule fidelity + per-code partition on real ClinGen ERepo (Tier A).

Runs on the REAL ClinGen Evidence Repository export, restricted to the DISCERN
bleeding-cluster genes. It reports FOUR falsifiable metrics (the earlier "no-inflation
rate" was a tautology - it summed the same points twice and returned 100% for any input -
and has been removed):

  1. acmg_combining_exact - given the experts' OWN applied codes, does the Tavtigian point
     engine reproduce the VCEP's bottom-line label? This validates the point values and
     banding ARITHMETIC only; it does NOT validate code assignment (the codes are given)
     nor the disease-variant coupling (no phenotype is present here).
  2. partition_coverage - every applied code is recognised and routed to exactly one
     owning factor (target: 1.0, zero unknown). This CAN fail - it flags any ACMG code the
     partition vocabulary does not cover - which is what makes it a real check.
  3. at-risk population + points routed out of the genetic stream - how many real variants
     carry codes owned by non-genetic factors (PP4/PS3/PP1/PM3/...) and how many points
     those move out of P(E_geno|V). This is exactly the evidence a tool consuming the
     VCEP's bottom-line label would double-count.
  4. inflation_prevented - among the at-risk variants, how many a naive all-codes genetic
     score would push into a HIGHER pathogenicity band than the intrinsic-only score (i.e.
     how often the moved evidence was band-determining). Expect this to be modest; that is
     the result. (The coupling's empirical payoff is unprovable on ERepo - there is
     no paired phenotype to double-count against - and is a cohort result, not this one.)

Run (on the VM, where the data lives):
    python3 -m eval.erepo_reconstruction data/raw/erepo/erepo_classifications.tab
"""
from __future__ import annotations

import csv
import sys
from collections import Counter

from rules.acmg_codes import code_points
from rules.point_engine import BANDS, Classification
from rules.vcep.partition import base_code, owner

# DISCERN cluster genes (the bleeding/platelet panel).
CLUSTER_GENES = {
    "ITGA2B", "ITGB3", "FERMT3", "RASGRP2", "ITGB2",
    "VWF", "GP1BA", "GP1BB", "GP9",
    "MYH9", "TUBB1", "ACTN1", "FLNA",
    "NBEAL2", "GFI1B", "LYST", "HPS1", "HPS3", "HPS4", "HPS5", "HPS6", "AP3B1",
    "RUNX1", "ETV6", "ANKRD26",
    "F8", "F9", "F11", "F7", "F10", "F5", "F13A1", "F13B", "FGA", "FGB", "FGG", "F2",
}

# factors whose codes belong to streams OTHER than the variant-intrinsic genetic factor.
NON_GENETIC_FACTORS = {"disease_pp4", "functional", "segregation", "phasing", "denovo"}

_ASSERT = {
    "pathogenic": Classification.P, "likely pathogenic": Classification.LP,
    "uncertain significance": Classification.VUS, "likely benign": Classification.LB,
    "benign": Classification.B,
}
_ORDER = [Classification.B, Classification.LB, Classification.VUS, Classification.LP, Classification.P]


def _band(pts: float, ba1: bool) -> Classification:
    if ba1:
        return Classification.B
    for thr, cls in BANDS:
        if pts >= thr:
            return cls
    return Classification.B


def run(path: str) -> dict:
    exact = adjacent = total = 0
    at_risk = inflated = 0
    total_codes = recognized_codes = 0
    unknown_codes: Counter = Counter()
    gene_counts: Counter = Counter()
    code_strengths: Counter = Counter()
    routed_out: Counter = Counter()          # non-genetic factor -> points moved out

    with open(path, encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        gi = header.index("HGNC Gene Symbol")
        ci = header.index("Applied Evidence Codes (Met)")
        ai = header.index("Assertion")
        for row in reader:
            if len(row) <= max(gi, ci, ai):
                continue
            gene = row[gi].strip()
            if gene not in CLUSTER_GENES:
                continue
            expert = _ASSERT.get(row[ai].strip().lower())
            if expert is None:
                continue
            codes = [c.strip() for c in row[ci].split(",") if c.strip()]

            all_pts = intrinsic_pts = 0.0
            ba1 = has_owned = False
            for c in codes:
                base = base_code(c)
                fac = owner(c)
                total_codes += 1
                if fac is None:
                    unknown_codes[base] += 1
                else:
                    recognized_codes += 1
                pts, is_ba1 = code_points(c)
                if is_ba1:
                    ba1 = True
                code_strengths[base] += 1
                all_pts += pts                                  # naive: a tool using ALL codes
                if fac == "variant_intrinsic":
                    intrinsic_pts += pts                        # the genetic factor uses only these
                elif fac in NON_GENETIC_FACTORS:
                    has_owned = True
                    routed_out[fac] += pts                      # routed out of the genetic stream

            total += 1
            gene_counts[gene] += 1

            # (1) ACMG combining fidelity: the experts' full code set -> band vs the assertion.
            naive_band = _band(all_pts, ba1)
            di, dp = _ORDER.index(expert), _ORDER.index(naive_band)
            if di == dp:
                exact += 1
            if abs(di - dp) <= 1:
                adjacent += 1

            # (3)/(4) at-risk + inflation prevented.
            if has_owned:
                at_risk += 1
                intrinsic_band = _band(intrinsic_pts, ba1)
                if _ORDER.index(naive_band) > _ORDER.index(intrinsic_band):
                    inflated += 1

    return {
        "bleeding_variants": total,
        "genes_covered": len(gene_counts),
        # (1) arithmetic fidelity, not reasoning, not coupling
        "acmg_combining_exact": round(exact / total, 4) if total else None,
        "acmg_combining_within_one_bin": round(adjacent / total, 4) if total else None,
        # (2) partition coverage - can fail
        "partition_coverage": round(recognized_codes / total_codes, 4) if total_codes else None,
        "unknown_codes": dict(unknown_codes),
        # (3) at-risk population + points routed out of the genetic stream
        "variants_with_owned_codes": at_risk,
        "at_risk_fraction": round(at_risk / total, 4) if total else None,
        "points_routed_out": dict(routed_out),
        # (4) inflation prevented (band-determining among at-risk)
        "inflation_prevented": inflated,
        "inflation_prevented_fraction": round(inflated / total, 4) if total else None,
        "top_genes": gene_counts.most_common(12),
        "per_code_strength_distribution": code_strengths.most_common(20),
    }


def main():
    s = run(sys.argv[1])
    print(f"Bleeding-gene VCEP variants: {s['bleeding_variants']} across {s['genes_covered']} genes")
    print(f"(1) ACMG combining EXACT:      {s['acmg_combining_exact']:.1%}  "
          f"(within-1-bin {s['acmg_combining_within_one_bin']:.1%}) - arithmetic only")
    print(f"(2) Partition coverage:        {s['partition_coverage']:.1%}  "
          f"(unknown codes: {s['unknown_codes'] or 'none'})")
    print(f"(3) Variants with owned codes: {s['variants_with_owned_codes']} "
          f"({s['at_risk_fraction']:.1%} at risk of double-count)")
    print(f"    Points routed out of the genetic stream: {s['points_routed_out']}")
    print(f"(4) Inflation prevented:       {s['inflation_prevented']} "
          f"({s['inflation_prevented_fraction']:.1%}) - naive all-codes band > intrinsic-only band")
    print("Top genes:", s["top_genes"][:8])
    print("Per-code strength distribution (top):", s["per_code_strength_distribution"][:12])


if __name__ == "__main__":
    main()
