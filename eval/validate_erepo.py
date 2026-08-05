"""Gate G1 - does the deterministic rule engine reproduce ClinGen eRepo calls?

Feeds each eRepo record's *expert-applied* evidence codes through the same
``point_engine`` the system uses, and compares the predicted ACMG class to the
expert's assertion. This validates the combining rule (points -> class) before any
adapter is trusted.

Run (on the VM, where the data lives):
    python -m eval.validate_erepo data/raw/erepo/erepo_classifications.tab
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter

from core.schemas import Classification, PointsLedger, Variant
from rules.acmg_codes import contributions_from_codes
from rules.point_engine import classify

# eRepo assertion text -> our Classification
_ASSERT = {
    "pathogenic": Classification.P,
    "likely pathogenic": Classification.LP,
    "uncertain significance": Classification.VUS,
    "likely benign": Classification.LB,
    "benign": Classification.B,
}
_ORDER = [Classification.B, Classification.LB, Classification.VUS,
          Classification.LP, Classification.P]


def _map_assertion(text: str):
    return _ASSERT.get(text.strip().lower())


def run(path: str, limit: int | None = None) -> dict:
    import pandas as pd

    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    code_col = "Applied Evidence Codes (Met)"
    assert_col = "Assertion"
    if limit:
        df = df.head(limit)

    total = len(df)
    evaluable = exact = adjacent = 0
    unmapped = 0
    confusion: Counter = Counter()       # (expert, predicted) -> n
    discordant_examples = []

    for _, row in df.iterrows():
        expert = _map_assertion(row.get(assert_col, ""))
        if expert is None:
            unmapped += 1
            continue
        codes = [c for c in row.get(code_col, "").split(",") if c.strip()]
        contribs, ba1 = contributions_from_codes(codes)
        led = PointsLedger(Variant("0", 0, "N", "N", gene=row.get("HGNC Gene Symbol", "")),
                           contributions=contribs)
        pred = Classification.B if ba1 else classify(led)
        evaluable += 1
        confusion[(expert.name, pred.name)] += 1
        di, dp = _ORDER.index(expert), _ORDER.index(pred)
        if di == dp:
            exact += 1
        if abs(di - dp) <= 1:
            adjacent += 1
        elif len(discordant_examples) < 15:
            discordant_examples.append({
                "variant": row.get("#Variation", "")[:60],
                "gene": row.get("HGNC Gene Symbol", ""),
                "expert": expert.name, "predicted": pred.name,
                "points": led.points, "codes": ",".join(codes),
            })

    summary = {
        "total_rows": total,
        "evaluable": evaluable,
        "unmapped_assertions": unmapped,
        "exact_concordance": round(exact / evaluable, 4) if evaluable else None,
        "within_one_bin": round(adjacent / evaluable, 4) if evaluable else None,
        "confusion_matrix": {f"{k[0]}->{k[1]}": v for k, v in sorted(confusion.items())},
        "discordant_examples": discordant_examples,
    }
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()
    s = run(args.path, args.limit)
    if args.json:
        print(json.dumps(s, indent=2))
        return
    print(f"eRepo records: {s['total_rows']}  evaluable: {s['evaluable']}  "
          f"(unmapped assertions: {s['unmapped_assertions']})")
    print(f"EXACT concordance:   {s['exact_concordance']:.1%}")
    print(f"Within-one-bin:      {s['within_one_bin']:.1%}")
    print("\nConfusion (expert -> predicted):")
    for k, v in s["confusion_matrix"].items():
        print(f"  {k:20s} {v}")
    if s["discordant_examples"]:
        print("\nSample >1-bin discordances:")
        for d in s["discordant_examples"][:10]:
            print(f"  {d['gene']:8s} expert={d['expert']:4s} pred={d['predicted']:4s} "
                  f"pts={d['points']:+.0f}  [{d['codes']}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
