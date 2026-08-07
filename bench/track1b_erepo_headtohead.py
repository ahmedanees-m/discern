"""TRACK 1' - eRepo-primary + time-split variant head-to-head (Phase A2 Part A).

Re-runs the Track-1 comparison on the FDA-recognized ClinGen eRepo expert-panel surface instead of
raw ClinVar - removing the circularity caveat. Two surfaces:
  - eRepo-primary : all eRepo P/B variants in the cluster genes (expert labels, not raw ClinVar).
  - time-split    : eRepo variants expert-approved AFTER a comparator's bundled ClinVar snapshot
                    (default 2021-05-01, InterVar's clinvar_20210501) - no tool could have memorised
                    the label. The gold circularity guard.

Tools scored here (disk-free, from the GeneBe eRepo annotation cache + local DISCERN):
  DISCERN (calibrated), REVEL, AlphaMissense, and GeneBe as the circularity EXHIBIT only (live
  ClinVar -> cannot be time-split). InterVar (legacy anchor) is added by the VM step
  (track1b_intervar.py) once its hg38 DBs are staged; its prior ClinVar-set number is cited meanwhile.

Adds the per-ACMG-code Cohen's kappa of DISCERN's applied codes vs eRepo's Applied Evidence Codes -
the partition made visible against an expert surface.

Run:  python -m bench.track1b_erepo_headtohead   (needs bench/data/genebe_erepo.jsonl)
"""
from __future__ import annotations

import json
import os

from sklearn.metrics import brier_score_loss, cohen_kappa_score

from bench.track1_variant_headtohead import (
    _cons,
    _cont_metrics,
    _disc_prob_oof,
    _ece,
)
from rules.variant_scoring import Annotations, score_variant

# The applied-code vocabulary is owned by rules.vcep.partition, so this analysis and the
# partition-coverage analysis cannot drift apart.
from rules.vcep.partition import applied_codes as _erepo_codes
from rules.vcep.partition import applied_codes_with_strength as _erepo_codes_with_strength
from rules.vcep.partition import default_strength

HERE = os.path.dirname(__file__)
CACHE = os.path.join(HERE, "data", "genebe_erepo.jsonl")
OUT_JSON = os.path.join(HERE, "track1b_erepo_metrics.json")
OUT_CSV = os.path.join(HERE, "track1b_erepo_headtohead.csv")
TIMESPLIT_AFTER = "2021-05-01"   # InterVar clinvar_20210501 bundle date

PATHL = {"P", "LP"}
BENL = {"B", "LB"}
def _label(a):
    return 1 if a in PATHL else 0 if a in BENL else None


def load_rows():
    rows = []
    for line in open(CACHE, encoding="utf-8"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        r = json.loads(line)
        cons = _cons(r["effect"])
        af = r.get("gnomad_exomes_af") or r.get("gnomad_genomes_af")
        rev = float(r["revel_score"]) if r.get("revel_score") is not None else None
        sv = score_variant(r["q_gene"], r["q_gene"], Annotations(
            af=float(af) if af is not None else None, revel=rev, consequence=cons))
        rows.append({
            "gene": r["q_gene"], "assertion": r["q_assertion"], "label": _label(r["q_assertion"]),
            "approval_date": r.get("q_approval_date", ""), "effect": r["effect"],
            "is_missense": r["effect"] == "missense_variant",
            "discern_points": sv.points, "discern_class": sv.classification.name,
            "discern_codes": {c.split("_")[0] for c in sv.codes},   # base ACMG code (strip strength suffix)
            "discern_strength": {c.split("_")[0]: (c.split("_", 1)[1].replace("_", " ") if "_" in c
                                                   else default_strength(c.split("_")[0]))
                                 for c in sv.codes},
            "erepo_codes": _erepo_codes(r.get("q_codes_met", "")),
            "erepo_strength": _erepo_codes_with_strength(r.get("q_codes_met", "")),
            "genebe_class": r.get("acmg_classification"), "genebe_score": r.get("acmg_score"),
            "revel": r.get("revel_score"), "alphamissense": r.get("alphamissense_score"),
        })
    return rows


def _surface_metrics(rows):
    mis = [r for r in rows if r["is_missense"] and r["label"] is not None]
    y = [r["label"] for r in mis]
    n_path = sum(y)
    out = {"n_missense_PB": len(mis), "n_path": n_path, "n_benign": len(mis) - n_path}
    out["discrimination"] = {
        "DISCERN": _cont_metrics(mis, "discern_points", y),
        "REVEL": _cont_metrics(mis, "revel", y),
        "AlphaMissense": _cont_metrics(mis, "alphamissense", y),
        "GeneBe_exhibit": _cont_metrics(mis, "genebe_score", y),
    }
    # calibration (continuous tools)
    disc_oof = _disc_prob_oof([r["discern_points"] for r in mis], y)
    cal = {"DISCERN_isotonic_oof": {"ece": round(_ece(disc_oof, y), 4),
                                    "brier": round(brier_score_loss(y, disc_oof), 4)}}
    for name, key in (("REVEL", "revel"), ("AlphaMissense", "alphamissense")):
        sc = [(r[key], yy) for r, yy in zip(mis, y, strict=False) if r[key] is not None]
        p = [float(a) for a, _ in sc]
        yy = [b for _, b in sc]
        cal[name + "_raw"] = {"ece": round(_ece(p, yy), 4), "brier": round(brier_score_loss(yy, p), 4)}
    out["calibration"] = cal
    return out


def _per_code_kappa(rows):
    """Agreement per criterion, at two levels.

    Criterion level asks whether the criterion was applied at all; strength level asks, among the
    variants where both applied it, whether they applied it at the same ClinGen strength. A single
    kappa conflates the two, and the second is the more informative comparison for criteria whose
    default strength the panels routinely modify.
    """
    pb = [r for r in rows if r["label"] is not None]
    codes = ["PVS1", "PS1", "PS3", "PS4", "PM1", "PM2", "PM5", "PP3", "BA1", "BS1", "BS2", "BP4", "BP7"]
    out = {}
    for c in codes:
        d = [1 if c in r["discern_codes"] else 0 for r in pb]
        e = [1 if c in r["erepo_codes"] else 0 for r in pb]
        if sum(e) == 0 and sum(d) == 0:
            continue
        k = cohen_kappa_score(e, d) if len(set(e)) > 1 and len(set(d)) > 1 else None
        both = [r for r in pb if c in r["discern_codes"] and c in r["erepo_codes"]]
        same = sum(1 for r in both
                   if r["erepo_strength"].get(c) == r["discern_strength"].get(c))
        out[c] = {"erepo_applied": sum(e), "discern_applied": sum(d),
                  "kappa": round(k, 3) if k is not None else None,
                  "both_applied": len(both),
                  "same_strength": same,
                  "strength_agreement": round(same / len(both), 3) if both else None}
    return out


def run():
    rows = load_rows()
    primary = _surface_metrics(rows)
    ts_rows = [r for r in rows if (r["approval_date"] or "") > TIMESPLIT_AFTER]
    timesplit = _surface_metrics(ts_rows)
    out = {
        "set": "eRepo cluster-gene expert-panel surface (GRCh38)",
        "n_total": len(rows), "n_timesplit": len(ts_rows), "timesplit_after": TIMESPLIT_AFTER,
        "eRepo_primary": primary, "time_split": timesplit,
        "per_code_kappa_vs_erepo": _per_code_kappa(rows),
        "note_genebe": "GeneBe shown as circularity exhibit only (queries live ClinVar; cannot be time-split).",
        "note_intervar": "InterVar (legacy anchor) added by the VM step once hg38 DBs staged; prior ClinVar-set missense AUROC 0.811.",
    }
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    return out


def main():
    o = run()
    print("== TRACK 1' : eRepo-primary + time-split (missense, expert-panel truth) ==")
    for surf in ("eRepo_primary", "time_split"):
        s = o[surf]
        n = "" if surf == "eRepo_primary" else f" (approved > {o['timesplit_after']})"
        print(f"\n-- {surf}{n}: missense P/B={s['n_missense_PB']} (path={s['n_path']} benign={s['n_benign']}) --")
        for t, m in s["discrimination"].items():
            if m:
                tag = " [exhibit: live ClinVar]" if t == "GeneBe_exhibit" else ""
                print(f"   {t:16} AUROC={m['auroc']:.3f}  AUPRC={m['auprc']:.3f}  sens@90={m['sens@90spec']:.3f}{tag}")
        print("   calibration ECE:", {k: v["ece"] for k, v in s["calibration"].items()})
    print("\n-- per-ACMG-code kappa (DISCERN applied vs eRepo applied) --")
    for c, m in o["per_code_kappa_vs_erepo"].items():
        print(f"   {c:5} eRepo={m['erepo_applied']:4} DISCERN={m['discern_applied']:4} kappa={m['kappa']}")
    print(f"\nwrote {os.path.relpath(OUT_JSON)}")


if __name__ == "__main__":
    main()
