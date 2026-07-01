"""TRACK 3 - trustworthiness layer (DISCERN_Benchmark_Execution_Plan_v1).

Foregrounds the three contributions no competitor ACMG/diagnosis paper reports together:
  3a calibration       - diagnosis-posterior reliability + ECE + confidently-wrong rate (variant-side
                         calibration is read from the Track-1 metrics).
  3b safety interlock  - hard-stop SENSITIVITY (fires when the planned tx is contraindicated) and
                         SPECIFICITY (stays silent when the planned tx is harmless).
  3c abstention        - risk-coverage curve: accuracy on the retained cases as a function of the
                         abstention rate (the operational payoff of calibration).

All local, no external tools. Run:  python -m bench.track3_trustworthiness
"""
from __future__ import annotations

import csv
import json
import os

import numpy as np

from diseases.ontology import cluster_for
from eval.curated_case_benchmark import _evidence, load_cases
from eval.reader_study import load_vignettes
from jointdx.abstain import decide
from jointdx.factorgraph import joint
from jointdx.infer import marginal_disease
from jointdx.orchestrate import diagnose

HERE = os.path.dirname(__file__)
OUT_JSON = os.path.join(HERE, "track3_metrics.json")
RC_CSV = os.path.join(HERE, "track3_risk_coverage.csv")
HARMLESS_TX = "tranexamic_acid"   # not a contraindication token in any cluster -> hard-stop must stay silent


def _ece(conf, correct, bins=10):
    conf, correct = np.asarray(conf, float), np.asarray(correct, float)
    edges = np.linspace(0, 1, bins + 1)
    e = 0.0
    for i in range(bins):
        hi = edges[i + 1] if i < bins - 1 else 1.0001
        m = (conf >= edges[i]) & (conf < hi)
        if m.sum():
            e += (m.sum() / len(conf)) * abs(correct[m].mean() - conf[m].mean())
    return float(e)


def _case_predictions():
    """Per curated case: leading confidence (decide.p) + whether the leading call is correct."""
    out = []
    for c in load_cases():
        cl = cluster_for(c["cluster"])
        ev = _evidence(c)
        md = marginal_disease(joint(cl, ev))
        lead = max(md, key=md.get)
        d = decide(cl, ev, n_mc=60)
        out.append({"id": c["id"], "true": c["true_dx"], "lead": lead,
                    "correct": int(lead == c["true_dx"]), "conf": float(d.p), "decided": bool(d.decided)})
    return out


def calibration_diagnosis(preds):
    conf = [p["conf"] for p in preds]
    correct = [p["correct"] for p in preds]
    ece = _ece(conf, correct)
    conf_wrong = sum(1 for p in preds if p["conf"] >= 0.80 and not p["correct"])
    return {"n": len(preds), "ece": round(ece, 4),
            "mean_confidence": round(float(np.mean(conf)), 4),
            "accuracy": round(float(np.mean(correct)), 4),
            "confidently_wrong_at_0.8": conf_wrong,
            "confidently_wrong_rate": round(conf_wrong / len(preds), 4)}


def safety_interlock():
    """Hard-stop sensitivity (real contraindicated tx) + specificity (harmless tx)."""
    vigs = load_vignettes()
    tier_s = [v for v in vigs if v.tier == "S" and v.harmful_tx]
    fired_real = silent_harmless = 0
    rows = []
    for v in tier_s:
        rec_real = diagnose(v.ev, planned_tx=v.harmful_tx, n_mc=80)
        rec_safe = diagnose(v.ev, planned_tx=HARMLESS_TX, n_mc=80)
        hs_real = any("HARD STOP" in f.message for f in (rec_real.safety_flags if rec_real else []))
        hs_safe = any("HARD STOP" in f.message for f in (rec_safe.safety_flags if rec_safe else []))
        fired_real += hs_real
        silent_harmless += (not hs_safe)
        rows.append({"id": v.id, "harmful_tx": v.harmful_tx, "hardstop_on_real": hs_real,
                     "hardstop_on_harmless": hs_safe})
    n = len(tier_s)
    return {"n_scenarios": n,
            "hardstop_sensitivity": round(fired_real / n, 4),
            "hardstop_specificity": round(silent_harmless / n, 4),
            "fired_when_contraindicated": f"{fired_real}/{n}",
            "silent_when_harmless": f"{silent_harmless}/{n}",
            "rows": rows}


def risk_coverage(preds):
    """Sweep a confidence threshold; report (coverage, accuracy_on_retained). Monotone => calibration pays off."""
    s = sorted(preds, key=lambda p: -p["conf"])
    curve = []
    for k in range(1, len(s) + 1):
        kept = s[:k]
        cov = k / len(s)
        acc = float(np.mean([p["correct"] for p in kept]))
        curve.append((round(cov, 4), round(acc, 4), round(s[k - 1]["conf"], 4)))
    # DISCERN's own decide() operating point
    dec = [p for p in preds if p["decided"]]
    op = {"coverage": round(len(dec) / len(preds), 4),
          "accuracy_on_decided": round(float(np.mean([p["correct"] for p in dec])), 4) if dec else None,
          "abstention_rate": round(1 - len(dec) / len(preds), 4)}
    # monotonicity: accuracy at low coverage (most confident) >= accuracy at full coverage
    acc_full = curve[-1][1]
    acc_top50 = next(a for c, a, _ in curve if c >= 0.5)
    return {"operating_point_decide": op, "accuracy_full_coverage": acc_full,
            "accuracy_top50pct_confident": acc_top50,
            "monotone_payoff": acc_top50 >= acc_full, "curve": curve}


def run():
    preds = _case_predictions()
    out = {"calibration_diagnosis": calibration_diagnosis(preds),
           "safety_interlock": safety_interlock(),
           "abstention_risk_coverage": risk_coverage(preds)}
    # variant-side calibration from Track 1 (if present)
    t1 = os.path.join(HERE, "track1_metrics.json")
    if os.path.exists(t1):
        c = json.load(open(t1, encoding="utf-8"))["calibration"]
        out["calibration_variant_from_track1"] = {
            "DISCERN_isotonic_ece": c["DISCERN_isotonic_oof"]["ece"],
            "REVEL_raw_ece": c["REVEL_raw_score"]["ece"],
            "AlphaMissense_raw_ece": c["AlphaMissense_raw_score"]["ece"]}
    with open(RC_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["coverage", "accuracy_on_retained", "confidence_threshold"])
        w.writerows(out["abstention_risk_coverage"]["curve"])
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    return out


def main():
    out = run()
    cd = out["calibration_diagnosis"]
    print("== TRACK 3: trustworthiness ==")
    print(f"\n3a calibration (diagnosis, n={cd['n']}): ECE={cd['ece']:.3f}  accuracy={cd['accuracy']:.0%}  "
          f"confidently-wrong(>=0.8)={cd['confidently_wrong_at_0.8']} ({cd['confidently_wrong_rate']:.0%})")
    if "calibration_variant_from_track1" in out:
        v = out["calibration_variant_from_track1"]
        print(f"   variant calibration (Track1): DISCERN ECE={v['DISCERN_isotonic_ece']} | "
              f"REVEL {v['REVEL_raw_ece']} | AlphaMissense {v['AlphaMissense_raw_ece']}")
    s = out["safety_interlock"]
    print(f"\n3b safety interlock (n={s['n_scenarios']} treatment-divergent scenarios):")
    print(f"   hard-stop SENSITIVITY = {s['hardstop_sensitivity']:.0%} (fired {s['fired_when_contraindicated']} when contraindicated)")
    print(f"   hard-stop SPECIFICITY = {s['hardstop_specificity']:.0%} (silent {s['silent_when_harmless']} when tx harmless)")
    rc = out["abstention_risk_coverage"]
    op = rc["operating_point_decide"]
    print(f"\n3c abstention / risk-coverage (n={cd['n']} cases):")
    print(f"   full-coverage accuracy={rc['accuracy_full_coverage']:.0%} -> top-50%-confident accuracy={rc['accuracy_top50pct_confident']:.0%}  "
          f"(monotone payoff: {rc['monotone_payoff']})")
    print(f"   DISCERN decide() operating point: coverage={op['coverage']:.0%}  accuracy-on-decided="
          f"{op['accuracy_on_decided']:.0%}  abstention={op['abstention_rate']:.0%}")
    print(f"\nwrote {os.path.relpath(OUT_JSON)} + {os.path.relpath(RC_CSV)}")


if __name__ == "__main__":
    main()
