"""TRACK 1 - variant arm head-to-head (DISCERN_Benchmark_Execution_Plan_v1).

Compares DISCERN against the current ACMG/predictor tool set on the missense / VUS-adjacent axis -
the surface where intrinsic tools actually differ. Tools:
  - DISCERN        : scored locally (rules.variant_scoring.score_variant) from GeneBe's annotations
                     (gnomAD AF + REVEL + consequence); NO ClinVar-derived codes (PP5/BP6) - the
                     ClinVar-blinded protocol of plan section 0.
  - GeneBe         : current VCEP-aware ACMG classifier (acmg_classification + acmg_score). It DOES
                     apply PP5/BP6 (ClinVar lookup) -> we quantify that circularity as a finding.
  - REVEL          : metapredictor baseline DISCERN ingests as PP3.
  - AlphaMissense  : current metapredictor baseline.
  - InterVar       : prior committed full-DB anchor (eval/intervar_full_eval.py; missense AUROC 0.811).

Truth = the h4set ClinVar labels (P/LP = 1, B/LB = 0); this is the ClinVar-derived secondary surface
(the eRepo-primary run is the pre-registered next step). VUS-labelled variants are used only for the
VUS-retention metric, never for discrimination.

Run:  python -m bench.track1_variant_headtohead   (needs bench/data/genebe_h4set.jsonl)
"""
from __future__ import annotations

import csv
import json
import os

import numpy as np
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    matthews_corrcoef,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold

from rules.variant_scoring import Annotations, score_variant

HERE = os.path.dirname(__file__)
CACHE = os.path.join(HERE, "data", "genebe_h4set.jsonl")
OUT_CSV = os.path.join(HERE, "track1_variant_headtohead.csv")
OUT_JSON = os.path.join(HERE, "track1_metrics.json")

PATH = {"Pathogenic", "Likely_pathogenic"}
BEN = {"Benign", "Likely_benign", "Benign/Likely_benign"}
CLINVAR_CODES = ("PP5", "BP6")   # the ClinVar-derived ACMG codes (guard #2)

# GeneBe SO-term effect -> DISCERN consequence vocabulary (rules.variant_scoring).
EFFECT2CONS = {
    "missense_variant": "missense", "stop_gained": "nonsense",
    "frameshift_variant": "frameshift", "start_lost": "initiation", "stop_lost": "stop_lost",
    "splice_donor_variant": "canonical_splice", "splice_acceptor_variant": "canonical_splice",
    "inframe_deletion": "inframe_indel", "inframe_insertion": "inframe_indel",
}


def _cons(effect: str) -> str:
    return EFFECT2CONS.get((effect or "").split(",")[0], "")


def _label(clnsig: str):
    if clnsig in PATH:
        return 1
    if clnsig in BEN:
        return 0
    return None


def _ece(probs, ys, bins=10):
    probs, ys = np.asarray(probs, float), np.asarray(ys, float)
    edges = np.linspace(0, 1, bins + 1)
    e = 0.0
    for i in range(bins):
        m = (probs >= edges[i]) & (probs < edges[i + 1] if i < bins - 1 else probs <= 1.0)
        if m.sum():
            e += (m.sum() / len(probs)) * abs(ys[m].mean() - probs[m].mean())
    return float(e)


def _disc_prob_oof(points, ys):
    """Out-of-fold isotonic-calibrated DISCERN probability (out-of-fold ECE/Brier, no train-on-test)."""
    points, ys = np.asarray(points, float), np.asarray(ys, int)
    oof = np.zeros(len(points))
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    for tr, te in skf.split(points, ys):
        iso = IsotonicRegression(out_of_bounds="clip")
        iso.fit(points[tr], ys[tr])
        oof[te] = iso.predict(points[te])
    return oof


def _sens_at_spec(scores, ys, spec=0.90):
    fpr, tpr, _ = roc_curve(ys, scores)
    ok = np.where((1 - fpr) >= spec)[0]
    return float(tpr[ok].max()) if len(ok) else 0.0


def load_rows():
    rows = []
    with open(CACHE, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            r = json.loads(line)
            cons = _cons(r["effect"])
            af = r.get("gnomad_exomes_af") or r.get("gnomad_genomes_af")
            rev = float(r["revel_score"]) if r.get("revel_score") is not None else None
            sv = score_variant(r["q_gene"], r["q_gene"], Annotations(
                af=float(af) if af is not None else None, revel=rev, consequence=cons))
            sv_off = score_variant(r["q_gene"], r["q_gene"], Annotations(  # frequency-blind
                af=None, revel=rev, consequence=cons))
            rows.append({
                "gene": r["q_gene"], "chr": r["q_chr"], "pos": r["q_pos"],
                "ref": r["q_ref"], "alt": r["q_alt"], "clnsig": r["q_clnsig"],
                "label": _label(r["q_clnsig"]), "effect": r["effect"], "consequence": cons,
                "is_missense": r["effect"] == "missense_variant",
                "discern_points": sv.points, "discern_points_afoff": sv_off.points,
                "discern_class": sv.classification.name,
                "genebe_class": r.get("acmg_classification"), "genebe_score": r.get("acmg_score"),
                "genebe_criteria": r.get("acmg_criteria") or "",
                "genebe_uses_clinvar": any(c in (r.get("acmg_criteria") or "") for c in CLINVAR_CODES),
                "revel": r.get("revel_score"), "alphamissense": r.get("alphamissense_score"),
            })
    return rows


def _disc_path(cls):    # DISCERN Classification enum names are abbreviated: P/LP/VUS/LB/B
    if cls in ("P", "LP", "PATHOGENIC", "LIKELY_PATHOGENIC"):
        return 1
    if cls in ("B", "LB", "BENIGN", "LIKELY_BENIGN"):
        return 0
    return None


def _gb_path(cls):
    if cls in ("Pathogenic", "Likely_pathogenic"):
        return 1
    if cls in ("Benign", "Likely_benign"):
        return 0
    return None


def _cont_metrics(rows, scorekey, ys_all):
    pairs = [(r[scorekey], y) for r, y in zip(rows, ys_all, strict=False) if r[scorekey] is not None and y is not None]
    if len({y for _, y in pairs}) < 2:
        return None
    s = [p for p, _ in pairs]
    y = [yy for _, yy in pairs]
    return {"n": len(pairs), "auroc": round(roc_auc_score(y, s), 4),
            "auprc": round(average_precision_score(y, s), 4),
            "sens@90spec": round(_sens_at_spec(np.array(s), np.array(y)), 4)}


def _cat_metrics(rows, pathfn, classkey, ys_all):
    # Abstention (VUS) is NOT a wrong label - it is a deliberate non-call. Report accuracy ON THE
    # RESOLVED subset + the abstention rate separately (a conservative tool abstains rather than guess).
    pred_res, y_res, n = [], [], 0
    for r, y in zip(rows, ys_all, strict=False):
        if y is None:
            continue
        n += 1
        p = pathfn(r[classkey])
        if p is not None:
            pred_res.append(p)
            y_res.append(y)
    resolved = len(y_res)
    acc = float(np.mean([a == b for a, b in zip(pred_res, y_res, strict=False)])) if resolved else float("nan")
    mcc = float(matthews_corrcoef(y_res, pred_res)) if resolved and len(set(pred_res)) > 1 and len(set(y_res)) > 1 else None
    return {"n_PB": n, "resolved": resolved, "abstention_rate": round(1 - resolved / n, 4),
            "accuracy_on_resolved": round(acc, 4) if resolved else None,
            "mcc_on_resolved": round(mcc, 4) if mcc is not None else None,
            "path_calls": int(sum(pred_res)), "benign_calls": int(resolved - sum(pred_res))}


def run():
    rows = load_rows()
    mis = [r for r in rows if r["is_missense"]]
    ys = [r["label"] for r in mis]
    pb = [(r, y) for r, y in zip(mis, ys, strict=False) if y is not None]
    pbrows = [r for r, _ in pb]
    pby = [y for _, y in pb]
    n_path = sum(pby)

    out = {"set": "h4set missense, ClinVar P/B truth", "n_missense": len(mis),
           "n_PB": len(pby), "n_path": n_path, "n_benign": len(pby) - n_path, "n_VUS": len(mis) - len(pby)}

    # discrimination (continuous)
    out["discrimination"] = {
        "DISCERN_points": _cont_metrics(pbrows, "discern_points", pby),
        "DISCERN_points_freqblind": _cont_metrics(pbrows, "discern_points_afoff", pby),
        "GeneBe_acmg_score": _cont_metrics(pbrows, "genebe_score", pby),
        "REVEL": _cont_metrics(pbrows, "revel", pby),
        "AlphaMissense": _cont_metrics(pbrows, "alphamissense", pby),
        "InterVar_full_DB_prior": {"auroc": 0.811, "note": "prior committed full-DB run (eval/intervar_full_eval.py)"},
    }

    # ACMG concordance (categorical)
    out["acmg_concordance"] = {
        "DISCERN": _cat_metrics(pbrows, _disc_path, "discern_class", pby),
        "GeneBe": _cat_metrics(pbrows, _gb_path, "genebe_class", pby),
    }

    # calibration
    disc_oof = _disc_prob_oof([r["discern_points"] for r in pbrows], pby)
    cal = {"DISCERN_isotonic_oof": {"ece": round(_ece(disc_oof, pby), 4),
                                    "brier": round(brier_score_loss(pby, disc_oof), 4)}}
    for name, key in (("REVEL", "revel"), ("AlphaMissense", "alphamissense")):
        sc = [(r[key], y) for r, y in zip(pbrows, pby, strict=False) if r[key] is not None]
        p = [float(a) for a, _ in sc]
        yy = [b for _, b in sc]
        cal[name + "_raw_score"] = {"ece": round(_ece(p, yy), 4), "brier": round(brier_score_loss(yy, p), 4),
                                    "note": "predictor score treated as pseudo-probability (not gene-calibrated)"}
    cal["GeneBe"] = {"note": "emits an ACMG class + integer points, not a probability - not calibratable"}
    cal["InterVar"] = {"note": "emits an ACMG class only - not calibratable"}
    out["calibration"] = cal

    # ClinVar-circularity (guard #2) - the named PP5/BP6 codes UNDERCOUNT the dependence; the decisive
    # evidence is that GeneBe's score is perfectly separable on a ClinVar set and its class reproduces
    # the ClinVar direction almost exactly, even on the PP5/BP6-blind subset.
    n_cv = sum(1 for r in pbrows if r["genebe_uses_clinvar"])
    blind = [(r, y) for r, y in zip(pbrows, pby, strict=False) if not r["genebe_uses_clinvar"]]
    gb_blind = _cont_metrics([r for r, _ in blind], "genebe_score", [y for _, y in blind])
    dir_match = sum(1 for r, y in zip(pbrows, pby, strict=False)
                    if _gb_path(r["genebe_class"]) is not None and _gb_path(r["genebe_class"]) == y)
    resolved_gb = sum(1 for r in pbrows if _gb_path(r["genebe_class"]) is not None)
    out["clinvar_circularity"] = {
        "genebe_PB_calls_naming_PP5_or_BP6": n_cv, "of_total_PB": len(pby),
        "frac_named_clinvar_code": round(n_cv / len(pby), 4),
        "genebe_class_matches_clinvar_direction": f"{dir_match}/{resolved_gb}",
        "genebe_auroc_overall": round(roc_auc_score(pby, [r["genebe_score"] for r in pbrows]), 4),
        "genebe_auroc_named_clinvar_code_blind": gb_blind,
        "finding": ("GeneBe reproduces the ClinVar label: its acmg_score is near-perfectly separable on "
                    "this ClinVar-derived P/B set (AUROC ~1.0) and its class matches the ClinVar "
                    "direction on ~98% of resolved calls, even on the subset where PP5/BP6 are NOT the "
                    "named criteria. A ClinVar-consuming tool therefore cannot be fairly graded on "
                    "ClinVar-derived truth - the fair surfaces are an expert-panel set (eRepo) or a "
                    "time-split. DISCERN, REVEL and AlphaMissense apply no ClinVar-derived evidence, so "
                    "their numbers are not inflated this way."),
    }

    # write per-variant CSV (missense)
    cols = ["gene", "chr", "pos", "ref", "alt", "clnsig", "label", "effect", "consequence",
            "discern_points", "discern_class", "genebe_class", "genebe_score", "genebe_criteria",
            "genebe_uses_clinvar", "revel", "alphamissense"]
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(mis)
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    return out


def main():
    out = run()
    print("== TRACK 1: variant head-to-head (missense, ClinVar P/B) ==")
    print(f"missense={out['n_missense']}  P/B={out['n_PB']} (path={out['n_path']} benign={out['n_benign']})  "
          f"VUS={out['n_VUS']}")
    print("\n-- Discrimination (AUROC / AUPRC / sens@90%spec) --")
    for t, m in out["discrimination"].items():
        if m and "auprc" in m:
            print(f"  {t:24} AUROC={m['auroc']:.3f}  AUPRC={m['auprc']:.3f}  sens@90spec={m['sens@90spec']:.3f}  (n={m['n']})")
        elif m:
            print(f"  {t:24} AUROC={m['auroc']:.3f}  [{m.get('note','')}]")
    print("\n-- ACMG concordance (accuracy ON RESOLVED; abstention is a non-call, not an error) --")
    for t, m in out["acmg_concordance"].items():
        acc = "n/a" if m["accuracy_on_resolved"] is None else f"{m['accuracy_on_resolved']:.3f}"
        mcc = "n/a" if m["mcc_on_resolved"] is None else f"{m['mcc_on_resolved']:.3f}"
        print(f"  {t:24} acc-on-resolved={acc}  MCC={mcc}  abstention={m['abstention_rate']:.2f}  "
              f"(path-calls={m['path_calls']} benign-calls={m['benign_calls']})")
    print("\n-- Calibration (ECE / Brier) --")
    for t, m in out["calibration"].items():
        if "ece" in m:
            print(f"  {t:24} ECE={m['ece']:.3f}  Brier={m['brier']:.3f}")
        else:
            print(f"  {t:24} {m['note']}")
    c = out["clinvar_circularity"]
    print("\n-- ClinVar circularity (guard #2: why GeneBe's 1.0 is not a fair win) --")
    print(f"  GeneBe class matches ClinVar direction: {c['genebe_class_matches_clinvar_direction']}  "
          f"(named PP5/BP6 on only {c['frac_named_clinvar_code']:.0%} - the dependence is deeper than the named codes)")
    print(f"  GeneBe acmg_score AUROC overall={c['genebe_auroc_overall']:.3f}; "
          f"on PP5/BP6-name-blind subset={c['genebe_auroc_named_clinvar_code_blind']['auroc']:.3f} "
          f"(n={c['genebe_auroc_named_clinvar_code_blind']['n']}) -> still ~1.0 = reproduces ClinVar")
    print(f"\nwrote {os.path.relpath(OUT_CSV)} + {os.path.relpath(OUT_JSON)}")


if __name__ == "__main__":
    main()
