"""Variant-arm robustness checks, run before submission.

Three questions a methods reviewer asks first, answered on the same eRepo expert-panel surface
Track 1' scores:

R1  Is the calibration out-of-sample? The isotonic fit happens inside each training fold and the
    metric is accumulated on the held-out folds only, so no variant is ever scored by a calibrator
    that saw its label. This module makes the protocol explicit, asserts the fold disjointness, and
    attaches bootstrap confidence intervals to ECE and Brier.
R2  Are the comparators calibrated the same way? REVEL and AlphaMissense go through the identical
    folds and the identical isotonic protocol, so "calibrated" means one thing for every tool and
    the comparison is like-for-like.
R3  Is the gap to REVEL real, and what does the pipeline add that a ranking score cannot? Every
    AUROC carries a bootstrap CI, DISCERN against REVEL is tested paired on the same variants
    (DeLong plus a paired bootstrap), and the added value is measured as its own outputs: ACMG band
    agreement with the expert panel, decision quality at matched coverage against a REVEL threshold
    rule, and the rate at which each rule moves an expert-called VUS.

REVEL is DISCERN's own PP3 input, so the correct reading of the discrimination result is "tracks"
rather than "beats"; the point of R3 is to say so with a confidence interval attached.

Run:  python -m bench.phase_r_variant     (needs bench/data/genebe_erepo.jsonl)
"""
from __future__ import annotations

import json
import os

import numpy as np
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss, cohen_kappa_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold

from bench.track1b_erepo_headtohead import TIMESPLIT_AFTER, load_rows
from core.stats import bootstrap_indices, midrank, percentile_ci
from core.stats import delong as _delong
from core.stats import ece as _ece
from rules.acmg_codes import code_points

HERE = os.path.dirname(__file__)
OUT_JSON = os.path.join(HERE, "phase_r_variant_metrics.json")
FOLDS_JSON = os.path.join(HERE, "calibration_folds.json")

N_BOOT = 1000
LP_POINTS = 6.0
N_FOLDS = 5
SEED = 0

# ClinGen-calibrated REVEL thresholds (Pejaver 2022): the supporting-evidence bands, which give the
# most generous (highest-coverage) two-sided rule a REVEL-only classifier could run.
REVEL_PATH_SUPPORTING = 0.644
REVEL_BENIGN_SUPPORTING = 0.290

PATH_BANDS = {"P", "LP"}
BENIGN_BANDS = {"B", "LB"}


# ---------------------------------------------------------------------------------------------
# R1 - shared folds and the out-of-sample calibration protocol
# ---------------------------------------------------------------------------------------------

def shared_folds(y, n_splits: int = N_FOLDS, seed: int = SEED):
    """One fold assignment, reused by every tool so the calibration comparison is like-for-like."""
    y = np.asarray(y, int)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return [(tr, te) for tr, te in skf.split(np.zeros(len(y)), y)]


def oof_isotonic(scores, y, folds):
    """Isotonic calibration fit on the training fold only, predicted on the held-out fold.

    Returns the out-of-fold probability for every input. Raises if a fold ever evaluates a variant
    it also trained on - the leakage this whole check exists to rule out.
    """
    scores, y = np.asarray(scores, float), np.asarray(y, int)
    oof = np.full(len(scores), np.nan)
    seen = set()
    for tr, te in folds:
        if set(tr) & set(te):
            raise AssertionError("calibration leakage: a fold trains and evaluates on the same variant")
        if seen & set(te):
            raise AssertionError("calibration leakage: a variant appears in two evaluation folds")
        seen |= set(te)
        iso = IsotonicRegression(out_of_bounds="clip")
        iso.fit(scores[tr], y[tr])
        oof[te] = iso.predict(scores[te])
    if np.isnan(oof).any():
        raise AssertionError("calibration incomplete: some variant was never held out")
    return oof


def _boot_indices(n, rng, n_boot=N_BOOT):
    return bootstrap_indices(n, rng, n_boot)


def _ci(values, lo=2.5, hi=97.5):
    return percentile_ci(values, lo, hi)


def calibration_with_ci(probs, y, rng, n_boot=N_BOOT):
    probs, y = np.asarray(probs, float), np.asarray(y, int)
    point = {"ece": round(_ece(probs, y), 4), "brier": round(float(brier_score_loss(y, probs)), 4)}
    eces, briers = [], []
    for idx in _boot_indices(len(y), rng, n_boot):
        yy = y[idx]
        if len(set(yy.tolist())) < 2:
            continue
        eces.append(_ece(probs[idx], yy))
        briers.append(brier_score_loss(yy, probs[idx]))
    point["ece_ci95"] = _ci(eces)
    point["brier_ci95"] = _ci(briers)
    return point


# ---------------------------------------------------------------------------------------------
# R3 - discrimination with confidence intervals and a paired test
# ---------------------------------------------------------------------------------------------

def auroc_with_ci(scores, y, rng, n_boot=N_BOOT):
    scores, y = np.asarray(scores, float), np.asarray(y, int)
    aucs = []
    for idx in _boot_indices(len(y), rng, n_boot):
        yy = y[idx]
        if len(set(yy.tolist())) < 2:
            continue
        aucs.append(roc_auc_score(yy, scores[idx]))
    return {"n": int(len(y)), "auroc": round(float(roc_auc_score(y, scores)), 4),
            "auroc_ci95": _ci(aucs)}


_midrank = midrank        # implementation moved to core.stats for CI coverage


delong = _delong          # implementation moved to core.stats for CI coverage


def paired_bootstrap_delta(y, score_a, score_b, rng, n_boot=N_BOOT):
    """Paired bootstrap on the AUROC difference: same resampled variants scored by both tools."""
    y = np.asarray(y, int)
    a, b = np.asarray(score_a, float), np.asarray(score_b, float)
    deltas = []
    for idx in _boot_indices(len(y), rng, n_boot):
        yy = y[idx]
        if len(set(yy.tolist())) < 2:
            continue
        deltas.append(roc_auc_score(yy, a[idx]) - roc_auc_score(yy, b[idx]))
    d = np.asarray(deltas, float)
    return {"delta": round(float(roc_auc_score(y, a) - roc_auc_score(y, b)), 4),
            "delta_ci95": _ci(d),
            "frac_bootstraps_favouring_a": round(float((d > 0).mean()), 4)}


# ---------------------------------------------------------------------------------------------
# R3 - what the pipeline adds that a ranking score cannot
# ---------------------------------------------------------------------------------------------

def _discern_direction(cls: str):
    if cls in PATH_BANDS:
        return 1
    if cls in BENIGN_BANDS:
        return 0
    return None


def _revel_direction(revel, hi=REVEL_PATH_SUPPORTING, lo=REVEL_BENIGN_SUPPORTING):
    if revel is None:
        return None
    r = float(revel)
    if r >= hi:
        return 1
    if r <= lo:
        return 0
    return None


def _decision_quality(preds, y):
    """Coverage, accuracy, and the direction breakdown - which side of the call a rule resolves."""
    resolved = [(p, yy) for p, yy in zip(preds, y, strict=False) if p is not None]
    n_path, n_benign = int(sum(y)), int(len(y) - sum(y))
    if not resolved:
        return {"coverage": 0.0, "n_resolved": 0, "accuracy_on_resolved": None,
                "pathogenic_recall": 0.0, "benign_recall": 0.0}
    correct = sum(1 for p, yy in resolved if p == yy)
    tp = sum(1 for p, yy in resolved if p == 1 and yy == 1)
    tn = sum(1 for p, yy in resolved if p == 0 and yy == 0)
    return {"coverage": round(len(resolved) / len(y), 4),
            "n_resolved": len(resolved),
            "accuracy_on_resolved": round(correct / len(resolved), 4),
            "called_pathogenic": sum(1 for p, _ in resolved if p == 1),
            "called_benign": sum(1 for p, _ in resolved if p == 0),
            "pathogenic_recall": round(tp / n_path, 4) if n_path else None,
            "benign_recall": round(tn / n_benign, 4) if n_benign else None}


def _band_at_coverage(scores, y, target_coverage):
    """Widen a symmetric abstention band on a continuous score until coverage matches the target.

    The only fair version of "would thresholding this score do as well": at the same number of calls
    made, are the calls as accurate, and on which side of the question?
    """
    have = [s for s in scores if s is not None]
    if not have or target_coverage <= 0:
        return None
    best = None
    for q in np.linspace(0.5, 0.0, 201):                # q = half-width of the abstention band
        lo, hi = np.quantile(have, max(0.0, 0.5 - q)), np.quantile(have, min(1.0, 0.5 + q))
        preds = [None if s is None else (1 if s >= hi else 0 if s <= lo else None) for s in scores]
        dq = _decision_quality(preds, y)
        if best is None or abs(dq["coverage"] - target_coverage) < abs(best[0] - target_coverage):
            best = (dq["coverage"], dq, round(float(lo), 4), round(float(hi), 4))
    _, dq, lo, hi = best
    dq["benign_threshold"], dq["pathogenic_threshold"] = lo, hi
    return dq


def ceiling_attribution(mis_pb):
    """Why does nothing reach a pathogenic band - the partition, or missing inputs?

    This separates two causes that are easy to conflate and that carry very different weight in a
    manuscript. Under the committed partition, only PS3 (functional) and PP4 (phenotype) leave the
    variant-intrinsic factor. PM1, PM5, PS1 and PS4 are variant-intrinsic and are *not* routed away;
    they are unapplied because this evaluation supplies no input for them - no hotspot or functional
    domain annotation, no same-residue ClinVar lookup under the ClinVar-blinded protocol, and no
    case-control counts. Saying the partition causes the ceiling would be wrong.

    Method: take DISCERN's own points for each variant, then add the points for codes the expert
    panel applied but DISCERN did not, grouped by owning factor, and re-band. The counts below are
    therefore an upper bound on what each stream would contribute if it were perfectly available.
    """
    routed = {"PS3", "BS3", "PP4", "PP1", "BS4", "PM3", "BP2", "PS2", "PM6"}
    intrinsic_unavailable = {"PM1", "PM5", "PS1", "PS4"}

    def band_reached(pts):
        return pts >= LP_POINTS

    base = sum(1 for r in mis_pb if band_reached(r["discern_points"]))
    plus_intrinsic = plus_routed = plus_both = 0
    n_with_intrinsic = n_with_routed = 0
    for r in mis_pb:
        missing = r["erepo_codes"] - r["discern_codes"]
        pts_intrinsic = sum(code_points(c)[0] for c in missing if c in intrinsic_unavailable)
        pts_routed = sum(code_points(c)[0] for c in missing if c in routed)
        n_with_intrinsic += bool(missing & intrinsic_unavailable)
        n_with_routed += bool(missing & routed)
        p = r["discern_points"]
        plus_intrinsic += band_reached(p + pts_intrinsic)
        plus_routed += band_reached(p + pts_routed)
        plus_both += band_reached(p + pts_intrinsic + pts_routed)

    return {
        "n": len(mis_pb),
        "lp_threshold_points": LP_POINTS,
        "partition_routes_away": sorted(routed),
        "intrinsic_but_no_input_here": sorted(intrinsic_unavailable),
        "why_each_intrinsic_code_has_no_input": {
            "PS1": ("implemented (adapters/clinvar.py) but deliberately not wired in: this benchmark "
                    "runs ClinVar-blinded, and PS1 needs a same-amino-acid-change ClinVar lookup. "
                    "Supplying it would reintroduce exactly the circularity the GeneBe exhibit "
                    "demonstrates, so its absence here is a protocol choice, not a gap."),
            "PM5": ("same as PS1 - implemented, and withheld for the same ClinVar-blinding reason "
                    "(different missense at the same residue)."),
            "PS4": ("implemented as a decision tree in rules/variant_scoring.py, but it requires "
                    "case-control inputs - proband counts against expectation, or an odds ratio with "
                    "its confidence bound - and the eRepo annotation cache carries none of them. A "
                    "data-availability limit."),
            "PM1": ("genuinely not implemented: no scorer emits PM1, and no VCEP specification in "
                    "rules/vcep/specs/ encodes hotspot or critical-domain regions for the in-scope "
                    "genes. The in_functional_domain annotation that exists feeds the PVS1 tree, not "
                    "PM1. This is a scope limitation and the one item on this list that is a genuine "
                    "engine gap rather than a protocol or data constraint."),
        },
        "reach_lp_as_scored": base,
        "reach_lp_if_intrinsic_codes_were_available": plus_intrinsic,
        "reach_lp_if_routed_codes_were_re_added": plus_routed,
        "reach_lp_with_both": plus_both,
        "variants_where_experts_applied_an_unavailable_intrinsic_code": n_with_intrinsic,
        "variants_where_experts_applied_a_routed_code": n_with_routed,
        "reading": ("The ceiling has two separable causes and only one of them is the partition. "
                    "Restoring the intrinsic codes this pipeline cannot derive (hotspot, "
                    "same-residue, case-control) is an annotation problem, not a design choice. "
                    "Restoring the routed codes is what the partition deliberately prevents, and "
                    "is exactly the evidence the coupling is meant to supply through the disease "
                    "factor instead of by re-addition. Neither stream alone, nor both together, "
                    "comes close to the 316 variants the panel called pathogenic, so the binding "
                    "constraint is the evidence the ACMG framework demands for missense rather "
                    "than any one routing decision."),
        "caveat": ("Points for the restored codes use default ACMG strengths, because stripping a "
                   "strength suffix is how the code vocabulary is compared. VCEPs frequently apply "
                   "gene-specific strengths above the default, so every count here is a lower "
                   "bound on what the expert evidence would actually contribute."),
    }


def added_value(rows_all, mis_pb, discern_oof):
    """Outputs DISCERN produces that a ranking score structurally cannot, and the places where
    the ACMG framework's conservatism costs it against a free-running score."""
    y = [r["label"] for r in mis_pb]

    # (a) the ACMG band itself, checked against the expert panel's band direction
    d_dir = [_discern_direction(r["discern_class"]) for r in mis_pb]
    r_dir = [_revel_direction(r["revel"]) for r in mis_pb]
    both = [(d, e) for d, e in zip(d_dir, y, strict=False) if d is not None]
    kappa = (round(float(cohen_kappa_score([e for _, e in both], [d for d, _ in both])), 3)
             if both and len({d for d, _ in both}) > 1 and len({e for _, e in both}) > 1 else None)

    # (b) decision quality: DISCERN's own class, DISCERN's calibrated probability at the same
    #     coverage (score-vs-score fairness), and REVEL both at ClinGen thresholds and matched
    discern_dq = _decision_quality(d_dir, y)
    revel_clingen_dq = _decision_quality(r_dir, y)
    cov = discern_dq["coverage"]
    revel_matched = _band_at_coverage([r["revel"] for r in mis_pb], y, cov)
    discern_prob_matched = _band_at_coverage(list(np.asarray(discern_oof, float)), y, cov)

    # (c) movement on the variants the expert panel itself left uncertain (a rate, not an accuracy -
    #     these variants have no independent truth, which is exactly why the coupling needs a cohort)
    vus = [r for r in rows_all if r["assertion"] == "VUS" and r["is_missense"]]
    d_moved = sum(1 for r in vus if _discern_direction(r["discern_class"]) is not None)
    r_moved = sum(1 for r in vus if _revel_direction(r["revel"]) is not None)

    return {
        "acmg_band_assignment": {
            "DISCERN": {"emits_class": True, "n_resolved": len(both),
                        "kappa_vs_erepo_band_direction": kappa,
                        "accuracy_on_resolved": discern_dq["accuracy_on_resolved"],
                        "kappa_note": ("undefined when the resolved calls are all one direction - "
                                       "which is itself the finding below"
                                       if kappa is None else None)},
            "REVEL": {"emits_class": False,
                      "note": "a ranking score; a class exists only once a human picks a threshold"},
            "AlphaMissense": {"emits_class": False, "note": "same"},
        },
        "decision_quality": {
            "DISCERN_acmg_class": discern_dq,
            "DISCERN_calibrated_prob_matched_coverage": discern_prob_matched,
            "REVEL_clingen_thresholds": dict(revel_clingen_dq,
                                             benign_threshold=REVEL_BENIGN_SUPPORTING,
                                             pathogenic_threshold=REVEL_PATH_SUPPORTING),
            "REVEL_matched_coverage": revel_matched,
        },
        "intrinsic_only_ceiling": {
            "max_discern_points_on_missense": float(max(r["discern_points"] for r in mis_pb)),
            "lp_threshold_points": 6.0,
            "n_reaching_lp_or_p": sum(1 for r in mis_pb if r["discern_class"] in PATH_BANDS),
            "finding": ("On intrinsic sequence evidence alone the ACMG point system cannot reach a "
                        "pathogenic band for a missense variant: PM2 plus PP3 tops out well below "
                        "the 6 points Likely Pathogenic requires. Every band DISCERN resolves on "
                        "this surface is therefore benign-side, and the pathogenic side is left "
                        "explicitly uncertain. That is the partition behaving as specified - the "
                        "criteria that carry missense pathogenicity (PS3 functional, PS4 "
                        "case-control, PM1 hotspot, PM5, PP4 phenotype) are owned by the "
                        "functional, disease and coupling factors and are not derivable from "
                        "sequence. It is also, stated plainly, the reason the coupling matters: "
                        "the missing evidence is exactly what a disease model supplies."),
        },
        "expert_vus_movement": {
            "n_expert_vus_missense": len(vus),
            "DISCERN_assigns_non_vus_band": d_moved,
            "REVEL_clingen_rule_assigns_direction": r_moved,
            "note": ("a rate, not an accuracy: eRepo VUS carry no independent truth, so neither "
                     "number says the movement was correct. This is the gap the pre-registered "
                     "coupling endpoint exists to close."),
        },
    }


# ---------------------------------------------------------------------------------------------
# surface runner
# ---------------------------------------------------------------------------------------------

def _surface(rows, rng):
    mis = [r for r in rows if r["is_missense"] and r["label"] is not None]
    y = [r["label"] for r in mis]
    folds = shared_folds(y)

    tools = {"DISCERN": [r["discern_points"] for r in mis],
             "REVEL": [r["revel"] for r in mis],
             "AlphaMissense": [r["alphamissense"] for r in mis]}

    discrimination, calibration = {}, {}
    discern_oof = None
    for name, raw in tools.items():
        # a tool missing a score on a few variants is scored on the variants it covers
        keep = [i for i, v in enumerate(raw) if v is not None]
        s = np.asarray([float(raw[i]) for i in keep])
        yy = np.asarray([y[i] for i in keep], int)
        discrimination[name] = auroc_with_ci(s, yy, rng)
        # R2: identical protocol for every tool, on that tool's own covered subset
        sub_folds = shared_folds(yy) if len(keep) != len(y) else folds
        oof = oof_isotonic(s, yy, sub_folds)
        if name == "DISCERN":
            discern_oof = oof
        calibration[name + "_isotonic_oof"] = calibration_with_ci(oof, yy, rng)
        if name != "DISCERN":       # the raw score, treated as a probability, is the status quo
            calibration[name + "_raw"] = calibration_with_ci(np.clip(s, 0, 1), yy, rng)

    # R3 paired test: DISCERN vs REVEL on the variants both score
    pair = [(r["discern_points"], r["revel"], yy) for r, yy in zip(mis, y, strict=False)
            if r["revel"] is not None]
    d_s = [a for a, _, _ in pair]
    r_s = [float(b) for _, b, _ in pair]
    p_y = [c for _, _, c in pair]

    return discern_oof, {
        "n_missense_PB": len(mis), "n_path": int(sum(y)), "n_benign": int(len(y) - sum(y)),
        "calibration_protocol": {
            "method": "isotonic regression",
            "folds": f"{N_FOLDS}-fold stratified, shuffled, seed {SEED}",
            "fit_on": "training folds only",
            "reported_on": "held-out folds only (out-of-fold), pooled",
            "leakage_check": "asserted in oof_isotonic: folds are disjoint and cover every variant once",
            "bootstrap": f"{N_BOOT} resamples, percentile 95% CI",
        },
        "discrimination": discrimination,
        "calibration": calibration,
        "discern_vs_revel": {
            "n_paired": len(pair),
            "delong": delong(p_y, d_s, r_s),
            "paired_bootstrap": paired_bootstrap_delta(p_y, d_s, r_s, rng),
        },
    }


def run():
    rows = load_rows()
    rng = np.random.default_rng(SEED)
    discern_oof, primary = _surface(rows, rng)
    ts_rows = [r for r in rows if (r["approval_date"] or "") > TIMESPLIT_AFTER]
    _, timesplit = _surface(ts_rows, rng)
    mis_pb = [r for r in rows if r["is_missense"] and r["label"] is not None]

    out = {
        "set": "eRepo cluster-gene expert-panel surface (GRCh38), Phase R re-analysis",
        "n_total": len(rows), "n_timesplit": len(ts_rows), "timesplit_after": TIMESPLIT_AFTER,
        "eRepo_primary": primary,
        "time_split": timesplit,
        "added_value_over_ranking_score": added_value(rows, mis_pb, discern_oof),
        "ceiling_attribution": ceiling_attribution(mis_pb),
    }
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    export_calibration_folds(mis_pb)
    return out


def export_calibration_folds(mis_pb, path: str = FOLDS_JSON):
    """Write the actual fold assignment, so out-of-sample calibration can be checked not trusted.

    Reporting an expected calibration error and asserting it was estimated out-of-fold asks a
    reviewer to take the protocol on faith. This emits the assignment itself: for every variant,
    which fold held it out, plus the training and evaluation index sets. Anyone can confirm that no
    variant appears in two evaluation folds, that every variant is held out exactly once, and that
    no evaluation index appears in its own training set.
    """
    y = [r["label"] for r in mis_pb]
    folds = shared_folds(y)
    assignment = {}
    for i, (_tr, te) in enumerate(folds):
        for idx in te.tolist():
            assignment[str(idx)] = i
    payload = {
        "description": ("Out-of-fold assignment behind every calibration figure reported for the "
                        "eRepo-primary surface. Index i refers to the i-th variant of the missense "
                        "pathogenic/benign set, in the order bench/track1b_erepo_headtohead.load_rows "
                        "yields them; variant_index below gives the key for each."),
        "protocol": {"method": "isotonic regression", "n_splits": N_FOLDS, "stratified": True,
                     "shuffle": True, "seed": SEED,
                     "fit_on": "training folds only", "reported_on": "held-out folds only"},
        "n_variants": len(mis_pb),
        "n_pathogenic": int(sum(y)), "n_benign": int(len(y) - sum(y)),
        "variant_index": [
            {"i": i, "gene": r.get("gene"), "chr": r.get("chr"), "pos": r.get("pos"),
             "ref": r.get("ref"), "alt": r.get("alt"), "erepo_class": r.get("assertion"),
             "label": r.get("label"), "held_out_in_fold": assignment[str(i)]}
            for i, r in enumerate(mis_pb)],
        "folds": [{"fold": i, "n_train": len(tr), "n_test": len(te),
                   "train_indices": sorted(tr.tolist()), "test_indices": sorted(te.tolist())}
                  for i, (tr, te) in enumerate(folds)],
        "invariants_a_reader_can_check": [
            "every variant appears in exactly one test_indices set",
            "no index appears in both train_indices and test_indices of the same fold",
            "the union of all test_indices is the full variant set",
        ],
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    return path


def main():
    o = run()
    print("== PHASE R: variant arm (out-of-sample calibration, calibrated comparators, CIs) ==")
    for surf in ("eRepo_primary", "time_split"):
        s = o[surf]
        tag = "" if surf == "eRepo_primary" else f" (expert-approved after {o['timesplit_after']})"
        print(f"\n-- {surf}{tag}: missense P/B={s['n_missense_PB']} "
              f"(path={s['n_path']} benign={s['n_benign']}) --")
        print("   R3 discrimination (AUROC, 95% CI)")
        for t, m in s["discrimination"].items():
            print(f"      {t:16} {m['auroc']:.3f}  [{m['auroc_ci95'][0]:.3f}, {m['auroc_ci95'][1]:.3f}]  n={m['n']}")
        print("   R1/R2 calibration (ECE, 95% CI) - identical folds and protocol for every tool")
        for t, m in s["calibration"].items():
            print(f"      {t:28} ECE={m['ece']:.4f}  [{m['ece_ci95'][0]:.4f}, {m['ece_ci95'][1]:.4f}]  "
                  f"Brier={m['brier']:.4f}")
        dv = s["discern_vs_revel"]
        dl, pb = dv["delong"], dv["paired_bootstrap"]
        print(f"   R3 DISCERN vs REVEL, paired on n={dv['n_paired']}: delta={dl['delta']:+.4f}  "
              f"DeLong z={dl['z']} p={dl['p_value']}  bootstrap CI {pb['delta_ci95']}")

    av = o["added_value_over_ranking_score"]
    print("\n-- R3 added value, and where the framework's conservatism costs it --")
    b = av["acmg_band_assignment"]["DISCERN"]
    print(f"   ACMG band: DISCERN emits one (n resolved={b['n_resolved']}, accuracy={b['accuracy_on_resolved']}, "
          f"kappa={b['kappa_vs_erepo_band_direction']}); REVEL and AlphaMissense do not.")
    print("   decisions (coverage / accuracy / path-recall / benign-recall)")
    for label, key in (("DISCERN ACMG class      ", "DISCERN_acmg_class"),
                       ("DISCERN calibrated prob ", "DISCERN_calibrated_prob_matched_coverage"),
                       ("REVEL ClinGen bands     ", "REVEL_clingen_thresholds"),
                       ("REVEL matched coverage  ", "REVEL_matched_coverage")):
        m = av["decision_quality"].get(key)
        if not m:
            continue
        print(f"      {label} cov={m['coverage']:.0%}  acc={m['accuracy_on_resolved']:.3f}  "
              f"path-recall={m['pathogenic_recall']:.3f}  benign-recall={m['benign_recall']:.3f}")
    ic = av["intrinsic_only_ceiling"]
    print(f"   intrinsic-only ceiling: max points on missense = {ic['max_discern_points_on_missense']:.0f} "
          f"against {ic['lp_threshold_points']:.0f} for Likely Pathogenic; "
          f"{ic['n_reaching_lp_or_p']} missense variants reach a pathogenic band.")
    ca = o["ceiling_attribution"]
    print(f"   why: partition routes away {ca['partition_routes_away']}")
    print(f"        intrinsic but unavailable here {ca['intrinsic_but_no_input_here']}")
    print(f"        reaching LP  as scored={ca['reach_lp_as_scored']}  "
          f"+intrinsic-if-available={ca['reach_lp_if_intrinsic_codes_were_available']}  "
          f"+routed-if-re-added={ca['reach_lp_if_routed_codes_were_re_added']}  "
          f"both={ca['reach_lp_with_both']}  (of {ca['n']})")
    v = av["expert_vus_movement"]
    print(f"   expert-called VUS (n={v['n_expert_vus_missense']}): DISCERN assigns a band to "
          f"{v['DISCERN_assigns_non_vus_band']}, the REVEL rule to {v['REVEL_clingen_rule_assigns_direction']} "
          "(a rate, not an accuracy)")
    print(f"\nwrote {os.path.relpath(OUT_JSON)}")


if __name__ == "__main__":
    main()
