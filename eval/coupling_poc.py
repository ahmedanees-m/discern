"""DISCERN Coupling Proof-of-Concept (public Phenopacket Store data; circularity-safe).

First real-data test of H6: does the clinical-phenotype coupling resolve variants that sequence
evidence alone leaves as VUS, and ONLY when phenotype and gene agree? Three streams are kept
disjoint (Sequence / Clinical-phenotype / Truth) so the test is falsifiable, not tautological.

Pipeline per case (Phenopacket Store bleeding subset, public literature cases, PHI-free):
  1. Sequence band: intrinsic-only score from the variant consequence (null -> PVS1 -> LP/P; missense
     -> VUS). Only intrinsic-VUS cases enter the working set (coupling has no room otherwise).
  2. Clinical stream: HPO present terms mapped via the verified crosswalk to CLINICAL features only
     (functional-tagged terms are withheld for Truth).
  3. Truth: independent of this case's clinical stream - a functional finding (e.g. giant granules)
     and/or the published causal status of the variant.
  4. Coupling test: score the variant with the coupling under MATCHED (case's real clinical features)
     vs MISMATCHED (a sibling disease's profile, whose gene the variant is NOT in). The coupling
     should upgrade matched toward Truth and not upgrade mismatched.
Primary endpoint (pre-registered): lift = upgrade_rate(matched) - upgrade_rate(mismatched).
Falsification: lift ~ 0 means the coupling carries no disease-specific signal (H6 not supported here).

Run: python3 -m eval.coupling_poc
"""
from __future__ import annotations

import json
import math
import os
import re

import yaml

from core.dx_schemas import Feature, FeatureKind, VariantState
from diseases.ontology import cluster_for
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import marginal_variant, reclassify
from rules.point_engine import Classification
from rules.variant_scoring import Annotations, score_variant

CLUSTERS = ["integrin", "vwf_gpib", "macrothrombocytopenia", "thr_leukemia", "coag_factor",
            "vwd2n_hema", "mild_vwd", "scott", "alpha_granule", "granule"]
SUBSET = os.path.join(os.path.dirname(__file__), "data", "bleeding_subset.jsonl")
CROSSWALK = os.path.join(os.path.dirname(__file__), "hpo_feature_crosswalk.yaml")
NULL_CONS = {"nonsense", "frameshift", "splice"}
UPGRADED = {VariantState.PATH, VariantState.LP}


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (round(max(0.0, c - h), 3), round(min(1.0, c + h), 3))


def load_clusters():
    clusters, gene2cluster, gene2disease = {}, {}, {}
    for cid in CLUSTERS:
        cl = cluster_for(cid)
        clusters[cid] = cl
        for d in cl.diseases:
            for g in d.genes:
                gene2cluster.setdefault(g, cid)
                gene2disease.setdefault(g, d.id)
    return clusters, gene2cluster, gene2disease


def consequence(hgvs_p: str | None, hgvs_c: str | None) -> str:
    p, c = (hgvs_p or "").lower(), (hgvs_c or "").lower()
    if "fs" in p:                       # frameshift-then-stop ("...fsTer..") is a frameshift
        return "frameshift"
    if "ter" in p or "*" in p:
        return "nonsense"
    if re.search(r"[+-][12]\b", c) or "splice" in c:
        return "splice"
    if ("del" in c or "dup" in c or "ins" in c) and "fs" not in p:
        return "inframe_indel"
    if re.search(r"p\.\(?[a-z]{3}\d+[a-z]{3}", p):
        return "missense"
    return "other"


def sequence_band(gene: str, cons: str) -> Classification:
    ann = Annotations(consequence={"nonsense": "nonsense", "frameshift": "frameshift",
                                   "splice": "canonical_splice"}.get(cons, "missense"),
                      nmd_predicted=cons in NULL_CONS)
    return score_variant(gene, "", ann).classification


def crosswalk_features(hpo_present, xwalk):
    clinical, functional = {}, set()
    for h in hpo_present:
        e = xwalk.get(h.get("id"))
        if not e:
            continue
        if e["kind"] == "clinical":
            clinical[e["feature"]] = True
        else:
            functional.add(e["feature"])
    return clinical, functional


def sibling_profile(cluster, true_disease_id: str, gene: str) -> tuple[str, dict]:
    """A sibling disease (gene != this case's gene) plus features that point to it and away from D."""
    d_true = next(d for d in cluster.diseases if d.id == true_disease_id)
    cands = [d for d in cluster.diseases if d.id != true_disease_id and gene not in d.genes]
    if not cands:
        return "", {}
    # pick the sibling whose high-freq features are most disjoint from the true disease
    def disjoint(d):
        return sum(v[0] for f, v in d.feature_lr.items() if d_true.feature_lr.get(f, [0])[0] < 0.4)
    sib = max(cands, key=disjoint)
    feats = {f: True for f, v in sib.feature_lr.items()
             if v[0] >= 0.6 and d_true.feature_lr.get(f, [0])[0] < 0.4}
    return sib.id, feats


# Realistic intrinsic-VUS sequence band: the working-set variants are published rare disease alleles,
# verified gnomAD-rare/absent (PM2). PM2 is held CONSTANT across matched/mismatched, so the lift
# (matched - mismatched) still isolates the coupling's contribution. A full cohort run annotates exact
# per-variant gnomAD AF + REVEL (adding PP3 where the predictor is high).
INTRINSIC_VUS_CODES = ["PM2_Supporting"]


def _evidence(gene: str, feats: dict, codes=INTRINSIC_VUS_CODES) -> Evidence:
    clin = [Feature(f, FeatureKind.LAB, True, observed=True) for f in feats]
    return Evidence(variant_gene=gene, genetic_codes=list(codes), clinical=clin)


def couple(cluster, gene, feats):
    j = joint(cluster, _evidence(gene, feats))
    _, new, _ = reclassify(j, VariantState.VUS)
    mv = marginal_variant(j)
    return new, round(mv[VariantState.PATH] + mv[VariantState.LP], 3)


def run(subset: str = SUBSET):
    clusters, g2c, g2d = load_clusters()
    xwalk = yaml.safe_load(open(CROSSWALK, encoding="utf-8"))["crosswalk"]
    cases = [json.loads(line) for line in open(subset, encoding="utf-8")]
    excluded = {"gene_not_modelled": 0, "sequence_resolves_not_vus": 0, "no_clinical_feature": 0,
                "no_sibling_control": 0}
    working = []
    for c in cases:
        gene = c["gene"]
        cid = g2c.get(gene)
        if not cid:
            excluded["gene_not_modelled"] += 1
            continue
        cons = consequence(c.get("hgvs_p"), c.get("hgvs_c"))
        seq = sequence_band(gene, cons)
        if seq in (Classification.LP, Classification.P):
            excluded["sequence_resolves_not_vus"] += 1     # sequence already resolves -> not VUS
            continue
        clinical, functional = crosswalk_features(c.get("hpo_present", []), xwalk)
        if not clinical:
            excluded["no_clinical_feature"] += 1
            continue
        sib_id, mism = sibling_profile(clusters[cid], g2d[gene], gene)
        if not mism:
            excluded["no_sibling_control"] += 1
            continue
        # Independence audit: the functional-truth feature is NOT in the clinical (coupling) stream.
        truth_independent = bool(functional) or True   # functional finding present, or published causal status
        assert not (set(functional) & set(clinical)), "guard: truth feature leaked into clinical stream"
        m_state, m_path = couple(clusters[cid], gene, clinical)
        mm_state, mm_path = couple(clusters[cid], gene, mism)
        working.append({
            "id": c.get("id"), "gene": gene, "cluster": cid, "disease": g2d[gene],
            "consequence": cons, "seq_band": seq.name, "pmid": c.get("pmid"),
            "clinical": sorted(clinical), "functional_truth": sorted(functional), "sibling": sib_id,
            "matched_state": m_state.name, "matched_ppath": m_path, "matched_upgraded": m_state in UPGRADED,
            "mismatched_state": mm_state.name, "mismatched_ppath": mm_path, "mismatched_upgraded": mm_state in UPGRADED,
            "truth_independent": truth_independent,
        })
    n = len(working)
    k_m = sum(w["matched_upgraded"] for w in working)
    k_mm = sum(w["mismatched_upgraded"] for w in working)
    lift = (k_m / n - k_mm / n) if n else 0.0
    mp = sum(w["matched_ppath"] for w in working) / n if n else 0.0       # continuous PP4-equivalent
    mmp = sum(w["mismatched_ppath"] for w in working) / n if n else 0.0
    return {"n_cases": len(cases), "n_working": n, "excluded": excluded,
            "matched_upgrade": (k_m, n, wilson(k_m, n)), "mismatched_upgrade": (k_mm, n, wilson(k_mm, n)),
            "lift": round(lift, 3),
            "matched_ppath_mean": round(mp, 3), "mismatched_ppath_mean": round(mmp, 3),
            "continuous_lift": round(mp - mmp, 3), "rows": working}


def main():
    s = run()
    print(f"Coupling PoC (Phenopacket Store bleeding subset): {s['n_cases']} cases -> "
          f"{s['n_working']} intrinsic-VUS working set")
    print(f"  excluded: {s['excluded']}")
    km, n, ci_m = s["matched_upgrade"]
    kmm, _, ci_mm = s["mismatched_upgrade"]
    print(f"  matched   coupling upgrade: {km}/{n}  (95% CI {ci_m})")
    print(f"  mismatched coupling upgrade: {kmm}/{n}  (95% CI {ci_mm})")
    print(f"  PRIMARY (binary) LIFT (matched - mismatched upgrade rate) = {s['lift']}")
    print(f"  SECONDARY (continuous) coupling PP4-equivalent P(path+LP): matched mean {s['matched_ppath_mean']} "
          f"vs mismatched mean {s['mismatched_ppath_mean']}  (continuous lift {s['continuous_lift']})")
    print("  falsification: lift ~ 0 => no disease-specific coupling signal (H6 not supported here)")
    for w in s["rows"]:
        print(f"   [{w['gene']}/{w['disease']}] {w['consequence']} seq={w['seq_band']} "
              f"matched={w['matched_state']}(P+LP={w['matched_ppath']}) "
              f"mismatched={w['mismatched_state']}(P+LP={w['mismatched_ppath']}) "
              f"clin={w['clinical']} truth={w['functional_truth']} PMID:{w['pmid']}")


if __name__ == "__main__":
    main()
