"""PHASE R - sensitivity of the diagnosis arm to the strength of the P(G | D) gene term.

Phase R found that the gene the variant was sequenced in had no effect on the disease posterior:
P(V | G, D) normalises over the five variant states, so the gene cancelled out on marginalisation
and a variant in F8 argued no more for haemophilia A than for haemophilia B. Adding the missing
P(G | D) fixes that, but its strength is a modelling choice, so this module reports what that
choice actually buys and costs instead of asserting a value.

Three quantities are swept together, because they trade off against each other:
  diagnosis    - curated Top-1, pooled and split by whether the gene already determines the answer
  safety       - hard-stop sensitivity and specificity on the treatment-divergent scenarios
  decidability - whether a cluster's deciding observation can still overturn the gene. If the gene
                 out-weighs the sharpest assay in the knowledge base, no laboratory result can ever
                 change the answer and the value-of-information layer stops meaning anything.

The committed value is the largest one that leaves decidability intact, which is a constraint from
the architecture rather than from the benchmark.

Run:  python -m bench.phase_r_gene_term_sensitivity
"""
from __future__ import annotations

import json
import os

from core.dx_schemas import Feature, FeatureKind
from diseases.ontology import cluster_for
from eval.gene_only_baseline import _stratum, gene_lookup_rank
from jointdx import factorgraph
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import marginal_disease
from nextobs.whatif import whatif

HERE = os.path.dirname(__file__)
OUT_JSON = os.path.join(HERE, "phase_r_gene_term_sensitivity.json")

# on-gene / off-gene pairs, from "the gene is inert" upward
GRID = [(0.50, 0.50), (0.60, 0.40), (0.70, 0.30), (0.80, 0.20), (0.90, 0.10), (0.95, 0.05)]


def _decidable():
    """Can the deciding assay still overturn the gene, on the path the clinician actually sees?

    This is the what-if the recommender shows before a test is ordered: for a VWF case with
    enhanced low-dose RIPA, a mixing study that comes back platelet-origin must point at
    platelet-type von Willebrand disease, whose gene is GP1BA, and plasma-origin at type 2B.
    RIPA mixing is the sharpest discriminator in the knowledge base, so if the gene out-weighs it
    here it out-weighs every assay everywhere, and the value-of-information layer can no longer
    predict a switch. Checked through whatif() rather than a hand-built posterior, because that is
    the code path the recommendation is drawn from.
    """
    cluster = cluster_for("vwf_gpib")
    ev = Evidence(variant_gene="VWF", clinical=[
        Feature("ripa_low_dose_enhanced", FeatureKind.LAB, True, observed=True)])
    shifts = whatif(cluster, ev, "ripa_mixing")
    if not shifts:
        return False
    return (shifts.get("platelet_origin", (None,))[0] == "ptvwd"
            and shifts.get("plasma_origin", (None,))[0] == "vwd2b")


def _safety():
    from bench.track3_trustworthiness import safety_interlock
    s = safety_interlock()
    return s["hardstop_sensitivity"], s["hardstop_specificity"]


def _diagnosis():
    from eval.curated_case_benchmark import _evidence, load_cases
    cases = load_cases()
    by_stratum: dict[str, list[int]] = {}
    hits = []
    for c in cases:
        cl = cluster_for(c["cluster"])
        md = marginal_disease(joint(cl, _evidence(c)))
        ok = int(max(md, key=md.get) == c["true_dx"])
        hits.append(ok)
        by_stratum.setdefault(_stratum(c), []).append(ok)
    return (round(sum(hits) / len(hits), 4),
            {k: round(sum(v) / len(v), 4) for k, v in sorted(by_stratum.items())})


def run() -> dict:
    from eval.curated_case_benchmark import load_cases
    cases = load_cases()
    gene_lookup_top1 = round(
        sum(gene_lookup_rank(c)[0] == c["true_dx"] for c in cases) / len(cases), 4)

    on0, off0 = factorgraph.ON_GENE, factorgraph.OFF_GENE
    rows = []
    try:
        for on, off in GRID:
            factorgraph.ON_GENE, factorgraph.OFF_GENE = on, off
            top1, by_stratum = _diagnosis()
            sens, spec = _safety()
            rows.append({"on_gene": on, "off_gene": off, "likelihood_ratio": round(on / off, 2),
                         "curated_top1": top1, "curated_top1_by_stratum": by_stratum,
                         "hardstop_sensitivity": sens, "hardstop_specificity": spec,
                         "deciding_assay_can_overturn_gene": _decidable(),
                         "committed": (on, off) == (on0, off0)})
    finally:
        factorgraph.ON_GENE, factorgraph.OFF_GENE = on0, off0

    return {
        "gene_lookup_top1_baseline": gene_lookup_top1,
        "committed_value": {"on_gene": on0, "off_gene": off0, "likelihood_ratio": round(on0 / off0, 2)},
        "selection_rule": ("the largest likelihood ratio that still lets a cluster's deciding "
                           "observation overturn the gene; the sharpest assay in the knowledge "
                           "base carries a ratio near 11, so the gene term is held well below it"),
        "grid": rows,
    }


def main():
    o = run()
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(o, fh, indent=2)
    print("== PHASE R: how much does the gene term matter, and what does it cost? ==")
    print(f"gene-lookup baseline Top-1 = {o['gene_lookup_top1_baseline']:.0%}   "
          f"committed LR = {o['committed_value']['likelihood_ratio']}")
    print(f"\n{'LR':>6}  {'Top-1':>6}  {'shared':>7}  {'unique':>7}  {'sens':>5}  {'spec':>5}  "
          f"{'assay wins':>10}")
    for r in o["grid"]:
        st = r["curated_top1_by_stratum"]
        mark = " <- committed" if r["committed"] else ""
        print(f"{r['likelihood_ratio']:>6}  {r['curated_top1']:>6.0%}  "
              f"{st.get('shared_gene', 0):>7.0%}  {st.get('unique_gene', 0):>7.0%}  "
              f"{r['hardstop_sensitivity']:>5.2f}  {r['hardstop_specificity']:>5.2f}  "
              f"{str(r['deciding_assay_can_overturn_gene']):>10}{mark}")
    print(f"\nrule: {o['selection_rule']}")
    print(f"wrote {os.path.relpath(OUT_JSON)}")


if __name__ == "__main__":
    main()
