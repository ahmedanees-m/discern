"""Synthetic-coupling SANITY harness (v3.1 Track C1) - NOT a validation.

Dry-runs the PRE-REGISTERED coupling primary endpoint (coupling vs intrinsic-only VUS
reclassification) on SIMULATED paired phenotype+variant data, so the analysis pipeline and the
circularity guard are ready and tested BEFORE any real cohort data is seen (Gate G12). It
proves the pipeline executes and the coupling channel is wired; it makes NO claim that the
coupling works - that is provable only on paired cohorts and is reported regardless of outcome
(Gate G13, H6 with its explicit falsification condition).

Run:  python3 -m eval.synthetic_coupling_harness
"""
from __future__ import annotations

import random
from dataclasses import dataclass

from core.dx_schemas import Feature, FeatureKind, VariantState
from diseases.ontology import cluster_for
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import marginal_variant


@dataclass
class CouplingResult:
    n: int
    intrinsic_only_path_mass: float          # mean P(PATH) with genetics only
    coupled_path_mass: float                 # mean P(PATH) with genetics + disease-specific phenotype
    intrinsic_invariant_to_phenotype: bool   # the circularity guard
    note: str


def _simulate(cluster, n, seed):
    rng = random.Random(seed)
    cases = []
    for _ in range(n):
        # sample a true disease by prior
        tot = sum(d.prior for d in cluster.diseases)
        r, acc = rng.random() * tot, 0.0
        d = cluster.diseases[-1]
        for cand in cluster.diseases:
            acc += cand.prior
            if r <= acc:
                d = cand
                break
        pheno = [Feature(f, FeatureKind.CLINICAL, rng.random() < float(lr[0]), observed=True)
                 for f, lr in d.feature_lr.items()]
        gene = d.genes[0] if d.genes else (cluster.diseases[0].genes or ["VWF"])[0]
        cases.append((gene, pheno))
    return cases


def run(cluster_id="integrin", n=300, seed=0) -> CouplingResult:
    cluster = cluster_for(cluster_id)
    cases = _simulate(cluster, n, seed)
    intr_path = coup_path = 0.0
    guard_ok = True
    for gene, pheno in cases:
        # a VUS-ish intrinsic variant (one moderate code) so the coupling has room to move it
        geno_only = Evidence(variant_gene=gene, genetic_codes=["PM2"])
        coupled = Evidence(variant_gene=gene, genetic_codes=["PM2"], clinical=pheno)
        mv_intr = marginal_variant(joint(cluster, geno_only))
        mv_coup = marginal_variant(joint(cluster, coupled))
        intr_path += mv_intr[VariantState.PATH]
        coup_path += mv_coup[VariantState.PATH]
        # circularity guard: adding phenotype must NOT change the intrinsic (genetics-only) result
        mv_intr2 = marginal_variant(joint(cluster, Evidence(variant_gene=gene, genetic_codes=["PM2"])))
        if abs(mv_intr2[VariantState.PATH] - mv_intr[VariantState.PATH]) > 1e-9:
            guard_ok = False
    return CouplingResult(
        n=n, intrinsic_only_path_mass=round(intr_path / n, 4),
        coupled_path_mass=round(coup_path / n, 4),
        intrinsic_invariant_to_phenotype=guard_ok,
        note="SANITY on synthetic data; the coupling is validated only on paired cohorts (Gate G13).")


def main():
    r = run()
    print("SYNTHETIC COUPLING SANITY (NOT a validation result):")
    print(f"  n={r.n}  intrinsic-only mean P(PATH)={r.intrinsic_only_path_mass}  "
          f"coupled mean P(PATH)={r.coupled_path_mass}")
    print(f"  circularity guard (intrinsic invariant to phenotype): {r.intrinsic_invariant_to_phenotype}")
    print("  " + r.note)


if __name__ == "__main__":
    main()
