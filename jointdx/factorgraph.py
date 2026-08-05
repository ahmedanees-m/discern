"""The coupled disease x variant joint model.

    P(D, V | E) proportional-to
        P(E_pheno | D) - P(E_geno | V) - P(E_func | D, V) - P(V | G, D) - P(G | D) - P(D)

Each evidence stream enters **once**: phenotype -> D, variant-intrinsic genetics -> V,
functional -> both (D,V). `P(V | G, D)` is the calibrated PP4 - the disease->variant coupling
- so phenotype's pull on V flows through D and is never double-counted. `P(G | D)` is the
gene the variant was found in, which is genetic evidence about disease identity that no other
factor owns; splitting it out is what lets a variant in F8 argue for haemophilia A over
haemophilia B. The cluster is small, so we enumerate D x V exactly. This is the answer to the
circularity critique.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from core.dx_schemas import DiscriminationCluster, Disease, Feature, VariantState
from evidence.genetic import variant_intrinsic_likelihood
from evidence.lab import func_loglik
from evidence.phenotype_lr import phenotype_loglik
from rules.vcep.loader import get_spec


@dataclass
class Evidence:
    variant_gene: str = ""
    variant_id: str = "V"
    genetic_codes: list[str] = field(default_factory=list)   # applied ACMG codes
    clinical: list[Feature] = field(default_factory=list)     # CLINICAL features (+ pertinent negatives)
    functional: list[Feature] = field(default_factory=list)   # LAB/FUNCTIONAL features


def _couple_weights(pp: float) -> dict[VariantState, float]:
    """Soft P(V | D) prior tilted by pp = P(variant pathogenic | disease) - the PP4 coupling."""
    w = {
        VariantState.PATH: pp,
        VariantState.LP: 0.7 * pp + 0.1,
        VariantState.VUS: 0.3,
        VariantState.LB: 0.7 * (1 - pp) + 0.1,
        VariantState.BEN: (1 - pp),
    }
    z = sum(w.values()) or 1.0
    return {k: v / z for k, v in w.items()}


# P(variant pathogenic | D) when the variant sits outside this disease's gene: it is an incidental
# finding, so a pathogenic call argues against the disease rather than for it. Separate from the
# gene term below, which answers a different question (where the variant is, not what it does).
INCIDENTAL_PP = 0.05


def coupling_loglik(v: VariantState, disease: Disease, variant_gene: str) -> float:
    # On-gene: the variant can be the cause of D -> tilt by p_path_given_disease.
    # Off-gene: the variant is incidental to D -> a pathogenic variant argues against D.
    pp = disease.p_path_given_disease if variant_gene in disease.genes else INCIDENTAL_PP
    return math.log(max(_couple_weights(pp)[v], 1e-9))


# P(the sequenced variant lies in this disease's gene | D), as an on-gene/off-gene pair.
#
# The value is bounded from above by an architectural constraint, not fitted to a benchmark: the
# gene term must stay weaker than the discriminating power of a cluster's deciding observation.
# The sharpest assays in the knowledge base carry likelihood ratios around 11 (RIPA mixing
# separating type 2B from platelet-type von Willebrand disease at 0.92 against 0.08). If the gene
# out-weighed those, no laboratory result could ever overturn the gene, the value-of-information
# layer would never predict a switch, and the next-test recommendation would be decorative. A
# likelihood ratio of 4 is the largest value at which the deciding assay still wins; the sweep in
# bench/phase_r_gene_term_sensitivity.py reports the whole grid, including where it breaks.
ON_GENE, OFF_GENE = 0.80, 0.20


def gene_loglik(variant_gene: str, disease: Disease) -> float:
    """P(G | D) - the disease factor's genetic term, and the reason the gene is not inert.

    `coupling_loglik` is P(V | G, D): the variant's *state* given where it sits. It normalises
    over the five variant states, so on its own the gene cancels out when the joint is marginalised
    back to a disease posterior, and a variant found in F8 says nothing about haemophilia A versus
    B. This term is the missing P(G | D) that completes the factorisation. The gene enters the
    disease factor exactly once and no other factor claims it, so the partition holds.

    It is a likelihood, not a filter, on purpose. An off-gene competitor keeps enough posterior
    mass to trip the treatment-safety interlock - a platelet-type von Willebrand case must still
    hard-stop desmopressin, because type 2B has not been excluded - and a strong enough phenotype
    can still overrule the gene rather than being silently vetoed by it.
    """
    if not variant_gene:
        return 0.0                                          # nothing sequenced -> uninformative
    return math.log(ON_GENE if variant_gene in disease.genes else OFF_GENE)


def joint(cluster: DiscriminationCluster, ev: Evidence,
          gene_evidence: bool = True) -> dict[tuple[str, VariantState], float]:
    """Normalized joint P(D, V | E) as {(disease_id, VariantState): prob}.

    `gene_evidence=False` drops the P(G | D) term. The safety interlock uses that view on purpose:
    finding a variant in one gene is not the same as ruling out a disease of another gene, because
    the competitor gene may never have been sequenced or its variant may be benign. Absence of
    evidence in gene X is not evidence that disease X is absent, so the gene is allowed to order
    the differential but never to quietly retire a treatment-divergent competitor.
    """
    spec = get_spec(ev.variant_gene) if ev.variant_gene else get_spec("")
    geno_lik = variant_intrinsic_likelihood(ev.genetic_codes, spec)   # P(E_geno | V), once
    logtbl: dict[tuple[str, VariantState], float] = {}
    for d in cluster.diseases:
        lp_ph = phenotype_loglik(ev.clinical, d)                       # P(E_pheno | D), once
        lp_gn = gene_loglik(ev.variant_gene, d) if gene_evidence else 0.0   # P(G | D), once
        lp_prior = math.log(max(d.prior, 1e-9))
        for v in VariantState:
            lp_ge = math.log(max(geno_lik[v], 1e-9))
            lp_fn = func_loglik(ev.functional, d, v)                   # P(E_func | D,V), once
            lp_cp = coupling_loglik(v, d, ev.variant_gene)            # P(V | G,D) = coupled PP4
            logtbl[(d.id, v)] = lp_ph + lp_ge + lp_fn + lp_cp + lp_gn + lp_prior
    # exp-normalize (shift by max for numerical stability)
    m = max(logtbl.values())
    exp = {k: math.exp(v - m) for k, v in logtbl.items()}
    z = sum(exp.values()) or 1.0
    return {k: v / z for k, v in exp.items()}
