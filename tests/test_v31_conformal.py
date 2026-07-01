"""v3.1 Track B2: Mondrian split-conformal selective prediction.

Synthetic-sampling SANITY check of the coverage guarantee (not a clinical validation - real
coverage is cohort-gated, H5/G11). Cases are drawn from the cluster's own feature frequencies,
so conformal coverage must hold at >= 1 - alpha (minus finite-sample slack).
"""
import random

from core.dx_schemas import Feature, FeatureKind
from diseases.ontology import cluster_for
from jointdx.conformal import calibrate, empirical_coverage, prediction_set, selective_predict
from jointdx.factorgraph import Evidence


def _sample(cluster, rng):
    """Draw a disease by prior, then sample its features present w.p. freq (model-faithful)."""
    tot = sum(d.prior for d in cluster.diseases)
    r, acc = rng.random() * tot, 0.0
    chosen = cluster.diseases[-1]
    for d in cluster.diseases:
        acc += d.prior
        if r <= acc:
            chosen = d
            break
    clinical = [Feature(f, FeatureKind.LAB, rng.random() < float(lr[0]), observed=True)
                for f, lr in chosen.feature_lr.items()]
    gene = chosen.genes[0] if chosen.genes else ""
    return Evidence(variant_gene=gene, clinical=clinical), chosen.id


def _dataset(cluster_id, n, seed):
    rng = random.Random(seed)
    c = cluster_for(cluster_id)
    return c, [_sample(c, rng) for _ in range(n)]


def test_conformal_coverage_guarantee_holds():
    c, cal = _dataset("integrin", 500, seed=1)
    _, test = _dataset("integrin", 500, seed=2)
    alpha = 0.1
    q = calibrate(c, cal, alpha=alpha)
    cov = empirical_coverage(c, test, q)
    # marginal coverage must be at/above 1-alpha minus finite-sample slack
    assert cov["overall"] >= (1 - alpha) - 0.05, cov


def test_selective_commit_and_abstain():
    c = cluster_for("integrin")
    _, cal = _dataset("integrin", 400, seed=3)
    q = calibrate(c, cal, alpha=0.1)
    # strong, coherent GT evidence -> should commit to a singleton {gt}
    strong = Evidence(variant_gene="ITGB3", genetic_codes=["PVS1", "PS1", "PM2"],
                      clinical=[Feature("glanzmann_type_bleeding", FeatureKind.CLINICAL, True, observed=True),
                                Feature("aiib3_expression_absent", FeatureKind.LAB, True, observed=True),
                                Feature("leukocytosis", FeatureKind.LAB, False, observed=True)])
    r = selective_predict(c, strong, q)
    assert isinstance(r.pred_set, set) and r.leading == "gt"
    # ambiguous (no evidence) -> conformal set should not be a confident singleton -> abstain
    vague = selective_predict(c, Evidence(variant_gene=""), q)
    assert vague.committed is False


def test_prediction_set_contains_leading_under_strong_evidence():
    c = cluster_for("integrin")
    _, cal = _dataset("integrin", 400, seed=4)
    q = calibrate(c, cal, alpha=0.1)
    ev = Evidence(variant_gene="ITGB3", genetic_codes=["PVS1", "PM2"],
                  clinical=[Feature("glanzmann_type_bleeding", FeatureKind.CLINICAL, True, observed=True)])
    assert "gt" in prediction_set(c, ev, q)
