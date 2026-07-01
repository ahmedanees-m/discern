"""Gate G3 - VCEP reconstruction with no inflation (plan Phase 9.2 / Section A.1).

Feeding the joint model the SAME variant-intrinsic evidence the VCEP had must reconstruct
the VCEP's variant classification, with no inflation from the bundled codes (PP4/PS3/PP1/
PM3) that DISCERN owns in other factors. Passing this is the precondition for trusting any
reclassification claim - it proves the per-code partition prevents double-counting.
"""
from core.dx_schemas import Feature, FeatureKind, VariantState
from diseases.ontology import cluster_for
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import marginal_variant, variant_call


def _ev(codes, gene="ITGB3"):
    return Evidence(variant_gene=gene, genetic_codes=codes,
                    clinical=[Feature("glanzmann_type_bleeding", FeatureKind.CLINICAL, True, observed=True)])


def test_pathogenic_codes_reconstruct_pathogenic():
    # PVS1(8)+PS1(4)+PM2(1 in ITGB3 spec) = 13 -> Pathogenic band
    v, p = variant_call(joint(cluster_for("integrin"), _ev(["PVS1", "PS1", "PM2"])))
    assert v in (VariantState.PATH, VariantState.LP)


def test_benign_codes_reconstruct_benign():
    v, p = variant_call(joint(cluster_for("integrin"), _ev(["BA1", "BP4"])))
    assert v in (VariantState.BEN, VariantState.LB)


def test_no_inflation_from_bundled_codes():
    cluster = cluster_for("integrin")
    intrinsic = _ev(["PVS1", "PS1", "PM2"])
    # the VCEP's bottom-line label would also include PP4/PS3/PP1/PM3 - adding them as codes
    # must NOT inflate the variant marginal (they are owned by other factors).
    bundled = _ev(["PVS1", "PS1", "PM2", "PP4", "PS3", "PP1", "PM3"])
    a = marginal_variant(joint(cluster, intrinsic))
    b = marginal_variant(joint(cluster, bundled))
    for s in VariantState:
        assert abs(a[s] - b[s]) < 1e-12
