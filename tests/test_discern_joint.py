"""Coupled joint model and the circularity guarantee."""
from core.dx_schemas import Feature, FeatureKind, VariantState
from diseases.ontology import cluster_for
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import leading_disease, marginal_variant, reclassify


def _clin(fid, present):
    return Feature(fid, FeatureKind.CLINICAL, present, observed=present)


def _func(fid, value):
    return Feature(fid, FeatureKind.FUNCTIONAL, value, observed=True)


# ---- flagship discrimination: GT vs LAD-III ----
def test_gt_vs_lad3_no_leukocytosis_favors_gt():
    cluster = cluster_for("integrin")
    ev = Evidence(variant_gene="ITGB3", genetic_codes=["PVS1", "PM2"],
                  clinical=[_clin("glanzmann_type_bleeding", True),
                            _clin("leukocytosis", False),         # pertinent negative
                            _clin("recurrent_infections", False)])
    lead, p = leading_disease(joint(cluster, ev))
    assert lead == "gt" and p > 0.5          # absence of leukocytosis argues against LAD-III


def test_lad3_with_leukocytosis_and_infections():
    cluster = cluster_for("integrin")
    ev = Evidence(variant_gene="FERMT3", genetic_codes=["PVS1"],
                  clinical=[_clin("glanzmann_type_bleeding", True),
                            _clin("leukocytosis", True),
                            _clin("recurrent_infections", True)])
    lead, p = leading_disease(joint(cluster, ev))
    assert lead == "lad3"                     # leukocytosis + infections + FERMT3 variant


# ---- VUS reclassification via the disease model + functional ----
def test_functional_upgrades_vus():
    cluster = cluster_for("integrin")
    # weak genetic (PM2 only -> VUS-ish), but glanzmann bleeding + ABSENT alphaIIbbeta3 (functional)
    ev = Evidence(variant_gene="ITGB3", genetic_codes=["PM2"],
                  clinical=[_clin("glanzmann_type_bleeding", True), _clin("leukocytosis", False)],
                  functional=[_func("aiib3_expression_absent", "absent")])
    old, new, drivers = reclassify(joint(cluster, ev), old_state=VariantState.VUS)
    assert new in (VariantState.LP, VariantState.PATH)      # disease+functional move it up
    assert old is VariantState.VUS


# ---- the per-code circularity guarantee ----
def test_each_code_enters_once():
    cluster = cluster_for("integrin")
    base = Evidence(variant_gene="ITGB3", genetic_codes=["PVS1", "PM2"],
                    clinical=[_clin("glanzmann_type_bleeding", True)])
    # add PP4/PS3/PP1/PM3 as *codes* - they are owned by other factors, so the joint must
    # be IDENTICAL (no double-counting through the genetic factor).
    plus = Evidence(variant_gene="ITGB3",
                    genetic_codes=["PVS1", "PM2", "PP4", "PS3", "PP1", "PM3"],
                    clinical=[_clin("glanzmann_type_bleeding", True)])
    mv_base = marginal_variant(joint(cluster, base))
    mv_plus = marginal_variant(joint(cluster, plus))
    for s in VariantState:
        assert abs(mv_base[s] - mv_plus[s]) < 1e-12          # no inflation from re-added codes
