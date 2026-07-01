"""DISCERN Phase 0/1.1: extended schemas, VCEP loader, per-code partition, genetic factor."""
from core.dx_schemas import (
    DiscriminationCluster,
    Disease,
    Feature,
    FeatureKind,
    Observation,
    VariantState,
)
from evidence.genetic import variant_intrinsic_likelihood, variant_intrinsic_points
from rules.vcep.loader import get_spec
from rules.vcep.partition import is_variant_intrinsic, owner, partition


# ---- per-code partition (the circularity fix, Gate G3) ----
def test_partition_routes_each_code_to_one_factor():
    assert is_variant_intrinsic("PM2") and is_variant_intrinsic("PVS1") and is_variant_intrinsic("BP4")
    assert owner("PP4") == "disease_pp4"          # phenotype -> coupling, NOT genetic
    assert owner("PS3") == "functional"
    assert owner("BS3") == "functional"
    assert owner("PP1") == "segregation"
    assert owner("PM3") == "phasing"
    assert owner("PS2") == "denovo"
    # strength suffix is stripped
    assert is_variant_intrinsic("PM2_Supporting")


def test_partition_groups():
    g = partition(["PM2", "PP4", "PS3", "PP1", "PM3"])
    assert g["variant_intrinsic"] == ["PM2"]
    assert g["disease_pp4"] == ["PP4"] and g["functional"] == ["PS3"]


# ---- VCEP loader ----
def test_spec_loads_and_reduced_confidence():
    gt = get_spec("ITGB3")
    assert gt.covered and gt.disease == "Glanzmann thrombasthenia" and gt.inheritance == "AR"
    f8 = get_spec("F8")
    assert f8.covered and f8.inheritance == "XLR"
    runx1 = get_spec("RUNX1")
    assert runx1.covered and "RUNX1" in runx1.genes
    lad3 = get_spec("FERMT3")          # no VCEP spec -> reduced confidence (Gate G2)
    assert not lad3.covered


# ---- genetic factor uses ONLY variant-intrinsic codes ----
def test_genetic_factor_excludes_nonintrinsic_codes():
    spec = get_spec("ITGB3")
    full = ["PVS1", "PM2", "PP4", "PS3", "PP1", "PM3"]   # PP4/PS3/PP1/PM3 owned elsewhere
    intrinsic_only = ["PVS1", "PM2"]
    # the genetic factor must give the SAME points whether or not the other codes are present
    assert variant_intrinsic_points(full, spec) == variant_intrinsic_points(intrinsic_only, spec)


def test_genetic_likelihood_peaks_correctly():
    spec = get_spec("ITGB3")
    hi = variant_intrinsic_likelihood(["PVS1", "PS1", "PM2"], spec)   # strong pathogenic points
    assert abs(sum(hi.values()) - 1.0) < 1e-9
    assert max(hi, key=hi.get) in (VariantState.PATH, VariantState.LP)
    lo = variant_intrinsic_likelihood(["BA1", "BP4"], spec)            # benign points
    assert max(lo, key=lo.get) in (VariantState.BEN, VariantState.LB)


# ---- schemas round-trip ----
def test_dx_schemas_construct():
    d = Disease("gt", "Glanzmann", ["ITGA2B", "ITGB3"], "AR", "integrin_expression",
                prior=0.4, p_path_given_disease=0.9, treatment="antifibrinolytics",
                contraindications=[], vcep_spec="PD-VCEP-ITGA2B_ITGB3")
    o = Observation("flow", "flow CD41/CD61", "functional", informs=["D", "V"],
                    outcome_lr={"absent": {"gt": 20.0}}, changes_management=False, accessibility="high")
    c = DiscriminationCluster("integrin", "Integrin/aggregation", [d],
                              discriminating_features=["leukocytosis"], next_observations=[o])
    f = Feature("leukocytosis", FeatureKind.CLINICAL, True, observed=False)   # pertinent negative
    assert c.diseases[0].genes == ["ITGA2B", "ITGB3"] and not f.observed
