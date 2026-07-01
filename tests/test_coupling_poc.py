"""Tests for the coupling proof-of-concept (eval/coupling_poc.py) and its HPO crosswalk."""
from __future__ import annotations

import yaml

from diseases.ontology import cluster_for
from eval import coupling_poc as cp
from rules.point_engine import Classification


def _all_feature_lr_keys():
    keys = set()
    for cid in cp.CLUSTERS:
        for d in cluster_for(cid).diseases:
            keys.update(d.feature_lr.keys())
    return keys


def test_crosswalk_features_are_real_cluster_features():
    xwalk = yaml.safe_load(open(cp.CROSSWALK, encoding="utf-8"))["crosswalk"]
    valid = _all_feature_lr_keys()
    for hp, e in xwalk.items():
        assert hp.startswith("HP:"), f"bad HPO id {hp}"
        assert e["kind"] in ("clinical", "functional")
        assert e["feature"] in valid, f"crosswalk feature {e['feature']} not in any cluster feature_lr"


def test_consequence_classification():
    assert cp.consequence("NP_x:p.(Gln96Ter)", "c.286C>T") == "nonsense"
    assert cp.consequence("NP_x:p.(Met2053IlefsTer31)", "c.6159_6160del") == "frameshift"
    assert cp.consequence("NP_x:p.(Cys258Arg)", "c.772T>C") == "missense"


def test_sequence_band_null_vs_missense():
    # null reaches LP/P intrinsically (PVS1); missense stays VUS intrinsically
    assert cp.sequence_band("LYST", "nonsense") in (Classification.LP, Classification.P)
    assert cp.sequence_band("LYST", "missense") == Classification.VUS


def test_coupling_distinguishes_matched_from_mismatched():
    """Mechanism check: a matched phenotype (points to the variant's-gene disease) must lift the
    variant's P(path+LP) above a mismatched one (points to a sibling disease). This validates the
    coupling + the specificity control independent of the (tiny) real working set."""
    cl = cluster_for("integrin")
    # ITGB3 = Glanzmann (gt). Matched = Glanzmann-type bleeding; mismatched = LAD-III profile.
    _, matched_p = cp.couple(cl, "ITGB3", {"glanzmann_type_bleeding": True})
    _, mismatched_p = cp.couple(cl, "ITGB3", {"leukocytosis": True, "recurrent_infections": True})
    assert matched_p > mismatched_p, (matched_p, mismatched_p)


def test_run_structure_and_independence_guard():
    s = cp.run()
    for k in ("n_cases", "n_working", "lift", "matched_ppath_mean", "mismatched_ppath_mean", "rows"):
        assert k in s
    assert s["n_working"] >= 0
    # the independence guard (functional-truth feature never in the clinical stream) holds on every row
    for w in s["rows"]:
        assert not (set(w["functional_truth"]) & set(w["clinical"]))
