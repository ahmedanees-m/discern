"""v3.1 Track B1/B3: new + enhanced discrimination clusters and their per-cluster safety map.

Includes a CI guard for the non-negotiable provenance rule (every feature_lr value carries a
PMID), plus discrimination + contraindication checks for C4/C5/C7/C10.
"""
from core.dx_schemas import Feature, FeatureKind
from diseases.ontology import all_clusters, cluster_for
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import leading_disease, marginal_disease
from safety.interlock import flags


def _f(fid, kind=FeatureKind.LAB, present=True):
    return Feature(fid, kind, present, observed=present)


def test_every_feature_lr_is_sourced():
    """Non-negotiable: no fabricated LRs - every feature_lr value = [freq, n, PMID-digits]."""
    for cid, c in all_clusters().items():
        for d in c.diseases:
            for feat, lr in d.feature_lr.items():
                assert len(lr) == 3, f"{cid}/{d.id}/{feat}: expected [freq,n,pmid]"
                freq, n, pmid = lr
                assert 0.0 <= float(freq) <= 1.0, f"{cid}/{d.id}/{feat}: freq out of range"
                assert str(pmid).isdigit() and len(str(pmid)) >= 6, \
                    f"{cid}/{d.id}/{feat}: PMID '{pmid}' not a real identifier"


def test_all_ten_clusters_present():
    cl = all_clusters()
    for cid in ("integrin", "vwf_gpib", "macrothrombocytopenia", "granule", "thr_leukemia",
                "coag_factor", "vwd2n_hema", "alpha_granule", "mild_vwd", "scott"):
        assert cid in cl, f"missing cluster {cid}"


def test_c4_splenectomy_hard_stop_for_inherited_thrombocytopenia():
    """C4: planned splenectomy must hard-stop when an inherited thrombocytopenia is plausible
    (avoid the ITP-misdiagnosis trap: needless splenectomy + missed surveillance)."""
    c = cluster_for("thr_leukemia")
    ev = Evidence(variant_gene="RUNX1", genetic_codes=["PM2"],
                  clinical=[_f("family_history_malignancy", FeatureKind.CLINICAL),
                            _f("lifelong_thrombocytopenia", FeatureKind.CLINICAL),
                            _f("ad_inheritance_pattern", FeatureKind.CLINICAL)])
    fl = flags(c, joint(c, ev), planned_tx="splenectomy")
    assert any("HARD STOP" in f.message and "splenectomy" in f.message for f in fl)


def test_c5_vwf_fviii_binding_points_to_2n():
    c = cluster_for("vwd2n_hema")
    ev = Evidence(variant_gene="VWF", clinical=[_f("reduced_vwf_fviii_binding"),
                                                _f("low_fviii_isolated")])
    md = marginal_disease(joint(c, ev))
    assert md["vwd2n"] > md["mild_hema"]            # the binding defect separates 2N from hem A


def test_c7_quebec_contraindicates_platelet_transfusion():
    c = cluster_for("alpha_granule")
    ev = Evidence(variant_gene="PLAU", clinical=[_f("urokinase_overexpression")])
    j = joint(c, ev)
    assert leading_disease(j)[0] == "quebec"
    fl = flags(c, j, planned_tx="platelet_transfusion")
    assert any("HARD STOP" in f.message and f.competitor_id == "quebec" for f in fl)


def test_c10_ps_exposure_detects_scott():
    c = cluster_for("scott")
    ev = Evidence(variant_gene="ANO6", clinical=[_f("reduced_ps_exposure"),
                                                 _f("abnormal_prothrombinase_assay")])
    md = marginal_disease(joint(c, ev))
    assert md["scott"] > md["normal_workup"]        # the PS-exposure assay is decisive
