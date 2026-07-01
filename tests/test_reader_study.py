"""Tests for the reader-study vignette bank + harness (coupling-proof plan, Track 4).

These are structural/validity guards (not accuracy assertions): every vignette must reference a
real cluster, a real disease in it, a valid deciding-test observation id, and a valid
contraindication token - the same source-validity discipline applied to the curated case bank.
The DISCERN-side accuracy itself is an empirical result reported in the docs, not asserted here.
"""
from __future__ import annotations

import pytest

from diseases.ontology import cluster_for, route_clusters
from eval.reader_study import load_vignettes, score_discern

VIGNETTES = load_vignettes()


def test_bank_nonempty_and_tiers():
    assert len(VIGNETTES) >= 10
    assert {v.tier for v in VIGNETTES} <= {"S", "T", "D"}
    # at least one of each meaningful tier
    assert any(v.tier == "S" for v in VIGNETTES)   # safety-divergent
    assert any(v.tier == "D" for v in VIGNETTES)   # acquired distractor


@pytest.mark.parametrize("v", VIGNETTES, ids=[v.id for v in VIGNETTES])
def test_vignette_targets_are_valid_kb_entries(v):
    if v.tier == "D":
        assert not v.gene, "distractor vignettes must be gene-less (acquired mimic)"
        return
    # validate against the vignette's DECLARED cluster (VWF routes to several clusters; the
    # runtime picks best-fit, but the source-validity check must use the intended one).
    assert v.cluster in {c.id for c in route_clusters([v.gene])}, \
        f"{v.gene} does not route to declared cluster {v.cluster}"
    cl = cluster_for(v.cluster)
    disease_ids = {d.id for d in cl.diseases}
    assert v.correct_disease_id in disease_ids, f"{v.correct_disease_id} not in {cl.id}"
    obs_ids = {o.id for o in cl.next_observations}
    if v.correct_next_obs_id is not None:
        assert v.correct_next_obs_id in obs_ids, f"bad next-obs {v.correct_next_obs_id} for {cl.id}"
    if v.harmful_tx is not None:
        contras = {c for d in cl.diseases for c in d.contraindications}
        assert v.harmful_tx in contras, f"{v.harmful_tx} is not a contraindication in {cl.id}"


def test_distractor_is_not_routed_and_safety_is_flagged():
    s = score_discern(VIGNETTES)
    # the gene-less acquired-ITP mimic must NOT be routed (correct: defer, do not over-call)
    assert s.n_skipped >= 1
    assert s.n_scored == s.n - s.n_skipped
    # every treatment-divergent safety vignette must fire a hard-stop / high-severity flag
    assert s.harmful_avoidance == 1.0
    assert 0.0 <= s.disease_accuracy <= 1.0
