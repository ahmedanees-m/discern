"""Phase 2 decision core: gap, mechanism-aware action mapping, VOI, recommend."""
from core.schemas import EvidenceContribution, PatientContext, PointsLedger, Strength, Variant
from engine.action_map import actions_for_code, expand_actions
from engine.gap import attainable_codes, detect_conflict, gap_to_target
from engine.recommend import recommend_next_action
from engine.voi import _H, lookahead, rank_actions, score_action
from rules.point_engine import Classification


def _vus_ledger(points_codes):
    v = Variant("17", 1, "C", "T", gene="ITGB3")
    contribs = [EvidenceContribution(c, s, True, "test") for c, s in points_codes]
    return PointsLedger(v, contributions=contribs)


# ---- VOI ----
def test_entropy_and_score_signs():
    assert _H(0.5) == 1.0
    assert _H(0.0) == 0.0
    s = score_action(3.0, _flow())
    assert s["eig"] >= 0 and 0.0 <= s["q"] <= 1.0
    assert s["p_pos"] > s["p"] > s["p_neg"]


def _flow():
    from core.schemas import Action
    return Action("flow_cd41_cd61", ["PS3"], 450, 7, 0.90, 0.92, "functional")


def test_rank_prefers_cheaper_higher_yield():
    from core.schemas import Action
    cheap = Action("cheap_assay", ["PS3"], 300, 5, 0.9, 0.92, "functional")
    pricey = Action("pricey_assay", ["PS3"], 3000, 30, 0.9, 0.92, "functional")
    ranked = rank_actions(3.0, [pricey, cheap])
    assert ranked[0][0].name == "cheap_assay"


# ---- mechanism-aware action mapping (the differentiator) ----
def test_mechanism_awareness_activation_vs_expression():
    expr = {a.name for a in actions_for_code("PS3", "integrin_expression")}
    activ = {a.name for a in actions_for_code("PS3", "integrin_activation")}
    assert "flow_cd41_cd61" in expr                  # Glanzmann: expression flow OK
    assert "flow_cd41_cd61" not in activ             # LAD-III: expression flow must NOT be offered
    assert "platelet_activation_pac1" in activ       # routes to activation assays instead


# ---- gap ----
def test_attainable_and_gap():
    led = _vus_ledger([("PM2", Strength.PM), ("PP3", Strength.PP)])   # +3
    assert led.points == 3.0
    assert gap_to_target(led, Classification.LP) == 3.0              # 6 - 3
    opps = {o.code for o in attainable_codes(led, PatientContext(), "integrin_expression")}
    assert "PS3" in opps                                            # assay exists for mechanism
    # PS2 needs parents -> surfaced with prerequisite
    o_ps2 = next(o for o in attainable_codes(led, PatientContext(), "integrin_expression")
                 if o.code == "PS2")
    assert "parents_available" in o_ps2.prerequisites


def test_conflict_detection():
    led = _vus_ledger([("PM3", Strength.PM), ("BS3", Strength.BS)])
    assert detect_conflict(led).in_conflict


# ---- end-to-end recommendation (design doc Glanzmann example) ----
def test_recommend_glanzmann_flow_cytometry():
    led = _vus_ledger([("PM2", Strength.PM), ("PP3", Strength.PP)])   # VUS at +3
    rec = recommend_next_action(led, PatientContext(), mechanism="integrin_expression")
    assert rec.current_class is Classification.VUS
    assert rec.ranked, "expected at least one ranked action"
    top = rec.ranked[0]
    assert top["modality"] == "functional"
    assert "PS3" in top["yields"]
    assert top["expected_post_class"] in ("LP", "P")                 # +3 +4 = +7 -> LP
    assert "Recommended" in rec.explanation


def test_lookahead_runs():
    val, seq = lookahead(2.0, expand_actions(["PS3", "PP1"], "integrin_expression"), depth=2)
    assert 0.0 <= val <= 1.0 and len(seq) <= 2
