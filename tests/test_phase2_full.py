"""Decision core completeness: risk gate, Pareto, conflict, recessive, audit."""
from core.audit import AuditLedger
from core.schemas import (
    Action,
    EvidenceContribution,
    PatientContext,
    PointsLedger,
    Strength,
    Variant,
)
from engine.recommend import recommend_next_action
from engine.voi import outcome_variance, pareto_frontier, rank_actions, risk_gate
from rules.point_engine import Classification


def _ledger(codes, gene="ITGB3"):
    v = Variant("17", 1, "C", "T", gene=gene)
    return PointsLedger(v, contributions=[EvidenceContribution(c, s, True, "t") for c, s in codes])


def _act(name, cost, days, sens=0.9, spec=0.92, mod="functional", yields=("PS3",)):
    return Action(name, list(yields), cost, days, sens, spec, mod)


# ---- Section 3.5 Pareto frontier ----
def test_pareto_frontier_is_nondominated():
    cheap = _act("cheap", 300, 5)                 # sens .9 / spec .92
    dup_costlier = _act("dup_costlier", 900, 20)  # identical accuracy, higher cost -> dominated
    front = pareto_frontier(3.0, [cheap, dup_costlier])
    names = [f["action"] for f in front]
    assert names[0] == "cheap"
    assert "dup_costlier" not in names            # equal EIG + higher cost -> excluded
    costs = [f["cost_usd"] for f in front]
    eigs = [f["eig_bits"] for f in front]
    assert costs == sorted(costs)                                          # cost-ordered
    assert all(b > a for a, b in zip(eigs, eigs[1:], strict=False))        # strictly increasing EIG


# ---- Section 3.6 risk gate + variance ----
def test_risk_gate_drops_low_utility():
    weak = _act("weak", 300, 5, sens=0.5, spec=0.5)
    strong = _act("strong", 300, 5, sens=0.95, spec=0.95)
    ranked = rank_actions(3.0, [weak, strong])
    gated = risk_gate(ranked, tau=0.05)
    kept = {t[0].name for t in gated}
    assert "strong" in kept
    assert all(t[2]["delta_utility"] >= 0.05 for t in gated)


def test_outcome_variance_nonnegative():
    assert outcome_variance(3.0, _act("a", 300, 5)) >= 0.0


# ---- Section 6 conflict -> orthogonal disambiguating action ----
def test_conflict_prefers_orthogonal():
    led = _ledger([("PM3", Strength.PM), ("BP4", Strength.BP)])   # opposing -> conflict
    rec = recommend_next_action(led, PatientContext(), mechanism="integrin_expression")
    assert rec.in_conflict
    assert "conflict" in rec.explanation.lower() or "VUS-by-conflict" in rec.explanation
    # all recommended actions provide orthogonal (independent) evidence
    assert all(r["modality"] in {"functional", "rna", "segregation"} for r in rec.ranked)


# ---- Section 3.1/Section 6 recessive case note ----
def test_recessive_case_note():
    led = _ledger([("PM2", Strength.PM), ("PP3", Strength.PP)])
    rec = recommend_next_action(led, PatientContext(), mechanism="integrin_expression",
                                inheritance="AR")
    assert "Recessive" in rec.case_note
    assert "Recessive" in rec.explanation


# ---- audit integration ----
def test_recommendation_audited():
    led = _ledger([("PM2", Strength.PM), ("PP3", Strength.PP)])
    audit = AuditLedger()
    rec = recommend_next_action(led, PatientContext(), mechanism="integrin_expression",
                                audit=audit)
    assert rec.audit_index == 0
    assert audit.verify()
    assert audit.entries[0].kind == "recommend"


# ---- risk gate inside recommend (high tau -> may empty the list, fallback) ----
def test_recommend_runs_with_target():
    led = _ledger([("PP3", Strength.PP)])
    rec = recommend_next_action(led, PatientContext(), mechanism="integrin_expression",
                                target=Classification.LP)
    assert rec.current_class is Classification.VUS
    assert rec.pareto_frontier         # frontier computed
