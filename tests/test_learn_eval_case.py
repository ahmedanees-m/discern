"""Stats (Phase 6.1), learning loop (Phase 5.2), and whole-odyssey policy (Phase 4)."""
from engine.case_policy import CaseState, decide_case
from eval.stats import bootstrap_ci, mcnemar, wilcoxon_signed_rank
from learn.outcome_store import Outcome, OutcomeStore
from learn.prior_update import BetaPrior, update_sensitivity


# ---- stats ----
def test_mcnemar_directionality():
    r = mcnemar(b=30, c=5)             # A resolves many that B does not
    assert r.p_value < 0.05
    assert mcnemar(b=10, c=10).p_value > 0.5     # symmetric -> non-significant


def test_wilcoxon_detects_cost_reduction():
    full = [100, 120, 90, 110, 95]
    tools_only = [150, 170, 140, 160, 145]       # consistently higher cost
    r = wilcoxon_signed_rank(full, tools_only)
    assert r.p_value < 0.1


def test_bootstrap_ci_brackets_mean():
    mean, lo, hi = bootstrap_ci([1, 2, 3, 4, 5], seed=1)
    assert lo <= mean <= hi and abs(mean - 3.0) < 1e-9


# ---- learning loop ----
def test_beta_update_is_attributable():
    p = BetaPrior()
    assert p.mean == 0.5
    p2 = p.update(True).update(True).update(False)
    assert update_sensitivity(1, 1, True) == (2, 1)
    assert p2.mean > 0.5 and p2.n == 3


def test_outcome_store_prior():
    s = OutcomeStore()
    for resolved in (True, True, True, False):
        s.record(Outcome("F8", "factor_activity", True, resolved, accepted=True))
    pr = s.assay_prior("F8", "factor_activity")
    assert pr.mean == 4 / 6          # alpha=1+3, beta=1+1 -> 4/6
    assert s.acceptance_rate() == 1.0
    assert OutcomeStore.from_json(s.to_json()).acceptance_rate() == 1.0


# ---- whole-odyssey policy ----
def test_case_policy_routes():
    assert decide_case(CaseState(has_vus_candidate=True, has_any_candidate=True)).route == "variant_voi"
    splice = decide_case(CaseState(False, True, suspected_splicing=True))
    assert splice.route == "modality_escalation" and splice.modality.name == "rna_seq"
    struct = decide_case(CaseState(False, True, suspected_structural=True))
    assert struct.modality.name in ("long_read_genome", "optical_genome_mapping")
    nores = decide_case(CaseState(False, False))
    assert nores.route in ("modality_escalation", "reanalysis")
    novel = decide_case(CaseState(False, True, novel_gene_candidate=True))
    assert novel.route == "matchmaking"
