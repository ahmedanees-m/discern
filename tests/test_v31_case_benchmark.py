"""v3.1 Track B4: curated published-case diagnosis benchmark (small-n; citations only)."""
from eval.curated_case_benchmark import load_cases, run

_FORBIDDEN = {"patient_id", "name", "dob", "mrn", "identifier", "initials"}


def test_benchmark_runs_top3_complete():
    s = run()
    assert s["n"] >= 8
    assert s["top3"] == 1.0                 # every curated case has its true dx in the top 3
    assert s["top1"] >= 0.7                 # most are Top-1; shared-feature pairs need the deciding test
    assert 0.0 <= s["abstention_rate"] <= 1.0


def test_cases_are_citation_only_no_identifiers():
    for c in load_cases():
        assert str(c.get("source_pmid", "")).isdigit(), f"{c['id']} lacks a PMID"
        assert not (_FORBIDDEN & set(c.keys())), f"{c['id']} carries a patient identifier (G7 violation)"
