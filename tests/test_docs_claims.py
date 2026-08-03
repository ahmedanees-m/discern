"""Guards on the documentation's claims, not its prose.

Phase R retired two claims, withdrew a third, and corrected one factual statement about the
partition. Each of those survived at least one correction pass somewhere in the archive, which is
why they are asserted here: a grep is a cheap way to stop a superseded sentence being reintroduced
by a copy-paste from an older document.
"""
from __future__ import annotations

import os
import pathlib

import pytest

DOCS = pathlib.Path(__file__).resolve().parents[1] / "docs"
README = pathlib.Path(__file__).resolve().parents[1] / "README.md"

# A superseded claim may still be quoted - that is how a retraction is written. What must never
# happen is the claim appearing as an assertion. So the rule is per line: the phrase is allowed
# only when the same line also marks it as retired, withdrawn or superseded.
SUPERSEDED = {
    "owned by the disease and coupling factors": (
        "PM1, PM5, PS1 and PS4 are variant-intrinsic; only PS3 and PP4 are routed away. This "
        "phrasing was the factual error Phase R corrected."),
    "only tool emitting a calibrated probability": (
        "retired by Phase R R2: under an identical isotonic protocol the comparators calibrate too, "
        "with overlapping intervals."),
    "monotone risk-coverage": "withdrawn by Phase R: the diagnosis arm is at ceiling at n=42.",
    "monotone payoff": "withdrawn by Phase R: the diagnosis arm is at ceiling at n=42.",
}
RETRACTION_MARKERS = ("retire", "Retire", "RETIRE", "withdraw", "Withdraw", "WITHDRAW",
                      "superseded", "Superseded", "SUPERSEDED", "Phase R", "was:", "no longer",
                      "not claim", "corrected", "uninformative")


def _markdown_files():
    return [f for f in [README, *sorted(DOCS.glob("*.md"))] if f.exists()]


def _offending_lines(text):
    """Flag a superseded phrase only when nothing nearby retracts it.

    Markdown hard-wraps prose, so a claim and its retraction routinely land on adjacent lines. The
    check therefore looks at a small window around the hit rather than the single line, which keeps
    the guard honest without forcing documents to be rewrapped around a test.
    """
    lines = text.splitlines()
    out = []
    for i, line in enumerate(lines):
        window = " ".join(lines[max(0, i - 1):i + 2])
        for phrase, why in SUPERSEDED.items():
            if phrase in line and not any(m in window for m in RETRACTION_MARKERS):
                out.append((i + 1, phrase, why, line.strip()[:160]))
    return out


@pytest.mark.parametrize("path", _markdown_files(), ids=lambda p: p.name)
def test_superseded_claims_are_never_asserted(path):
    """A retired claim may be quoted while being retired, never stated as fact."""
    bad = _offending_lines(path.read_text(encoding="utf-8"))
    assert not bad, "\n".join(
        f"{path.name}:{ln} asserts '{ph}' with no retraction marker - {why}\n    {snippet}"
        for ln, ph, why, snippet in bad)


def test_manuscript_when_available_asserts_no_superseded_claim():
    """The manuscript lives outside the repository; check it when the path is supplied."""
    ms = os.environ.get("DISCERN_MANUSCRIPT_PATH")
    if not ms or not os.path.exists(ms):
        pytest.skip("manuscript not available in this checkout (kept outside the public repo)")
    bad = _offending_lines(pathlib.Path(ms).read_text(encoding="utf-8"))
    assert not bad, "\n".join(
        f"manuscript:{ln} asserts '{ph}'\n    {snippet}" for ln, ph, _, snippet in bad)
