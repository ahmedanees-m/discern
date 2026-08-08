"""Guards on the documentation's claims, not its prose.

Three statements are not supported by the results and one is factually wrong about the evidence
partition. Each is easy to reintroduce by copying a sentence from an older draft, so each is
asserted against here: the phrase may appear while being disclaimed, never as an assertion.
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
        "PM1, PM5, PS1 and PS4 are variant-intrinsic; only PS3 and PP4 are routed away."),
    "only tool emitting a calibrated probability": (
        "under an identical isotonic protocol the comparators calibrate too, with overlapping "
        "intervals."),
    "monotone risk-coverage": "the diagnosis arm is at ceiling at n=42, so none is demonstrable.",
    "monotone payoff": "the diagnosis arm is at ceiling at n=42, so none is demonstrable.",
}
# Case-insensitive: a disclaimer is a disclaimer wherever it falls in a sentence.
RETRACTION_MARKERS = ("retire", "withdraw", "superseded", "was:", "no longer", "not claim",
                      "is not claimed", "no monotone", "corrected", "uninformative")


def _markdown_files():
    return [f for f in [README, *sorted(DOCS.glob("*.md"))] if f.exists()]


def _offending_lines(text):
    """Flag a superseded phrase only when nothing nearby retracts it.

    Markdown hard-wraps prose, so a claim and its retraction routinely land on adjacent lines. The
    check therefore looks at a small window around the hit rather than the single line, so
    a document need not be rewrapped to satisfy the guard.
    """
    lines = text.splitlines()
    out = []
    for i, line in enumerate(lines):
        window = " ".join(lines[max(0, i - 1):i + 2])
        for phrase, why in SUPERSEDED.items():
            low = window.lower()
            if phrase in line and not any(m in low for m in RETRACTION_MARKERS):
                out.append((i + 1, phrase, why, line.strip()[:160]))
    return out


@pytest.mark.parametrize("path", _markdown_files(), ids=lambda p: p.name)
def test_superseded_claims_are_never_asserted(path):
    """An unsupported claim may be quoted while being disclaimed, never stated as fact."""
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
