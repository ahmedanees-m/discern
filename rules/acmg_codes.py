"""Parse ACMG/AMP evidence-code strings into signed Tavtigian points.

ClinGen/eRepo writes applied criteria as ``CODE`` or ``CODE_Strength`` - e.g.
``PM2``, ``PP4_Moderate``, ``PVS1_Strong``, ``PM2_Supporting``, ``BS3_Supporting``,
``BA1``. A modifier overrides the code's default strength; otherwise the default
for the code's prefix applies.

Points (naturally-scaled, Tavtigian 2020): Supporting/Moderate/Strong/VeryStrong =
1/2/4/8, negated for benign. BA1 is a standalone benign pre-filter (flagged, not summed).
"""
from __future__ import annotations

import re

from core.schemas import EvidenceContribution, Strength

# strength label -> magnitude
_LABEL_POINTS = {
    "Standalone": 8,   # BA1 magnitude (handled as pre-filter, not summed)
    "VeryStrong": 8,
    "Strong": 4,
    "Moderate": 2,
    "Supporting": 1,
}

# default strength label by code prefix
_PREFIX_DEFAULT = {
    "PVS": "VeryStrong",
    "PS": "Strong",
    "PM": "Moderate",
    "PP": "Supporting",
    "BA": "Standalone",
    "BVS": "VeryStrong",
    "BS": "Strong",
    "BM": "Moderate",
    "BP": "Supporting",
}

# normalize free-text strength variants seen in the wild (eRepo writes e.g.
# "PM3_Very Strong" with a space, "Stand Alone", etc.). Keys are space/underscore-
# stripped, lowercased.
_LABEL_ALIASES = {
    "verystrong": "VeryStrong",
    "standalone": "Standalone",
    "strong": "Strong",
    "moderate": "Moderate",
    "supporting": "Supporting",
}


def _norm_label(s: str) -> str:
    return s.replace(" ", "").replace("_", "").lower()

_PREFIX_RE = re.compile(r"^(B?[A-Z]+?)(\d+)?$")   # PVS1 -> (PVS,1); PM2 -> (PM,2); BA1 -> (BA,1)


def _prefix(base: str) -> str:
    m = _PREFIX_RE.match(base)
    if not m:
        # fall back: leading alpha run
        return re.match(r"^[A-Za-z]+", base).group(0)
    return m.group(1)


def code_points(code: str) -> tuple[float, bool]:
    """Return (signed_points, is_ba1) for an eRepo-style code string.

    Unknown/empty codes return (0.0, False).
    """
    code = code.strip()
    if not code:
        return 0.0, False
    parts = code.split("_", 1)
    base = parts[0].strip()
    label = None
    if len(parts) == 2:
        label = _LABEL_ALIASES.get(_norm_label(parts[1]), parts[1].strip())
    prefix = _prefix(base).upper()
    if prefix not in _PREFIX_DEFAULT:
        return 0.0, False
    if label is None or label not in _LABEL_POINTS:
        label = _PREFIX_DEFAULT[prefix]
    is_ba1 = prefix == "BA"
    if is_ba1:
        return 0.0, True   # standalone benign: excluded from the sum (pre-filter)
    sign = 1 if prefix.startswith("P") else -1
    return float(sign * _LABEL_POINTS[label]), False


def _nearest_strength(points: float) -> Strength:
    mag = abs(points)
    label = min(_LABEL_POINTS, key=lambda k: abs(_LABEL_POINTS[k] - mag) if k != "Standalone" else 99)
    base = {"VeryStrong": "PVS", "Strong": "PS", "Moderate": "PM", "Supporting": "PP"}[label]
    name = base if points >= 0 else {"PVS": "BVS", "PS": "BS", "PM": "BM", "PP": "BP"}[base]
    return Strength[name]


def contributions_from_codes(codes: list[str]) -> tuple[list[EvidenceContribution], bool]:
    """Build EvidenceContributions from a list of code strings; also return ba1_met."""
    contribs: list[EvidenceContribution] = []
    ba1 = False
    for c in codes:
        pts, is_ba1 = code_points(c)
        if is_ba1:
            ba1 = True
            continue
        if pts == 0.0:
            continue
        contribs.append(
            EvidenceContribution(
                code=c.strip(),
                strength=_nearest_strength(pts),
                applied=True,
                source="eRepo",
                rationale="parsed from applied criteria",
            )
        )
    return contribs, ba1
