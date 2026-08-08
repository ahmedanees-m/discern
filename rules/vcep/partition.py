"""Per-code ACMG partition - the deep circularity fix.

A VCEP bottom-line label bundles PS3 (functional), PP4 (phenotype), PM3 (phasing),
PP1 (segregation) etc. - the very codes DISCERN owns in *other* factors. Consuming the
label would double-count them. Instead DISCERN decomposes the spec to per-code strengths
and routes **each code to exactly one owning factor**, so each code enters the joint
model once. VCEPs already organize into sub-teams of this exact shape, so the partition
is given.
"""
from __future__ import annotations

import re

# code (base, strength-suffix stripped) -> owning factor.
FACTOR_OF: dict[str, str] = {}


def _register(factor: str, codes: list[str]) -> None:
    for c in codes:
        FACTOR_OF[c] = factor


# variant-intrinsic: everything about the variant itself (frequency, in-silico, null,
# same-residue) - enters P(E_geno | V) in the genetic factor.
_register("variant_intrinsic",
          ["PM2", "BS1", "BA1", "PP3", "BP4", "PVS1", "PM5", "PS1", "PM4", "BP3", "BP7",
           "PM1", "PS4", "PP2", "BP1", "PP5", "BP6", "BS2", "BP5"])
# PP4 -> the disease->variant coupling P(V|D); never added in the genetic factor.
_register("disease_pp4", ["PP4"])
# functional -> P(E_func | D, V); touches both disease and variant, counted once.
_register("functional", ["PS3", "BS3"])
# segregation / phasing / de-novo -> next-observation factors; enter only when performed.
_register("segregation", ["PP1", "BS4"])
_register("phasing", ["PM3", "BP2"])
_register("denovo", ["PS2", "PM6"])


def base_code(code: str) -> str:
    """Strip a strength suffix: 'PM2_Supporting' -> 'PM2'."""
    return code.split("_", 1)[0].strip()


# The applied-code vocabulary lives here, beside base_code, so every analysis that reads applied
# criteria shares one parser. A second parser anchored on word boundaries would drop every
# strength-modified code, since "_" is a word character; a single owner keeps the vocabulary and
# the factor map in step.
CODE_RE = re.compile(
    r"\b(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])"
    r"(_(?:Very[ _]Strong|Strong|Moderate|Supporting|Stand[ _]?Alone))?")

# An unmodified code carries its criterion's default strength, so "PP3" and "PP3_Supporting" are
# the same assertion. Comparing the literal strings makes identical strengths look like a
# disagreement, which is the same normalization mistake one level down.
_DEFAULT_STRENGTH = {"PVS1": "Very Strong", "BA1": "Stand Alone"}
_TIER_DEFAULT = {"PS": "Strong", "PM": "Moderate", "PP": "Supporting",
                 "BS": "Strong", "BP": "Supporting"}


def default_strength(code: str) -> str:
    """The ClinGen default strength an unmodified code carries."""
    base = base_code(code)
    return _DEFAULT_STRENGTH.get(base) or _TIER_DEFAULT[base[:2]]


def applied_codes(s: str) -> set[str]:
    """Base ACMG criteria in an applied-code string, with any strength modifier stripped."""
    return {m.group(1) for m in CODE_RE.finditer(s or "")}


def applied_codes_with_strength(s: str) -> dict[str, str]:
    """Criterion -> applied strength, with an unmodified code resolved to its default."""
    out = {}
    for m in CODE_RE.finditer(s or ""):
        mod = (m.group(2) or "").lstrip("_").replace("_", " ").strip()
        out[m.group(1)] = mod or default_strength(m.group(1))
    return out


def owner(code: str) -> str | None:
    """The single factor that owns this ACMG code, or None if unknown."""
    return FACTOR_OF.get(base_code(code))


def is_variant_intrinsic(code: str) -> bool:
    return owner(code) == "variant_intrinsic"


def partition(codes: list[str]) -> dict[str, list[str]]:
    """Group codes by owning factor (for the circularity audit / VCEP reconstruction)."""
    out: dict[str, list[str]] = {}
    for c in codes:
        f = owner(c) or "unknown"
        out.setdefault(f, []).append(c)
    return out
