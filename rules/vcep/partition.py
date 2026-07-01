"""Per-code ACMG partition - the deep circularity fix (plan Phase 1.1 / Section A.1 / Gate G3).

A VCEP bottom-line label bundles PS3 (functional), PP4 (phenotype), PM3 (phasing),
PP1 (segregation) etc. - the very codes DISCERN owns in *other* factors. Consuming the
label would double-count them. Instead DISCERN decomposes the spec to per-code strengths
and routes **each code to exactly one owning factor**, so each code enters the joint
model once. VCEPs already organize into sub-teams of this exact shape, so the partition
is given.
"""
from __future__ import annotations

# code (base, strength-suffix stripped) -> owning factor.
FACTOR_OF: dict[str, str] = {}


def _register(factor: str, codes: list[str]) -> None:
    for c in codes:
        FACTOR_OF[c] = factor


# variant-intrinsic: everything about the variant itself (frequency, in-silico, null,
# same-residue) - enters P(E_geno | V) in the genetic factor (Phase 1.1).
_register("variant_intrinsic",
          ["PM2", "BS1", "BA1", "PP3", "BP4", "PVS1", "PM5", "PS1", "PM4", "BP3", "BP7",
           "PM1", "PS4", "PP2", "BP1", "PP5", "BP6", "BS2", "BP5"])
# PP4 -> the disease->variant coupling P(V|D) (Phase 3); never added in the genetic factor.
_register("disease_pp4", ["PP4"])
# functional -> P(E_func | D, V) (Phase 3); touches both disease and variant, counted once.
_register("functional", ["PS3", "BS3"])
# segregation / phasing / de-novo -> next-observation factors; enter only when performed.
_register("segregation", ["PP1", "BS4"])
_register("phasing", ["PM3", "BP2"])
_register("denovo", ["PS2", "PM6"])


def base_code(code: str) -> str:
    """Strip a strength suffix: 'PM2_Supporting' -> 'PM2'."""
    return code.split("_", 1)[0].strip()


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
