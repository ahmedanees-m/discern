"""Mechanism-aware code -> action mapping.

Loads a domain action catalog (e.g. action_catalog/bleeding.yaml). The same code
(PS3) maps to different assays depending on the gene's disease *mechanism* - e.g.
an activation-defect gene (FERMT3/LAD-III) has NORMAL expression, so expression
flow cytometry must NOT be offered for it; the catalog routes to activation assays.
"""
from __future__ import annotations

import os
from functools import cache

import yaml

from core.schemas import Action

_CATALOG_DIR = os.path.join(os.path.dirname(__file__), "action_catalog")


@cache
def load_catalog(domain: str = "bleeding") -> dict:
    path = os.path.join(_CATALOG_DIR, f"{domain}.yaml")
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _action(name: str, spec: dict) -> Action:
    return Action(
        name=name,
        yields_codes=list(spec.get("yields", [])),
        cost_usd=float(spec.get("cost_usd", 0)),
        turnaround_days=int(spec.get("days", 0)),
        sensitivity=float(spec.get("sens", 0.0)),
        specificity=float(spec.get("spec", 0.0)),
        modality=str(spec.get("modality", "")),
        burden=float(spec.get("burden", 0.0)),
    )


def actions_for_code(code: str, mechanism: str, domain: str = "bleeding") -> list[Action]:
    """Concrete, mechanism-correct actions that can yield ``code`` (e.g. 'PS3').

    Falls back to a domain-wide search by yielded code when no mechanism entry
    exists (e.g. for segregation/trio actions that are mechanism-independent).
    """
    cat = load_catalog(domain)
    actions = cat.get("actions", {})
    names: list[str] = []
    mech = cat.get("mechanisms", {}).get(mechanism, {})
    names.extend(mech.get(code, []))
    if not names:
        # mechanism-independent fallback: any action whose `yields` contains the code
        names = [n for n, s in actions.items() if code in s.get("yields", [])]
    seen, out = set(), []
    for n in names:
        if n in actions and n not in seen:
            seen.add(n)
            out.append(_action(n, actions[n]))
    return out


def expand_actions(codes: list[str], mechanism: str, domain: str = "bleeding") -> list[Action]:
    """All distinct actions that can yield any of the requested codes."""
    seen, out = set(), []
    for code in codes:
        for a in actions_for_code(code, mechanism, domain):
            if a.name not in seen:
                seen.add(a.name)
                out.append(a)
    return out
