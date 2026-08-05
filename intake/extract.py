"""Free-text -> HPO extractor with pertinent-negative capture.

The legitimate LLM soft task: turn a clinical narrative into HPO terms, capturing both
**present** and **explicitly excluded** features (pertinent negatives - "no leukocytosis"
becomes an absent term, which the joint model uses against LAD-III). Output is validated
against the HPO ontology and confirmed by a human before it informs the verdict; the rule
engine, not the LLM, assigns the evidence weight.

The extractor is injectable; the default routes through llm/gateway (cloud Nemotron / Qwen).
"""
from __future__ import annotations

from collections.abc import Callable

from core.dx_schemas import Feature, FeatureKind

# A single extracted item: (feature_id, present:bool, hpo_id:str)
ExtractedTerm = tuple[str, bool, str]


def _llm_extract(note: str) -> list[ExtractedTerm]:  # pragma: no cover - needs the LLM
    """Default extractor - prompts the gateway to emit present + explicitly-absent terms."""
    from llm.gateway import complete
    prompt = (
        "Extract clinical features from this note as JSON list of "
        '{"feature":..., "present":true/false, "hpo_id":...}. '
        "Include EXPLICITLY ABSENT findings (pertinent negatives) with present=false. "
        f"Note:\n{note}"
    )
    import json
    try:
        return [(d["feature"], bool(d["present"]), d.get("hpo_id", ""))
                for d in json.loads(complete("fast", prompt))]
    except Exception:
        return []


def extract_features(note: str,
                     extractor: Callable[[str], list[ExtractedTerm]] | None = None) -> list[Feature]:
    """Return Features (present + pertinent negatives) from a clinical note."""
    fn = extractor or _llm_extract
    out = []
    for fid, present, hpo_id in fn(note):
        out.append(Feature(id=fid, kind=FeatureKind.CLINICAL, value=present,
                           observed=present, source=f"note;hpo={hpo_id}"))
    return out
