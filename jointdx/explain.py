"""Templated explanation.

Composed deterministically from the joint result, the reclassification, the safety flags,
and the next observation. An LLM (Qwen via the gateway) may *phrase* this for readability
but never decides it - the verdict, the numbers, and the flags are all computed.
"""
from __future__ import annotations

from core.dx_schemas import DiscriminationCluster, SafetyFlag


def explain(cluster: DiscriminationCluster, p_disease: dict, reclass: tuple,
            flags: list[SafetyFlag], next_obs, decided: bool) -> str:
    names = {d.id: d.name for d in cluster.diseases}
    lead_id = max(p_disease, key=lambda k: p_disease[k][0]) if p_disease else ""
    lead_p = p_disease.get(lead_id, (0, 0, 0))
    parts = []
    if decided:
        parts.append(f"Leading: {names.get(lead_id, lead_id)} "
                     f"({lead_p[0]:.0%}; 95% CI {lead_p[1]:.0%}-{lead_p[2]:.0%}).")
    else:
        parts.append(f"Undecidable on current evidence (leading {names.get(lead_id, lead_id)} "
                     f"{lead_p[0]:.0%}) - acquire the deciding observation below.")
    old, new, _ = reclass
    if new != old:
        parts.append(f"Variant reclassified {old.name} -> {new.name} via the disease model.")
    else:
        parts.append(f"Variant remains {old.name} on current evidence.")
    for f in flags:
        parts.append(f"[flag] {f.message}")
    if next_obs is not None:
        parts.append(f"Cheapest decisive next step: {next_obs.name} "
                     f"({next_obs.kind}, accessibility={next_obs.accessibility}).")
    return " ".join(parts)
