"""Partial-input mode (plan Phase 6b) - the equity-load-bearing feature.

The settings that most need disambiguation have the least complete structured input, so
usefulness cannot assume clean inputs. DISCERN runs on any subset of
{genetics, phenotype, labs} (the joint already handles empty streams) and quantifies what
each *missing* modality would add - making the cheapest decisive next observation the
feature that survives incomplete data.
"""
from __future__ import annotations

from core.dx_schemas import DiscriminationCluster
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import leading_disease, marginal_disease
from nextobs.recommend import _entropy

_MODALITIES = ("genetic_codes", "clinical", "functional")


def present_modalities(ev: Evidence) -> list[str]:
    return [m for m in _MODALITIES if getattr(ev, m)]


def missing_modalities(ev: Evidence) -> list[str]:
    return [m for m in _MODALITIES if not getattr(ev, m)]


def current_uncertainty(cluster: DiscriminationCluster, ev: Evidence) -> float:
    """Shannon entropy of the disease marginal under the available evidence."""
    return _entropy(marginal_disease(joint(cluster, ev)))


def diagnose_partial(cluster: DiscriminationCluster, ev: Evidence) -> dict:
    """Run on whatever subset of evidence exists; report present/missing modalities."""
    j = joint(cluster, ev)
    lead, p = leading_disease(j)
    return {
        "leading": lead, "p": round(p, 3),
        "present": present_modalities(ev), "missing": missing_modalities(ev),
        "entropy": round(_entropy(marginal_disease(j)), 4),
    }
