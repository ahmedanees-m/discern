"""Genetic-similarity (ancestry) inference.

Infers a genetic-similarity group (NOT self-reported race/ethnicity) by projecting the
sample onto reference principal components, plus a confidence. The production path uses
``somalier``/``peddy`` over the VCF against a reference PC panel - these run inside a
tool image on the VM (Docker-only); this module defines the interface + result type and
a pure-Python projector used when PCs are already extracted, so the engine stays
testable without the tool.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

# gnomAD v4 genetic-ancestry group labels (genetic-similarity groups).
GROUPS = ("afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas", "remaining")


@dataclass
class AncestryResult:
    group: str
    confidence: float            # 0..1
    pcs: tuple[float, ...] = ()


def assign_group(sample_pcs: list[float], reference_centroids: dict[str, list[float]]) -> AncestryResult:
    """Nearest-centroid assignment in PC space + a softmax-style confidence.

    `reference_centroids` maps group -> centroid PCs (from the reference panel). This is
    the deterministic core; PC extraction from a VCF is done by somalier/peddy upstream.
    """
    def d2(a, b):
        return sum((x - y) ** 2 for x, y in zip(a, b, strict=False))

    dists = {g: d2(sample_pcs, c) for g, c in reference_centroids.items()}
    nearest = min(dists, key=dists.get)
    # confidence = relative closeness vs the runner-up (1 - softmin share)
    inv = {g: math.exp(-v) for g, v in dists.items()}
    total = sum(inv.values()) or 1.0
    return AncestryResult(group=nearest, confidence=inv[nearest] / total,
                          pcs=tuple(sample_pcs))


def infer_from_vcf(vcf_path: str) -> AncestryResult:  # pragma: no cover - requires tool image
    """Production entry point - runs somalier/peddy in a tool image, then assign_group.

    Not implemented in-process (Docker-only on the VM). Wire via deploy/remote.py to a
    container that emits PCs, then call assign_group().
    """
    raise NotImplementedError(
        "ancestry inference from a VCF runs in a somalier/peddy tool image (Docker-only); "
        "extract PCs there, then call assign_group(sample_pcs, reference_centroids)."
    )
