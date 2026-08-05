"""Sparse-LR uncertainty propagation.

Each disease-feature frequency is treated as Beta(freq-n+1, (1-freq)-n+1) with n = n_cases
(the provenance sample size). Monte-Carlo resampling of the frequencies propagates the
estimation uncertainty into a **credible interval** on P(D|E): ultra-rare diseases with
tiny n get wide intervals (-> abstention), so sparse data widens uncertainty rather than
fabricating confidence. (Hierarchical shrinkage is the rigorous form; MC is the v1.)
"""
from __future__ import annotations

import random
from dataclasses import replace

from core.dx_schemas import DiscriminationCluster
from jointdx.factorgraph import Evidence, joint
from jointdx.infer import marginal_disease


def _perturb(cluster: DiscriminationCluster, rng: random.Random) -> DiscriminationCluster:
    diseases = []
    for d in cluster.diseases:
        flr = {}
        for fid, entry in d.feature_lr.items():
            freq = float(entry[0])
            n = float(entry[1]) if len(entry) > 1 else 20.0
            a, b = freq * n + 1.0, (1 - freq) * n + 1.0
            flr[fid] = (rng.betavariate(a, b), n, entry[2] if len(entry) > 2 else "")
        diseases.append(replace(d, feature_lr=flr))
    return replace(cluster, diseases=diseases)


def disease_intervals(cluster: DiscriminationCluster, ev: Evidence, n_mc: int = 200,
                      seed: int = 0, alpha: float = 0.05) -> dict[str, tuple[float, float, float]]:
    """disease_id -> (mean, lo, hi) credible interval on P(D|E)."""
    rng = random.Random(seed)
    samples: dict[str, list[float]] = {d.id: [] for d in cluster.diseases}
    for _ in range(n_mc):
        md = marginal_disease(joint(_perturb(cluster, rng), ev))
        for k in samples:
            samples[k].append(md.get(k, 0.0))
    out = {}
    for k, vals in samples.items():
        vals.sort()
        n = len(vals)
        lo = vals[int(alpha / 2 * n)]
        hi = vals[min(n - 1, int((1 - alpha / 2) * n))]
        out[k] = (sum(vals) / n, lo, hi)
    return out
