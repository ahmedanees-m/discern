"""Disease ontology + discrimination-cluster loader (plan Phase 2).

Loads the provenance-tagged discrimination clusters (six YAML files) into
``DiscriminationCluster`` objects, and routes a patient to the relevant cluster(s) by
gene and/or phenotype signature. Every LR/discriminator in the YAML carries a PMID +
sample size, so the knowledge base is curatable and citable (the ISTH/ClinGen-partnership
target). Multi-cluster routing is allowed.
"""
from __future__ import annotations

import os
from functools import lru_cache

import yaml

from core.dx_schemas import DiscriminationCluster, Disease, Observation

_CLUSTER_DIR = os.path.join(os.path.dirname(__file__), "clusters")


def _disease(d: dict) -> Disease:
    return Disease(
        id=d["id"], name=d["name"], genes=list(d.get("genes", [])),
        inheritance=d.get("inheritance", ""), mechanism=d.get("mechanism", ""),
        prior=float(d.get("prior", 0.1)),
        feature_lr={k: tuple(v) for k, v in d.get("feature_lr", {}).items()},
        p_path_given_disease=float(d.get("p_path_given_disease", 0.5)),
        treatment=d.get("treatment", ""), contraindications=list(d.get("contraindications", [])),
        vcep_spec=d.get("vcep_spec"),
    )


def _observation(o: dict) -> Observation:
    return Observation(
        id=o["id"], name=o["name"], kind=o.get("kind", "lab"),
        informs=list(o.get("informs", [])), outcome_lr=o.get("outcome_lr", {}),
        changes_management=bool(o.get("changes_management", False)),
        accessibility=o.get("accessibility", "moderate"),
    )


def load_cluster(path: str) -> DiscriminationCluster:
    with open(path, encoding="utf-8") as fh:
        d = yaml.safe_load(fh)
    return DiscriminationCluster(
        id=d["id"], name=d["name"],
        diseases=[_disease(x) for x in d.get("diseases", [])],
        discriminating_features=list(d.get("discriminating_features", [])),
        next_observations=[_observation(o) for o in d.get("next_observations", [])],
    )


@lru_cache(maxsize=1)
def all_clusters() -> dict[str, DiscriminationCluster]:
    out: dict[str, DiscriminationCluster] = {}
    for fn in sorted(os.listdir(_CLUSTER_DIR)):
        if fn.endswith(".yaml"):
            c = load_cluster(os.path.join(_CLUSTER_DIR, fn))
            out[c.id] = c
    return out


def route_clusters(genes: list[str] | None = None) -> list[DiscriminationCluster]:
    """Clusters whose diseases include any of the given genes (gene-based routing)."""
    genes = set(genes or [])
    hits = []
    for c in all_clusters().values():
        cluster_genes = {g for d in c.diseases for g in d.genes}
        if genes & cluster_genes:
            hits.append(c)
    return hits


def cluster_for(cluster_id: str) -> DiscriminationCluster:
    return all_clusters()[cluster_id]
