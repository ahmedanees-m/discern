"""Vignette harness for with-and-without evaluation.

Loads the vignette bank in `eval/cases/reader_vignettes.yaml`, builds the card the engine would
present for a case (ranked diagnosis, recommended deciding test, any management hard stop), and
measures the engine-side ceiling: whether it identifies the disease, recommends the deciding
test, and raises the safety flag on treatment-divergent pairs. Vignettes are citation-only and
carry no patient identifiers.

Run:  python3 -m eval.reader_study
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import yaml

from core.dx_schemas import Feature, FeatureKind
from diseases.ontology import route_clusters
from jointdx.factorgraph import Evidence
from jointdx.orchestrate import diagnose

_BANK = os.path.join(os.path.dirname(__file__), "cases", "reader_vignettes.yaml")


@dataclass
class Vignette:
    id: str
    ev: Evidence
    correct_disease_id: str
    correct_next_obs_id: str | None = None
    harmful_tx: str | None = None        # a treatment that must be flagged/avoided
    tier: str = "T"                      # S=safety-divergent, T=deciding-test, D=acquired distractor
    gene: str = ""
    cluster: str = ""                    # the declared confusable cluster (for validity checks)
    scenario: str = ""
    source_pmid: str = ""


@dataclass
class ReaderScore:
    n: int                               # total vignettes in the bank
    n_scored: int                        # vignettes DISCERN routes (has a causal gene in a cluster)
    n_skipped: int                       # gene-less distractors DISCERN correctly does not route
    disease_accuracy: float              # over n_scored
    next_test_accuracy: float            # over scored vignettes that name a correct deciding test
    harmful_avoidance: float             # fraction of safety vignettes where the hazard was flagged
    rows: list = field(default_factory=list)


def _evidence(v: dict) -> Evidence:
    # Same present/absent convention as eval/curated_case_benchmark.py and evidence/phenotype_lr.py:
    # `observed` is the present(True)/pertinent-negative(False) flag the phenotype LR reads.
    clin = [Feature(fid, FeatureKind.LAB, bool(present), observed=bool(present))
            for fid, present in v.get("features", {}).items()]
    return Evidence(variant_gene=v.get("gene", "") or "", clinical=clin)


def load_vignettes(path: str = _BANK) -> list[Vignette]:
    raw = yaml.safe_load(open(path, encoding="utf-8"))["vignettes"]
    out = []
    for v in raw:
        out.append(Vignette(
            id=v["id"], ev=_evidence(v), correct_disease_id=v["true_dx"],
            correct_next_obs_id=v.get("correct_next_obs"), harmful_tx=v.get("harmful_tx"),
            tier=v.get("tier", "T"), gene=v.get("gene", "") or "", cluster=v.get("cluster", ""),
            scenario=" ".join((v.get("scenario") or "").split()), source_pmid=v.get("source_pmid", "")))
    return out


def discern_card(v: Vignette, n_mc: int = 80) -> dict | None:
    """The DISCERN-side card a reader sees in the aided arm. None when DISCERN does not route
    (gene-less acquired mimic) - the correct behaviour is to defer to the clinician, not over-call."""
    if not v.ev.variant_gene or not route_clusters([v.ev.variant_gene]):
        return None
    rec = diagnose(v.ev, planned_tx=v.harmful_tx, n_mc=n_mc)
    if rec is None:
        return None
    return {
        "id": v.id,
        "leading": rec.posterior.leading,
        "confidence": rec.posterior.confidence,
        "decided": rec.posterior.decided,
        "next_obs": rec.next_observation.id if rec.next_observation is not None else None,
        "n_high_flags": sum(1 for f in rec.safety_flags if f.severity == "high"),
        "flag_msgs": [f.message for f in rec.safety_flags],
    }


def score_discern(vignettes: list[Vignette], n_mc: int = 80) -> ReaderScore:
    dx = nt = harm = 0
    nt_total = harm_total = n_scored = 0
    rows = []
    for v in vignettes:
        card = discern_card(v, n_mc=n_mc)
        if card is None:                 # gene-less distractor: not routed (correct)
            rows.append({"id": v.id, "tier": v.tier, "routed": False})
            continue
        n_scored += 1
        dx_ok = card["leading"] == v.correct_disease_id
        dx += dx_ok
        nt_ok = None
        if v.correct_next_obs_id is not None:
            nt_total += 1
            nt_ok = card["next_obs"] == v.correct_next_obs_id
            nt += nt_ok
        harm_ok = None
        if v.harmful_tx is not None:
            harm_total += 1
            harm_ok = (card["n_high_flags"] > 0
                       or any(v.harmful_tx in m.lower() for m in card["flag_msgs"]))
            harm += harm_ok
        rows.append({"id": v.id, "tier": v.tier, "routed": True, "lead": card["leading"],
                     "dx_ok": dx_ok, "next_obs": card["next_obs"], "nt_ok": nt_ok,
                     "harmful_tx": v.harmful_tx, "harm_flagged": harm_ok})
    n = len(vignettes)
    return ReaderScore(
        n=n, n_scored=n_scored, n_skipped=n - n_scored,
        disease_accuracy=dx / n_scored if n_scored else 0.0,
        next_test_accuracy=nt / nt_total if nt_total else 0.0,
        harmful_avoidance=harm / harm_total if harm_total else 1.0,
        rows=rows)


def main():
    vs = load_vignettes()
    s = score_discern(vs)
    print(f"Reader-study vignette bank: n={s.n}  routed(scored)={s.n_scored}  "
          f"abstained-distractors={s.n_skipped}")
    print(f"DISCERN-side ceiling:  disease={s.disease_accuracy:.0%}  "
          f"next-test={s.next_test_accuracy:.0%}  harmful-avoidance={s.harmful_avoidance:.0%}")
    for r in s.rows:
        if not r["routed"]:
            print(f"  [D abstain ] {r['id']:20} (gene-less mimic; DISCERN does not route)")
            continue
        dxm = "OK " if r["dx_ok"] else "MISS"
        ntm = "-" if r["nt_ok"] is None else ("OK" if r["nt_ok"] else "no")
        hzm = "-" if r["harm_flagged"] is None else ("FLAGGED" if r["harm_flagged"] else "MISSED")
        print(f"  [{r['tier']} dx:{dxm}] {r['id']:20} lead={r['lead']:14} "
              f"next-test={ntm:2} ({r['next_obs']})  hazard={hzm}")


if __name__ == "__main__":
    main()
