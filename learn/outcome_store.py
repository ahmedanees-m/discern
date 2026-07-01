"""Outcome store (plan Step 5.2).

Records human accept/override decisions and realized assay outcomes, keyed by
(gene, assay), so per-institution priors can be updated auditably (each update
attributable to specific cases). In production this is the Postgres `db` service
(Phase 9); here it is an in-memory/JSON store with the same interface.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field

from learn.prior_update import BetaPrior


@dataclass
class Outcome:
    gene: str
    assay: str
    predicted_resolving: bool       # the engine recommended this action
    realized_resolved: bool         # the assay actually resolved the case
    accepted: bool                  # clinician accepted the recommendation
    case_id: str = ""


@dataclass
class OutcomeStore:
    outcomes: list[Outcome] = field(default_factory=list)

    def record(self, o: Outcome) -> None:
        self.outcomes.append(o)

    def assay_prior(self, gene: str, assay: str) -> BetaPrior:
        """Beta posterior on this (gene, assay)'s sensitivity from recorded outcomes."""
        prior = BetaPrior()
        for o in self.outcomes:
            if o.gene == gene and o.assay == assay:
                prior = prior.update(o.realized_resolved)
        return prior

    def acceptance_rate(self) -> float:
        return (sum(o.accepted for o in self.outcomes) / len(self.outcomes)
                if self.outcomes else 0.0)

    def to_json(self) -> str:
        return json.dumps([asdict(o) for o in self.outcomes])

    @classmethod
    def from_json(cls, s: str) -> OutcomeStore:
        return cls([Outcome(**d) for d in json.loads(s)])
