"""Load gene-specific VCEP specifications into a VCEPSpec (plan Step 0.2 / 0.4).

A VCEPSpec carries the gene-specific parameters the engine needs: allele-frequency
thresholds (PM2/BS1/BA1), the prior probability, disease mechanism (for action
routing), inheritance, and the BA1 standalone pre-filter. Specs are versioned YAML
in ``rules/specs/`` so any VCEP rule set is swappable without touching the engine.

Threshold provenance: gene-specific numbers should come from the published ClinGen
CSpec spec (e.g. F8 = GN071). Where a number is not yet extracted from the spec PDF
it is a clearly-labeled general default; the YAML records this in ``notes``.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from functools import cache

import yaml

from core.schemas import PointsLedger

_SPEC_DIR = os.path.join(os.path.dirname(__file__), "specs")


@dataclass
class VCEPSpec:
    gene: str
    version: str
    inheritance: str = ""
    mechanism: str = ""
    prior_p: float = 0.10
    af: dict = field(default_factory=lambda: {"ba1": 0.05, "bs1": 0.0015, "pm2": 0.0001})
    notes: str = ""

    # AF accessors take an optional gene arg for forward-compat with multi-gene specs.
    def ba1_threshold(self, gene: str | None = None) -> float:
        return float(self.af["ba1"])

    def bs1_threshold(self, gene: str | None = None) -> float:
        return float(self.af["bs1"])

    def pm2_threshold(self, gene: str | None = None) -> float:
        return float(self.af["pm2"])

    def ba1_met(self, ledger: PointsLedger) -> bool:
        """Standalone benign pre-filter: BA1 applied in the ledger (set by gnomAD adapter)."""
        return any(c.applied and c.code.split("_", 1)[0] == "BA1"
                   for c in ledger.contributions)


def load_spec(path: str) -> VCEPSpec:
    with open(path, encoding="utf-8") as fh:
        d = yaml.safe_load(fh)
    return VCEPSpec(
        gene=d["gene"], version=d["version"],
        inheritance=d.get("inheritance", ""), mechanism=d.get("mechanism", ""),
        prior_p=float(d.get("prior_p", 0.10)),
        af=d.get("af", {"ba1": 0.05, "bs1": 0.0015, "pm2": 0.0001}),
        notes=d.get("notes", ""),
    )


@cache
def get_spec(gene: str) -> VCEPSpec:
    """Load the gene's spec if present, else the general base ACMG defaults."""
    path = os.path.join(_SPEC_DIR, f"{gene}.yaml")
    if not os.path.exists(path):
        path = os.path.join(_SPEC_DIR, "base_acmg.yaml")
    return load_spec(path)
