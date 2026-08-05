"""The swappability contract.

Every external tool / data source sits behind an ``EvidenceAdapter`` so it can be
replaced without touching the engine. Contract tests guarantee this. No primary
result of DISCERN may hard-depend on a closed/proprietary API.
"""
from __future__ import annotations

from abc import ABC, abstractmethod

from core.schemas import EvidenceContribution, PatientContext, Variant


class EvidenceAdapter(ABC):
    """Maps an ACMG evidence code (or code group) to a concrete source."""

    code_group: tuple[str, ...]   # e.g. ("PM2", "BS1", "BA1")
    version: str

    @abstractmethod
    def evaluate(self, v: Variant, p: PatientContext) -> list[EvidenceContribution]:
        """Return the evidence contributions this source supports for ``v``."""

    @abstractmethod
    def health_check(self) -> bool:
        """True if the source is reachable / required data is present."""
