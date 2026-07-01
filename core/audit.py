"""Immutable, hash-chained audit ledger (plan Step 5.1, gate G6).

Every classification, recommendation, and LLM call is appended here with a hash
linking it to the previous entry, so the full reasoning trace is tamper-evident
and reproducible. No PHI is written to the ledger in hosted/public contexts.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from typing import Any


def _hash(prev: str, payload: str) -> str:
    return hashlib.sha256((prev + payload).encode("utf-8")).hexdigest()


@dataclass
class AuditEntry:
    index: int
    kind: str                 # "classify" | "recommend" | "llm_call" | ...
    payload: dict[str, Any]
    prev_hash: str
    hash: str


@dataclass
class AuditLedger:
    entries: list[AuditEntry] = field(default_factory=list)

    GENESIS = "0" * 64

    def append(self, kind: str, payload: dict[str, Any]) -> AuditEntry:
        prev = self.entries[-1].hash if self.entries else self.GENESIS
        body = json.dumps(payload, sort_keys=True, default=str)
        entry = AuditEntry(
            index=len(self.entries),
            kind=kind,
            payload=payload,
            prev_hash=prev,
            hash=_hash(prev, f"{kind}|{body}"),
        )
        self.entries.append(entry)
        return entry

    def verify(self) -> bool:
        """Re-derive the chain and confirm no entry was altered."""
        prev = self.GENESIS
        for e in self.entries:
            body = json.dumps(e.payload, sort_keys=True, default=str)
            if e.prev_hash != prev or e.hash != _hash(prev, f"{e.kind}|{body}"):
                return False
            prev = e.hash
        return True
