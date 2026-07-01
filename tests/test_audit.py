"""Hash-chained audit ledger integrity (plan Step 5.1)."""
from core.audit import AuditLedger


def test_chain_verifies():
    led = AuditLedger()
    led.append("classify", {"variant": "X-100-A-T", "class": "VUS", "points": 1})
    led.append("recommend", {"action": "flow_cd41_cd61", "cost": 450})
    assert led.verify()


def test_tamper_detected():
    led = AuditLedger()
    led.append("classify", {"points": 1})
    led.append("recommend", {"action": "rna_seq"})
    led.entries[0].payload["points"] = 99   # tamper after the fact
    assert not led.verify()
