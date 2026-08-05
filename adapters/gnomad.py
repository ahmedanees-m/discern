"""gnomAD adapter - PM2 / BS1 / BA1.

Ancestry-resolved population frequency, mapped to ACMG codes via the gene's VCEP
allele-frequency thresholds:
  * AF above BA1 threshold  -> BA1 (standalone benign pre-filter)
  * AF above BS1 threshold  -> BS1 (Strong benign)
  * absent / below PM2      -> PM2 (Supporting, per 2020 ClinGen SVI downgrade)

Frequency is the 95% filtering allele frequency popmax (grpmax) when available, the
worst case across genome + exome - this is the ClinGen-preferred statistic. The
network query (gnomAD GraphQL API) is injectable as ``af_lookup`` for offline tests.
"""
from __future__ import annotations

import json
import urllib.request
from collections.abc import Callable

from core.adapter import EvidenceAdapter
from core.schemas import EvidenceContribution, PatientContext, Strength, Variant

_API = "https://gnomad.broadinstitute.org/api"
_QUERY = """
query Variant($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    genome { af faf95 { popmax } }
    exome  { af faf95 { popmax } }
  }
}
"""


def _gnomad_faf(variant_id: str, dataset: str = "gnomad_r4", timeout: int = 30) -> float | None:
    """Return the worst-case 95% filtering AF popmax (fallback to af), or None if absent."""
    body = json.dumps({"query": _QUERY,
                       "variables": {"variantId": variant_id, "dataset": dataset}}).encode()
    req = urllib.request.Request(_API, data=body,
                                 headers={"Content-Type": "application/json",
                                          "User-Agent": "discern"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = json.load(r)
    v = (data.get("data") or {}).get("variant")
    if not v:
        return None
    vals = []
    for src in ("genome", "exome"):
        s = v.get(src)
        if not s:
            continue
        faf = (s.get("faf95") or {}).get("popmax")
        vals.append(faf if faf is not None else s.get("af"))
    vals = [x for x in vals if x is not None]
    return max(vals) if vals else None


class GnomadAdapter(EvidenceAdapter):
    code_group = ("PM2", "BS1", "BA1")
    version = "v4.1"

    def __init__(self, spec, af_lookup: Callable[[str], float | None] | None = None):
        self.spec = spec
        self._af_lookup = af_lookup or _gnomad_faf

    def health_check(self) -> bool:
        if self._af_lookup is not _gnomad_faf:
            return True
        try:
            # a well-known common variant should return a frequency
            return self._af_lookup("1-55051215-G-A") is not None
        except Exception:
            return False

    def evaluate(self, v: Variant, p: PatientContext) -> list[EvidenceContribution]:
        variant_id = f"{v.chrom.replace('chr', '')}-{v.pos}-{v.ref}-{v.alt}"
        try:
            af = self._af_lookup(variant_id)
        except Exception:
            return []
        src = f"gnomAD {self.version}"
        gene = v.gene
        if af is not None and af > self.spec.ba1_threshold(gene):
            return [EvidenceContribution("BA1", Strength.BVS, True, src,
                    rationale=f"grpmax FAF {af:.3g} > BA1 {self.spec.ba1_threshold(gene):.3g}")]
        if af is not None and af > self.spec.bs1_threshold(gene):
            return [EvidenceContribution("BS1", Strength.BS, True, src,
                    rationale=f"grpmax FAF {af:.3g} > BS1 {self.spec.bs1_threshold(gene):.3g}")]
        if af is None or af < self.spec.pm2_threshold(gene):
            shown = "absent" if af is None else f"{af:.3g}"
            return [EvidenceContribution("PM2", Strength.PP, True, src,
                    rationale=f"grpmax FAF {shown} < PM2 {self.spec.pm2_threshold(gene):.3g}")]
        return []
