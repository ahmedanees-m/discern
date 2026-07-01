"""Build a PS1/PM5 lookup index from the kept ClinVar variant_summary (plan Step 1.5).

PS1 = same amino-acid change as a known pathogenic variant.
PM5 = a *different* pathogenic missense at the *same* residue.

We extract (gene, residue, ref_aa, pos, alt_aa) -> aggregate clinical significance
for GRCh38 missense variants classified Pathogenic / Likely pathogenic, and write a
compact JSON keyed by gene:
    {gene: {"Arg2228": {"Gln": "Pathogenic", "Trp": "Likely pathogenic"}}}

Run on the VM (where the data lives):
    python3 -m data.sources.build_clinvar_index \
        data/raw/clinvar/variant_summary.txt.gz data/processed/clinvar_ps1_pm5.json
"""
from __future__ import annotations

import gzip
import json
import re
import sys
from collections import defaultdict

# p.Arg2228Gln  /  p.Arg2326Ter  /  p.Gly12=
_PROT = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*|=)")
_PATH = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}


def parse_protein_change(name: str):
    """Return (ref_aa, pos, alt_aa) for a missense, else None."""
    m = _PROT.search(name or "")
    if not m:
        return None
    ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
    if alt in ("Ter", "*", "=") or alt == ref:        # nonsense / synonymous -> not missense
        return None
    return ref, pos, alt


def build(src: str, out: str) -> dict:
    # column indices from the variant_summary header (1-based in docs; 0-based here)
    idx = {"Name": 2, "GeneSymbol": 4, "ClinicalSignificance": 6, "Assembly": 16}
    index: dict[str, dict[str, dict[str, str]]] = defaultdict(lambda: defaultdict(dict))
    n_rows = n_path = 0
    op = gzip.open if src.endswith(".gz") else open
    with op(src, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline()  # noqa: F841 (skip header)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= idx["Assembly"]:
                continue
            n_rows += 1
            if cols[idx["Assembly"]] != "GRCh38":
                continue
            sig = cols[idx["ClinicalSignificance"]].strip()
            if sig.lower() not in _PATH:
                continue
            pc = parse_protein_change(cols[idx["Name"]])
            if not pc:
                continue
            ref, pos, alt = pc
            gene = cols[idx["GeneSymbol"]].strip()
            if not gene:
                continue
            residue = f"{ref}{pos}"
            # keep the strongest assertion if duplicated
            prev = index[gene][residue].get(alt)
            if prev != "Pathogenic":
                index[gene][residue][alt] = sig
            n_path += 1
    # materialize plain dicts
    plain = {g: {r: dict(alts) for r, alts in res.items()} for g, res in index.items()}
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(plain, fh)
    print(f"rows scanned: {n_rows}; P/LP missense indexed: {n_path}; "
          f"genes: {len(plain)} -> {out}")
    return plain


if __name__ == "__main__":
    build(sys.argv[1], sys.argv[2])
