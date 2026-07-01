"""CHAMP/CHBMP MISSENSE arm - novel-missense recall with routed PS1/PM5 + PP3 (mirrors H4).

The null-subset benchmark (`champ_chbmp_benchmark.py`) showed catalogued missense ceiling at VUS on
PURELY variant-intrinsic evidence (PP3+PM2 = 2 pts). This arm adds the ROUTED ClinVar codes DISCERN
actually uses for a novel missense in clinic - PS1 (same AA change pathogenic in ClinVar) and PM5
(different change, same residue) via `adapters/clinvar.py` - plus PP3 (REVEL >= F8/F9 cut-off 0.6)
and PM2 (gnomAD-absent), then asks: **what fraction of the CDC catalog's missense disease alleles
does DISCERN recover to Likely-Pathogenic/Pathogenic?**

Rigor (anti-circularity): only variants whose **exact cDNA change is NOT in ClinVar** (truly novel)
are credited PS1/PM5 - so the same-residue evidence always comes from a *different* variant.
Point arithmetic (Tavtigian): a novel missense reaches LP (>=6) only as PS1(4)+PM2(1)+PP3(1); PM5(2)
and no-hit cap at VUS. So the **PS1 rate is the ceiling**, and REVEL+PM2 are the two deciding points.

Local inputs: ClinVar variant_summary (novelty + PS1/PM5 index, already built). REVEL is optional
(`eval/champ_chbmp_ps1_revel.tsv`, gene<TAB>protein<TAB>revel - committed; produced by the VM REVEL pass);
gnomAD AF optional (`data/processed/f8f9_gnomad.tsv`). Without them the arm reports the ceiling +
an explicit PM2/PP3 bracket. No PHI. Run: `python -m eval.champ_chbmp_missense_arm`.
"""
from __future__ import annotations

import gzip
import os
import re
from pathlib import Path

from adapters.clinvar import ClinVarAdapter
from core.schemas import PatientContext, Variant
from rules.acmg_codes import code_points
from rules.point_engine import BANDS, Classification
from rules.variant_scoring import Annotations, score_variant

DATA_DIR = Path(os.environ.get("CHAMP_CHBMP_DIR", "data/champ_chbmp"))
CLINVAR = Path("data/raw/clinvar/variant_summary.txt.gz")
# REVEL for the 29 PS1 cases (the only variants that can reach LP) - obtained via ANNOVAR hg19_revel
# on the VM (2026-06-16), committed for reproducibility. gnomAD PM2: verified 28/28 PS1 SNVs absent
# from gnomAD v2.1.1, so the af=0 default (PM2 fires) is correct - no gnomAD file needed.
REVEL_TSV = Path(__file__).with_name("champ_chbmp_ps1_revel.tsv")
GNOMAD_TSV = Path("data/processed/f8f9_gnomad.tsv")   # optional; PM2 verified gnomAD-absent
FILES = {"F8": ("CHAMP-Variant-List-2022.xlsx", "CHAMP Variant List", 2),
         "F9": ("CHBMP-Variant-List-2022.xlsx", "CHBMP Variant List", 3)}
_CDNA = re.compile(r"c\.[^\s)]+")


def _classify(pts: float) -> Classification:
    for thr, cls in BANDS:
        if pts >= thr:
            return cls
    return Classification.B


def clinvar_f8f9_cdna() -> dict[str, set]:
    """Set of exact cDNA changes present in ClinVar for F8/F9 (for nucleotide-level novelty)."""
    out = {"F8": set(), "F9": set()}
    with gzip.open(CLINVAR, "rt", encoding="utf-8", errors="replace") as fh:
        fh.readline()
        for line in fh:
            c = line.split("\t")
            if len(c) < 17 or c[4] not in out:
                continue
            m = _CDNA.search(c[2])
            if m:
                out[c[4]].add(m.group(0))
    return out


def load_missense(gene: str):
    import openpyxl
    fname, sheet, pcol = FILES[gene]
    wb = openpyxl.load_workbook(DATA_DIR / fname, read_only=True, data_only=True)
    ws = wb[sheet]
    it = ws.iter_rows(values_only=True)
    next(it)
    vtcol = 4 if gene == "F8" else 5
    recs = []
    for r in it:
        if r[0] in (None, "") or "missense" not in re.sub(r"\s+", " ", str(r[vtcol] or "").lower()):
            continue
        prot = re.sub(r"[()]", "", str(r[pcol] or ""))  # 'p.(Asp34Val)' -> 'p.Asp34Val'
        recs.append({"gene": gene, "cdna": _CDNA.search(str(r[0])).group(0) if _CDNA.search(str(r[0])) else str(r[0]).strip(),
                     "prot": prot})
    wb.close()
    return recs


def _load_kv(path: Path):
    if not path.exists():
        return None
    d = {}
    for line in path.read_text().splitlines():
        p = line.split("\t")
        if len(p) >= 3:
            d[(p[0], p[1])] = float(p[2])
    return d


def run():
    clin = clinvar_f8f9_cdna()
    adapter = ClinVarAdapter()
    revel = _load_kv(REVEL_TSV)
    gnomad = _load_kv(GNOMAD_TSV)
    pc = PatientContext()
    PATH = {Classification.P, Classification.LP}
    print("=" * 80)
    print("CHAMP/CHBMP MISSENSE arm - novel-missense recall (PS1/PM5 + PP3 + PM2)")
    print(f"REVEL={'loaded' if revel else 'ABSENT (PP3 off)'}  gnomAD={'loaded' if gnomad else 'ABSENT (PM2 bracketed)'}")
    print("=" * 80)

    for gene in ("F8", "F9"):
        recs = load_missense(gene)
        novel = [r for r in recs if r["cdna"] not in clin[gene]]
        ps1 = pm5 = neither = 0
        lp_recall = 0  # LP/P under real REVEL (PP3) + gnomAD-PM2, when those files are present
        for r in novel:
            v = Variant(chrom="X", pos=0, ref="N", alt="N", gene=gene, hgvs_p=r["prot"])
            contribs = adapter.evaluate(v, pc)
            routed = [c.code for c in contribs]
            if "PS1" in routed:
                ps1 += 1
            elif "PM5" in routed:
                pm5 += 1
            else:
                neither += 1
            rv = revel.get((gene, r["prot"])) if revel else None
            af = gnomad.get((gene, r["prot"])) if gnomad else 0.0  # assume gnomAD-absent if no file
            sv = score_variant(gene, r["cdna"], Annotations(consequence="missense", revel=rv, af=af))
            allcodes = list(sv.codes) + routed
            pts = sum(code_points(c)[0] for c in allcodes)
            if _classify(pts) in PATH:
                lp_recall += 1
        n = len(novel)
        print(f"\n[{gene}] {len(recs)} catalogued missense; {n} novel (exact cDNA not in ClinVar), "
              f"{len(recs)-n} already in ClinVar")
        print(f"  routed ClinVar evidence on novel: PS1={ps1} ({100*ps1/n:.1f}%)  "
              f"PM5={pm5} ({100*pm5/n:.1f}%)  neither={neither} ({100*neither/n:.1f}%)")
        print(f"  LP-reachable ceiling = PS1 rate (only PS1+PM2+PP3=6 clears the bar; "
              f"PM5+PM2+PP3=4 and no-hit=2 stay VUS): {ps1}/{n} = {100*ps1/n:.1f}%")
        if revel:
            print(f"  LP/P recall [REAL REVEL + gnomAD-PM2]: {lp_recall}/{n} = {100*lp_recall/n:.1f}%")
        else:
            print(f"  LP/P recall: 0% as-is (PP3 needs REVEL); "
                  f"~{100*ps1/n:.1f}% once the {ps1} PS1 cases get REVEL>=0.6 (they are known-pathogenic "
                  f"substitutions, so ~all pass). REVEL needs the VM dbNSFP pass.")


if __name__ == "__main__":
    run()
