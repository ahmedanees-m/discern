"""Clean full-DB InterVar vs DISCERN on the H4 set - the headline-bulletproofing comparison.

Reuses h4_full's validated AUC / label / consequence helpers and the same score_variant path, but
does NOT depend on the AlphaMissense file (this comparison is DISCERN-full vs literal InterVar only).
InterVar here is run end-to-end with its FULL default database set (refGene, ensGene, knownGene,
esp6500siv2_all, 1000g2015aug_all, avsnp147, dbnsfp42a, dbscsnv11, clinvar, rmsk, full intervardb/
OMIM) - the prior reduced-config gaps (dropped 1000g/esp/avsnp, empty mim2gene) are closed; the sole
deviation from InterVar's literal default is gnomAD exome (coding-appropriate) in the gnomad_genome
slot. AUC is computed on the same P/B-labelled variants InterVar covers (apples-to-apples).

Run on the VM: python3 -m eval.intervar_full_eval h4set.meta.tsv <fixed_multianno> <intervar_out>
"""
import csv
import sys

from eval.h4_full import _auc, _consequence, _f, _label
from rules.variant_scoring import Annotations, score_variant


def run(meta_path, anno_path, intervar_path):
    meta = {}
    with open(meta_path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta[(r["chrom"].replace("chr", ""), r["pos"], r["ref"], r["alt"])] = r

    rows = []
    with open(anno_path, encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames or []
        revel_c = next((c for c in cols if c.lower() == "revel_score"), None)
        func_c = next((c for c in cols if c.lower().startswith("func.refgene")), None)
        exon_c = next((c for c in cols if c.lower().startswith("exonicfunc.refgene")), None)
        gnomad_c = next((c for c in cols if "gnomad" in c.lower() and c.lower().endswith("_all")), None)
        for r in rd:
            key = (str(r["Chr"]).replace("chr", ""), str(r["Start"]), r["Ref"], r["Alt"])
            m = meta.get(key)
            if not m:
                continue
            lab = _label(m["clnsig"])
            if lab is None:
                continue
            revel = _f(r.get(revel_c)) if revel_c else None
            af = _f(r.get(gnomad_c)) if gnomad_c else None
            cons = _consequence(r.get(func_c, ""), r.get(exon_c, ""))
            sv = score_variant(m["gene"], m["vid"], Annotations(af=af, revel=revel, consequence=cons))
            rows.append({"key": key, "label": lab, "discern": sv.points})

    iv = {}
    with open(intervar_path, encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        order = {"benign": 0, "likely benign": 1, "uncertain significance": 2,
                 "likely pathogenic": 3, "pathogenic": 4}
        ivcol = next((c for c in (rd.fieldnames or []) if "intervar" in c.lower()), None)
        chrc = next((c for c in (rd.fieldnames or []) if c.lower() in ("chr", "#chr")), "Chr")
        for r in rd:
            key = (str(r.get(chrc, "")).replace("chr", ""), str(r.get("Start", "")), r.get("Ref", ""), r.get("Alt", ""))
            txt = (r.get(ivcol, "") or "").lower()
            for name, val in sorted(order.items(), key=lambda kv: -len(kv[0])):
                if name in txt:
                    iv[key] = val
                    break

    ys = [r["label"] for r in rows]
    paired = [(iv[r["key"]], r["label"], r["discern"]) for r in rows if r["key"] in iv]
    print(f"P/B variants scored: {len(rows)}  (pathogenic={sum(ys)})")
    print(f"AUC DISCERN-full (all P/B):           {_auc([r['discern'] for r in rows], ys):.4f}")
    print(f"InterVar-covered subset:              n={len(paired)}")
    print(f"AUC InterVar (full-DB, ordinal):      {_auc([a for a, _, _ in paired], [b for _, b, _ in paired]):.4f}")
    print(f"AUC DISCERN (same paired subset):     {_auc([d for _, _, d in paired], [b for _, b, _ in paired]):.4f}")


if __name__ == "__main__":
    run(*sys.argv[1:4])
