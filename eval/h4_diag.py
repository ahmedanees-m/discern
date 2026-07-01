"""Diagnostic: reconcile DISCERN AUC with/without frequency + full-DB InterVar, on one multianno."""
import csv
import sys
from collections import Counter

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
        print(f"columns: revel={revel_c} gnomad={gnomad_c} func={func_c} exon={exon_c}")
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
            sv_on = score_variant(m["gene"], m["vid"], Annotations(af=af, revel=revel, consequence=cons))
            sv_off = score_variant(m["gene"], m["vid"], Annotations(af=None, revel=revel, consequence=cons))
            rows.append({"key": key, "label": lab, "d_on": sv_on.points, "d_off": sv_off.points,
                         "revel": revel, "af": af, "is_mis": "missense" in cons})
    iv = {}
    dist = Counter()
    with open(intervar_path, encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        order = {"benign": 0, "likely benign": 1, "uncertain significance": 2, "likely pathogenic": 3, "pathogenic": 4}
        ivcol = next((c for c in (rd.fieldnames or []) if "intervar" in c.lower()), None)
        chrc = next((c for c in (rd.fieldnames or []) if c.lower() in ("chr", "#chr")), "Chr")
        for r in rd:
            key = (str(r.get(chrc, "")).replace("chr", ""), str(r.get("Start", "")), r.get("Ref", ""), r.get("Alt", ""))
            txt = (r.get(ivcol, "") or "").lower()
            for name, val in sorted(order.items(), key=lambda kv: -len(kv[0])):
                if name in txt:
                    iv[key] = val
                    dist[name] += 1
                    break
    ys = [r["label"] for r in rows]
    mis = [r for r in rows if r["is_mis"]]
    rv = [r for r in rows if r["revel"] is not None]
    print(f"\nn_PB={len(rows)} pathogenic={sum(ys)} | missense={len(mis)} | revel_present={len(rv)}")
    print(f"AUC DISCERN af-ON  (all): {_auc([r['d_on'] for r in rows], ys):.4f}")
    print(f"AUC DISCERN af-OFF (all): {_auc([r['d_off'] for r in rows], ys):.4f}")
    print(f"AUC REVEL-alone   (present): {_auc([r['revel'] for r in rv], [r['label'] for r in rv]):.4f} (n={len(rv)})")
    pj = [(iv[r["key"]], r) for r in rows if r["key"] in iv]
    print(f"\nInterVar covered n={len(pj)}; class dist={dict(dist)}")
    print(f"AUC InterVar (full-DB):   {_auc([a for a, _ in pj], [r['label'] for _, r in pj]):.4f}")
    print(f"AUC DISCERN af-ON (IV subset): {_auc([r['d_on'] for _, r in pj], [r['label'] for _, r in pj]):.4f}")
    # missense subset
    mj = [(iv[r['key']], r) for r in mis if r['key'] in iv]
    if mj:
        print(f"\n[missense n={len(mj)}] InterVar {_auc([a for a,_ in mj],[r['label'] for _,r in mj]):.4f} | "
              f"DISCERN-on {_auc([r['d_on'] for _,r in mj],[r['label'] for _,r in mj]):.4f} | "
              f"REVEL {_auc([r['revel'] for _,r in mj if r['revel'] is not None] or [0],[r['label'] for _,r in mj if r['revel'] is not None] or [0]):.4f}")


if __name__ == "__main__":
    run(*sys.argv[1:4])
