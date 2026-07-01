"""Full H4 benchmark with REAL predictors (v3.1 Track A2 / H4).

Joins the ANNOVAR-annotated variant set (REVEL + consequence from refGene/dbNSFP) with the
AlphaMissense lookup and the ClinVar binary label, runs DISCERN score_variant with the REAL
per-variant predictors (REVEL + AlphaMissense + PVS1-from-consequence), and compares AUC vs
ClinVar against REVEL-alone and AlphaMissense-alone. Optionally folds in an InterVar
classification column if present (InterVar's ACMG class -> ordinal score).

Run on the VM:
    python3 -m eval.h4_full /tmp/h4set.meta.tsv /tmp/am.tsv /tmp/anno.hg38_multianno.txt [intervar.txt]
"""
from __future__ import annotations

import csv
import sys

from rules.variant_scoring import Annotations, score_variant

_LOF = {"stopgain": "nonsense", "frameshift_deletion": "frameshift",
        "frameshift_insertion": "frameshift", "stoploss": "stop_lost",
        "startloss": "start_lost"}


def _label(clnsig: str):
    t = clnsig.lower()
    if "conflict" in t or "uncertain" in t:
        return None
    if "pathogenic" in t and "benign" not in t:
        return 1
    if "benign" in t and "pathogenic" not in t:
        return 0
    return None


def _consequence(func: str, exonic: str) -> str:
    f, e = (func or "").lower(), (exonic or "").lower()
    if "splic" in f:
        return "canonical_splice"
    for k, v in _LOF.items():
        if k in e:
            return v
    if "nonsynonymous" in e:
        return "missense"
    return e or f


def _auc(scores, ys):
    pos = [s for s, y in zip(scores, ys, strict=False) if y == 1]
    neg = [s for s, y in zip(scores, ys, strict=False) if y == 0]
    if not pos or not neg:
        return None
    wins = sum((p > q) + 0.5 * (p == q) for p in pos for q in neg)
    return round(wins / (len(pos) * len(neg)), 4)


def _f(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def run(meta_path, am_path, anno_path, intervar_path=None):
    meta = {}
    with open(meta_path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta[(r["chrom"].replace("chr", ""), r["pos"], r["ref"], r["alt"])] = r
    am = {}
    with open(am_path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            am[r["vid"]] = _f(r["am_pathogenicity"])

    # ANNOVAR multianno (tab); columns by name
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
            amv = am.get(m["vid"])
            sv = score_variant(m["gene"], m["vid"],
                               Annotations(af=af, revel=revel, alphamissense=amv, consequence=cons))
            rows.append({"key": key, "label": lab, "discern": sv.points,
                         "revel": revel, "am": amv})

    iv = {}
    if intervar_path:
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
    out = {"n": len(rows), "n_pathogenic": sum(ys),
           "auc_discern_full": _auc([r["discern"] for r in rows], ys)}
    rv = [r for r in rows if r["revel"] is not None]
    av = [r for r in rows if r["am"] is not None]
    out["auc_revel"] = _auc([r["revel"] for r in rv], [r["label"] for r in rv])
    out["n_revel"] = len(rv)
    out["auc_alphamissense"] = _auc([r["am"] for r in av], [r["label"] for r in av])
    out["n_am"] = len(av)
    # fair head-to-head: DISCERN-full restricted to the SAME variants each predictor covers
    out["auc_discern_on_revel_subset"] = _auc([r["discern"] for r in rv], [r["label"] for r in rv])
    out["auc_discern_on_am_subset"] = _auc([r["discern"] for r in av], [r["label"] for r in av])
    if iv:
        paired = [(iv[r["key"]], r["label"]) for r in rows if r["key"] in iv]
        out["auc_intervar"] = _auc([a for a, _ in paired], [b for _, b in paired])
        out["n_intervar"] = len(paired)
    return out


def main():
    iv = sys.argv[4] if len(sys.argv) > 4 else None
    s = run(sys.argv[1], sys.argv[2], sys.argv[3], iv)
    print(f"Variants with binary ClinVar label: {s['n']} ({s['n_pathogenic']} pathogenic)")
    print(f"  AUC DISCERN-full, all variants:        {s['auc_discern_full']}  (n={s['n']})")
    print(f"  -- on REVEL subset:  DISCERN {s['auc_discern_on_revel_subset']} vs REVEL {s['auc_revel']}  (n={s['n_revel']})")
    print(f"  -- on AM subset:     DISCERN {s['auc_discern_on_am_subset']} vs AlphaMissense {s['auc_alphamissense']}  (n={s['n_am']})")
    if "auc_intervar" in s:
        print(f"  AUC InterVar (ordinal): {s['auc_intervar']}  (n={s['n_intervar']})")


if __name__ == "__main__":
    main()
