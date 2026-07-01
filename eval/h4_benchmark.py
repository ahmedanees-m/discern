"""Novel-variant scoring benchmark vs ClinVar (v3.1 Track A2 / H4).

Compares, on real ClinVar-labelled spec-gene variants, the AUC of:
  (a) DISCERN gene-specific intrinsic score (CSpec freq thresholds + REVEL at the spec cut-off),
  (b) a GENERIC-ACMG baseline (the InterVar paradigm: generic ACMG freq/REVEL cut-offs, no
      gene-specific VCEP thresholds), and
  (c) REVEL alone.
Per-variant gnomAD AF and REVEL are the values the curators cite in the eRepo narratives
(no predictor download needed). Labels: ClinVar P/LP=1, B/LB=0 (VUS/conflicting dropped).

SCOPE: this isolates the value of the GENE-SPECIFIC thresholds over generic ACMG (which
is what InterVar applies) on the freq+REVEL evidence available in-text. The LITERAL InterVar
tool requires an ANNOVAR + humandb install (academic-registration-gated) and is a separate
follow-up; the generic-ACMG baseline here is the same paradigm, transparently reimplemented.

Run:  python3 -m eval.h4_benchmark data/raw/erepo/erepo_classifications.tab \
        data/raw/clinvar/variant_summary.txt.gz
"""
from __future__ import annotations

import csv
import gzip
import re
import sys

from eval.gnomad_freq_check import SPEC_GENES, _observed_af
from rules.variant_scoring import Annotations, score_variant

_REVEL = re.compile(r"REVEL[^.]{0,30}?(\d\.\d+)", re.I)


def _parse_revel(narr: str):
    m = _REVEL.search(narr)
    return float(m.group(1)) if m else None


def _clinvar_label(txt: str):
    t = txt.lower()
    if "conflict" in t or "uncertain" in t:
        return None
    if "pathogenic" in t and "benign" not in t:
        return 1
    if "benign" in t and "pathogenic" not in t:
        return 0
    return None


def _col(header, *cands):
    for c in cands:
        if c in header:
            return header.index(c)
    raise KeyError(cands)


def _generic_acmg_points(af, revel):
    """The InterVar paradigm: generic ACMG freq + REVEL cut-offs, no gene-specific thresholds."""
    pts = 0.0
    if af is not None:
        if af >= 0.05:
            return -8.0                       # generic BA1 (stand-alone benign)
        if af >= 0.005:
            pts += -4.0                       # generic BS1
        elif af < 0.0001:
            pts += 2.0                        # generic PM2 at Moderate (the non-VCEP default)
    if revel is not None:
        if revel >= 0.7:
            pts += 1.0                        # generic PP3
        elif revel <= 0.15:
            pts += -1.0                       # generic BP4
    return pts


def _auc(scores, ys):
    pos = [s for s, y in zip(scores, ys, strict=False) if y == 1]
    neg = [s for s, y in zip(scores, ys, strict=False) if y == 0]
    if not pos or not neg:
        return None
    wins = sum((p > q) + 0.5 * (p == q) for p in pos for q in neg)
    return round(wins / (len(pos) * len(neg)), 4)


def run(erepo_path, clinvar_path):
    feats = {}
    with open(erepo_path, encoding="utf-8", errors="replace") as fh:
        r = csv.reader(fh, delimiter="\t")
        h = next(r)
        gi, vi, si = (h.index("HGNC Gene Symbol"), h.index("ClinVar Variation Id"),
                      h.index("Summary of interpretation"))
        for row in r:
            if len(row) <= max(gi, vi, si):
                continue
            gene, vid = row[gi].strip(), row[vi].strip()
            if gene not in SPEC_GENES or not vid:
                continue
            feats[vid] = (gene, _observed_af(row[si]), _parse_revel(row[si]))

    rows = []
    opener = gzip.open if clinvar_path.endswith(".gz") else open
    seen = set()
    with opener(clinvar_path, "rt", encoding="utf-8", errors="replace") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = [c.lstrip("#") for c in next(r)]
        vid_i = _col(header, "VariationID")
        cls_i = _col(header, "Germline classification", "ClinicalSignificance")
        for row in r:
            if len(row) <= max(vid_i, cls_i):
                continue
            vid = row[vid_i].strip()
            if vid not in feats or vid in seen:
                continue
            lab = _clinvar_label(row[cls_i])
            if lab is None:
                continue
            seen.add(vid)
            gene, af, revel = feats[vid]
            sv = score_variant(gene, vid, Annotations(af=af, revel=revel))
            rows.append((sv.points, _generic_acmg_points(af, revel), revel, lab))

    ys = [r[3] for r in rows]
    n_rev = [r for r in rows if r[2] is not None]
    return {
        "n_labeled": len(rows), "n_pathogenic": sum(ys), "n_with_revel": len(n_rev),
        "auc_discern_gene_specific": _auc([r[0] for r in rows], ys),
        "auc_generic_acmg_intervar_paradigm": _auc([r[1] for r in rows], ys),
        "auc_revel_alone": _auc([r[2] for r in n_rev], [r[3] for r in n_rev]),
    }


def main():
    s = run(sys.argv[1], sys.argv[2])
    print(f"Spec-gene ClinVar-labelled variants: {s['n_labeled']} "
          f"({s['n_pathogenic']} pathogenic; REVEL in {s['n_with_revel']})")
    print(f"  AUC DISCERN gene-specific:            {s['auc_discern_gene_specific']}")
    print(f"  AUC generic-ACMG (InterVar paradigm): {s['auc_generic_acmg_intervar_paradigm']}")
    print(f"  AUC REVEL alone:                      {s['auc_revel_alone']}")


if __name__ == "__main__":
    main()
