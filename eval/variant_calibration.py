"""Variant-score calibration vs ClinVar labels (v3.1 Track A3; Gate G11 / H5 variant half).

Maps DISCERN's variant-intrinsic Tavtigian points -> calibrated P(pathogenic) and measures
calibration against ClinVar's binary label (P/LP=1, B/LB=0; VUS/conflicting dropped). Reports,
on a held-out split, the uncalibrated baseline vs Platt (logistic) vs isotonic (PAV):
Brier score, ECE (10-bin), and AUC. Dependency-light (pure stdlib) so it runs on the VM.

Run:
    python3 -m eval.variant_calibration data/raw/erepo/erepo_classifications.tab \
        data/raw/clinvar/variant_summary.txt.gz
"""
from __future__ import annotations

import csv
import gzip
import math
import random
import sys

from rules.acmg_codes import code_points
from rules.vcep.partition import owner


def _intrinsic_points(codes: list[str]) -> float:
    pts = 0.0
    for c in codes:
        p, is_ba1 = code_points(c)
        if is_ba1:
            return -8.0
        if owner(c) == "variant_intrinsic":
            pts += p
    return pts


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


def _sigmoid(z):
    if z < -60:
        return 0.0
    if z > 60:
        return 1.0
    return 1.0 / (1.0 + math.exp(-z))


def _platt_fit(xs, ys, iters=2000, lr=0.02):
    """Logistic P=sigmoid(a*x+b) by gradient descent on scaled points."""
    mu = sum(xs) / len(xs)
    sd = (sum((x - mu) ** 2 for x in xs) / len(xs)) ** 0.5 or 1.0
    zs = [(x - mu) / sd for x in xs]
    a = b = 0.0
    n = len(zs)
    for _ in range(iters):
        ga = gb = 0.0
        for z, y in zip(zs, ys, strict=False):
            e = _sigmoid(a * z + b) - y
            ga += e * z
            gb += e
        a -= lr * ga / n
        b -= lr * gb / n
    return lambda x: _sigmoid(a * (x - mu) / sd + b)


def _isotonic_fit(xs, ys):
    """Pool-adjacent-violators isotonic regression; returns predictor f(x) (non-decreasing)."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    sx = [xs[i] for i in order]
    sy = [float(ys[i]) for i in order]
    # PAV
    val = sy[:]
    wt = [1.0] * len(sy)
    bx = sx[:]
    blocks = [[val[k], wt[k], bx[k]] for k in range(len(sy))]
    merged = True
    while merged:
        merged = False
        out = []
        for blk in blocks:
            while out and out[-1][0] > blk[0]:
                last = out.pop()
                tw = last[1] + blk[1]
                blk = [(last[0] * last[1] + blk[0] * blk[1]) / tw, tw, max(last[2], blk[2])]
                merged = True
            out.append(blk)
        blocks = out
    thr = [b[2] for b in blocks]
    vals = [b[0] for b in blocks]

    def f(x):
        lo, hi = 0, len(thr)
        while lo < hi:
            mid = (lo + hi) // 2
            if thr[mid] < x:
                lo = mid + 1
            else:
                hi = mid
        return vals[min(lo, len(vals) - 1)]
    return f


def _brier(probs, ys):
    return sum((p - y) ** 2 for p, y in zip(probs, ys, strict=False)) / len(ys)


def _ece(probs, ys, bins=10):
    tot = len(ys)
    e = 0.0
    for b in range(bins):
        lo, hi = b / bins, (b + 1) / bins
        idx = [i for i, p in enumerate(probs) if (lo <= p < hi or (b == bins - 1 and p == 1.0))]
        if not idx:
            continue
        conf = sum(probs[i] for i in idx) / len(idx)
        acc = sum(ys[i] for i in idx) / len(idx)
        e += len(idx) / tot * abs(conf - acc)
    return e


def _auc(probs, ys):
    pos = [p for p, y in zip(probs, ys, strict=False) if y == 1]
    neg = [p for p, y in zip(probs, ys, strict=False) if y == 0]
    if not pos or not neg:
        return None
    wins = sum((p > q) + 0.5 * (p == q) for p in pos for q in neg)
    return wins / (len(pos) * len(neg))


def run(erepo_path, clinvar_path, seed=0):
    want = {}
    with open(erepo_path, encoding="utf-8", errors="replace") as fh:
        r = csv.reader(fh, delimiter="\t")
        h = next(r)
        vi, ci = h.index("ClinVar Variation Id"), h.index("Applied Evidence Codes (Met)")
        for row in r:
            if len(row) <= max(vi, ci):
                continue
            vid = row[vi].strip()
            codes = [c.strip() for c in row[ci].split(",") if c.strip()]
            if vid and codes:
                want[vid] = _intrinsic_points(codes)

    data = []
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
            if vid not in want or vid in seen:
                continue
            lab = _clinvar_label(row[cls_i])
            if lab is None:
                continue
            seen.add(vid)
            data.append((want[vid], lab))

    rng = random.Random(seed)
    rng.shuffle(data)
    half = len(data) // 2
    cal, test = data[:half], data[half:]
    cx, cy = [d[0] for d in cal], [d[1] for d in cal]
    tx, ty = [d[0] for d in test], [d[1] for d in test]

    def base(x):
        return _sigmoid(x / 4.0)              # naive Tavtigian-points sigmoid baseline

    platt = _platt_fit(cx, cy)
    iso = _isotonic_fit(cx, cy)
    out = {"n_labeled": len(data), "n_pathogenic": sum(d[1] for d in data),
           "n_test": len(test)}
    for name, f in (("uncalibrated", base), ("platt", platt), ("isotonic", iso)):
        probs = [f(x) for x in tx]
        out[name] = {"brier": round(_brier(probs, ty), 4), "ece": round(_ece(probs, ty), 4),
                     "auc": round(_auc(probs, ty), 4) if _auc(probs, ty) else None}
    return out


def main():
    s = run(sys.argv[1], sys.argv[2])
    print(f"Labeled (P/B) variants: {s['n_labeled']} ({s['n_pathogenic']} pathogenic); test n={s['n_test']}")
    for k in ("uncalibrated", "platt", "isotonic"):
        m = s[k]
        print(f"  {k:13} Brier={m['brier']:.4f}  ECE={m['ece']:.4f}  AUC={m['auc']}")


if __name__ == "__main__":
    main()
