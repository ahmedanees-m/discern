"""Annotate the eRepo bleeding SNV set via GeneBe, using the API and a local cache only.

Reuses the cached GeneBe client from Track 1. Resumable; caches to bench/data/genebe_erepo.jsonl.
Run:  python -m bench.genebe_erepo
"""
from __future__ import annotations

import csv
import json
import os
import time

from bench.genebe_client import KEEP, _post

HERE = os.path.dirname(__file__)
SRC = os.path.join(HERE, "data", "erepo_bleeding.tsv")
CACHE = os.path.join(HERE, "data", "genebe_erepo.jsonl")


def _cached():
    keys = set()
    if os.path.exists(CACHE):
        for line in open(CACHE, encoding="utf-8"):
            line = line.strip()
            if line and not line.startswith("#"):
                r = json.loads(line)
                keys.add((r["q_chr"], str(r["q_pos"]), r["q_ref"], r["q_alt"]))
    return keys


def fill(batch=60, sleep=1.0):
    rows = list(csv.DictReader(open(SRC, encoding="utf-8"), delimiter="\t"))
    done = _cached()
    todo = [r for r in rows if (r["chrom"], r["pos"], r["ref"], r["alt"]) not in done]
    print(f"erepo SNVs={len(rows)} cached={len(done)} todo={len(todo)}")
    with open(CACHE, "a", encoding="utf-8") as out:
        if not done:
            out.write("# GeneBe cache for eRepo bleeding SNVs; genome=hg38\n")
        for i in range(0, len(todo), batch):
            chunk = todo[i:i + batch]
            q = [{"chr": r["chrom"], "pos": int(r["pos"]), "ref": r["ref"], "alt": r["alt"]} for r in chunk]
            res = _post(q)
            for r, v in zip(chunk, res or [], strict=False):
                rec = {"q_chr": r["chrom"], "q_pos": int(r["pos"]), "q_ref": r["ref"], "q_alt": r["alt"],
                       "q_gene": r["gene"], "q_assertion": r["assertion"],
                       "q_approval_date": r["approval_date"], "q_codes_met": r["codes_met"]}
                rec.update({k: v.get(k) for k in KEEP})
                out.write(json.dumps(rec) + "\n")
            out.flush()
            print(f"  {min(i + batch, len(todo))}/{len(todo)}")
            time.sleep(sleep)
    print("eRepo GeneBe cache complete")


if __name__ == "__main__":
    fill()
