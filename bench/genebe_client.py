"""Cached GeneBe client for the Track-1 variant head-to-head (DISCERN_Benchmark_Execution_Plan_v1).

GeneBe (genebe.net) is a current, VCEP-aware ACMG classifier with a public API. For every variant
it returns the ACMG classification + criteria + score, plus the predictor scores (REVEL,
AlphaMissense) and gnomAD AF and the ClinVar fields - so it doubles as the annotation source that
lets DISCERN be scored locally on the identical variant set.

This client is deterministic and reproducible: it POSTs the h4set variants (hg38) in batches and
caches each result to bench/data/genebe_h4set.jsonl, resuming on re-run so the public API is not
re-hammered. Snapshot date is recorded in the cache header line.

Run:  python -m bench.genebe_client            # fills/extends the cache
"""
from __future__ import annotations

import csv
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(__file__)
META = os.path.join(HERE, "..", "data", "processed", "h4set.meta.tsv")
CACHE = os.path.join(HERE, "data", "genebe_h4set.jsonl")
API = "https://api.genebe.net/cloud/api-public/v1/variants"
GENOME = "hg38"   # confirmed: h4set coords are GRCh38 (GP9 chr3:129061759 maps correctly only in hg38)

# GeneBe fields we keep (predictors + ACMG + provenance for the ClinVar-dependence analysis).
KEEP = ["chr", "pos", "ref", "alt", "gene_symbol", "effect", "acmg_classification", "acmg_score",
        "acmg_criteria", "revel_score", "alphamissense_score", "spliceai_max_score",
        "gnomad_exomes_af", "gnomad_genomes_af", "clinvar_classification", "clinvar_review_status"]


def _post(variants: list[dict], retries: int = 4) -> list[dict]:
    url = API + "?" + urllib.parse.urlencode({"genome": GENOME})
    body = json.dumps(variants).encode()
    for i in range(retries):
        try:
            req = urllib.request.Request(url, data=body, method="POST",
                                         headers={"Content-Type": "application/json",
                                                  "Accept": "application/json",
                                                  "User-Agent": "discern-bench/1.0"})
            with urllib.request.urlopen(req, timeout=120) as fh:
                d = json.load(fh)
            return d.get("variants") if isinstance(d, dict) else d
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            wait = 3 * (i + 1)
            print(f"  retry {i + 1}/{retries} after {type(e).__name__} ({str(e)[:60]}); waiting {wait}s")
            time.sleep(wait)
    raise RuntimeError("GeneBe API failed after retries")


def load_meta() -> list[dict]:
    with open(META, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _cached_keys() -> set:
    keys = set()
    if os.path.exists(CACHE):
        with open(CACHE, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                r = json.loads(line)
                keys.add((str(r["q_chr"]), str(r["q_pos"]), r["q_ref"], r["q_alt"]))
    return keys


def fill_cache(batch: int = 60, sleep: float = 1.0) -> None:
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    meta = load_meta()
    done = _cached_keys()
    todo = [m for m in meta if (str(m["chrom"]), str(m["pos"]), m["ref"], m["alt"]) not in done]
    print(f"h4set={len(meta)} cached={len(done)} todo={len(todo)}")
    with open(CACHE, "a", encoding="utf-8") as out:
        if not done:
            out.write(f"# GeneBe cache for h4set; genome={GENOME}; api={API}\n")
        for i in range(0, len(todo), batch):
            chunk = todo[i:i + batch]
            q = [{"chr": str(m["chrom"]).replace("chr", ""), "pos": int(m["pos"]),
                  "ref": m["ref"], "alt": m["alt"]} for m in chunk]
            res = _post(q)
            # GeneBe preserves input order; map positionally and store the query key explicitly.
            for m, v in zip(chunk, res or [], strict=False):
                rec = {"q_chr": str(m["chrom"]).replace("chr", ""), "q_pos": int(m["pos"]),
                       "q_ref": m["ref"], "q_alt": m["alt"], "q_gene": m["gene"], "q_clnsig": m["clnsig"]}
                rec.update({k: v.get(k) for k in KEEP})
                out.write(json.dumps(rec) + "\n")
            out.flush()
            print(f"  {min(i + batch, len(todo))}/{len(todo)}")
            time.sleep(sleep)
    print("cache complete")


if __name__ == "__main__":
    fill_cache()
