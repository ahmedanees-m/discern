"""Per-variant AlphaMissense lookup for the H4 variant set (v3.1 full H4).

Reads the variant meta (chrom,pos,ref,alt) and the precomputed AlphaMissense_hg38.tsv.gz
(DeepMind, CC BY 4.0), emitting the max am_pathogenicity per variant. One streaming pass.

Run on the VM:
    python3 -m eval.lookup_alphamissense /tmp/h4set.meta.tsv \
        ~/predictors/AlphaMissense_hg38.tsv.gz /tmp/am.tsv
"""
import gzip
import sys


def run(meta_path, am_path, out_path):
    want = {}     # (chrom_no_chr, pos, ref, alt) -> vid
    with open(meta_path, encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            vid, chrom, pos, ref, alt = line.rstrip("\n").split("\t")[:5]
            want[(chrom.replace("chr", ""), pos, ref, alt)] = vid
    best = {}
    with gzip.open(am_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            key = (f[0].replace("chr", ""), f[1], f[2], f[3])
            vid = want.get(key)
            if vid is None:
                continue
            try:
                am = float(f[8])
            except ValueError:
                continue
            if vid not in best or am > best[vid]:
                best[vid] = am
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("vid\tam_pathogenicity\n")
        for vid, am in best.items():
            out.write(f"{vid}\t{am}\n")
    print(f"variants: {len(want)}; AlphaMissense scored: {len(best)}")


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], sys.argv[3])
