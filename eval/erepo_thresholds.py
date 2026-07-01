"""Triangulate VCEP frequency/computational thresholds from real ClinGen ERepo narratives.

Refined: only counts CURATOR-STATED thresholds (phrases like "BA1 cut-off of >=0.000333",
"BS1 criteria of MAF > 0.000167", "threshold of >0.0024"), not the variant's own observed
AF. This removes the noise from proximity matching. Confirms applied PM2 strength too.

Pure stdlib; run on the VM:
    python3 eval/erepo_thresholds.py data/raw/erepo/erepo_classifications.tab
"""
import csv
import re
import sys
from collections import Counter, defaultdict

GROUPS = {
    "GT (ITGA2B/ITGB3)": {"ITGA2B", "ITGB3"},
    "F8": {"F8"}, "F9": {"F9"}, "VWF": {"VWF"}, "RUNX1": {"RUNX1"},
    "BSS (GP1BA/GP1BB/GP9)": {"GP1BA", "GP1BB", "GP9"},
}
GENE_TO_GROUP = {g: name for name, gs in GROUPS.items() for g in gs}

NUM = r"(\d+\.\d+(?:[eE]-?\d+)?|\d+[eE]-?\d+)"
# curator-stated freq threshold: CODE ... (cut-off|criteria|threshold) ... (>|>=|<) NUM
THR = {c: re.compile(c + r"[^.]{0,60}?(?:cut[- ]?off|criteria|threshold)[^.]{0,25}?(>=|>/=|>|<|less than|greater than)?\s*" + NUM, re.I)
       for c in ("BA1", "BS1", "PM2")}
# computational thresholds stated as cut-offs
REVEL_THR = re.compile(r"REVEL[^.]{0,50}?(?:cut[- ]?off|threshold|>=|>|<|of)\s*(\d\.\d+)", re.I)
CADD_THR = re.compile(r"CADD[^.]{0,50}?(?:cut[- ]?off|threshold|>=|>|<|of)\s*(\d+(?:\.\d+)?)", re.I)
CODEMOD = re.compile(r"\b(PM2|BA1|BS1)(_[A-Za-z ]+)?")


def run(path):
    strengths = defaultdict(Counter)
    thr = defaultdict(lambda: defaultdict(Counter))
    revel = defaultdict(Counter)
    cadd = defaultdict(Counter)
    n = defaultdict(int)
    with open(path, encoding="utf-8", errors="replace") as fh:
        r = csv.reader(fh, delimiter="\t")
        h = next(r)
        gi, ci, si = h.index("HGNC Gene Symbol"), h.index("Applied Evidence Codes (Met)"), h.index("Summary of interpretation")
        for row in r:
            if len(row) <= max(gi, ci, si):
                continue
            grp = GENE_TO_GROUP.get(row[gi].strip())
            if not grp:
                continue
            n[grp] += 1
            for m in CODEMOD.finditer(row[ci]):
                if m.group(1) in {"PM2", "BA1", "BS1"}:
                    strengths[grp][m.group(0).strip()] += 1
            narr = row[si]
            for code, pat in THR.items():
                for m in pat.finditer(narr):
                    op = (m.group(1) or "").replace("/", "").strip()
                    thr[grp][code][f"{op}{m.group(2)}"] += 1
            for m in REVEL_THR.finditer(narr):
                revel[grp][m.group(1)] += 1
            for m in CADD_THR.finditer(narr):
                cadd[grp][m.group(1)] += 1
    for grp in GROUPS:
        if not n[grp]:
            continue
        print("=" * 72)
        print(f"{grp}: {n[grp]} variants | applied strengths: {dict(strengths[grp])}")
        for code in ("BA1", "BS1", "PM2"):
            if thr[grp][code]:
                print(f"  {code} curator-stated threshold (value:count): {thr[grp][code].most_common(6)}")
        if revel[grp]:
            print(f"  REVEL cut-off (value:count): {revel[grp].most_common(5)}")
        if cadd[grp]:
            print(f"  CADD cut-off (value:count): {cadd[grp].most_common(5)}")


if __name__ == "__main__":
    run(sys.argv[1])
