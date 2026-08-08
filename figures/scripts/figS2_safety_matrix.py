"""Figure S3: the treatment-divergence safety interlock, by scenario.

Rendered from the committed benchmark JSON, on the shared journal style in style.py.

Rows are labelled from the vignette file and the disease ontology rather than from the raw
scenario ids: the ids carry underscores and internal abbreviations that mean nothing to a reader.
Both columns are drawn as cells, so the matrix reads as a matrix; the earlier version left the
"stays silent" column blank and the panel looked like one solid block.
"""
from __future__ import annotations

import json
import os

import matplotlib.pyplot as plt
import yaml

from figures.scripts import style as st

# scripts/ -> figures/ -> repository root
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BENCH = os.path.join(ROOT, "bench")
EVAL = os.path.join(ROOT, "eval")

# Presentation only: the underscored treatment keys spelled the way a clinician writes them.
TX = {"ddavp": "desmopressin",
      "splenectomy": "splenectomy",
      "affected_related_donor_transplant": "affected-relative donor transplant",
      "platelet_transfusion": "platelet transfusion",
      "recombinant_fviii_monotherapy": "recombinant FVIII monotherapy"}


def _load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _scenario_labels():
    """scenario id -> (disease name, gene), read from the vignettes and the disease ontology."""
    from diseases.ontology import all_clusters

    names = {d.id: d.name.replace(" (pseudo)", "")
             for c in all_clusters().values() for d in c.diseases}
    with open(os.path.join(EVAL, "cases", "reader_vignettes.yaml"), encoding="utf-8") as fh:
        vignettes = yaml.safe_load(fh)["vignettes"]
    return {v["id"]: (names.get(v.get("true_dx"), v.get("true_dx", "")), v.get("gene", ""))
            for v in vignettes}


def _row_label(labels, row):
    """Disease, gene and planned management, without repeating a gene the name already carries."""
    disease, gene = labels.get(row["id"], ("", ""))
    head = disease if (not gene or gene in disease) else f"{disease} ({gene})"
    tx = TX.get(row["harmful_tx"], row["harmful_tx"].replace("_", " "))
    return f"{head} - planned: {tx}"


def build(outdir):
    st.apply()
    m = _load(os.path.join(BENCH, "track3_metrics.json"))["safety_interlock"]
    rows = m["rows"]
    labels = _scenario_labels()

    fig, ax = plt.subplots(figsize=(st.FULL_W, 0.46 * len(rows) + 1.25))
    fig.subplots_adjust(left=0.315, right=0.995, top=0.80, bottom=0.185)

    cols = [("planned management is contraindicated\n(the interlock must fire)", "hardstop_on_real",
             True),
            ("planned management is harmless\n(the interlock must stay silent)",
             "hardstop_on_harmless", False)]
    n = len(rows)
    for i, r in enumerate(rows):
        y = n - 1 - i
        for j, (_title, key, must_fire) in enumerate(cols):
            fired = bool(r[key])
            correct = fired == must_fire
            face = st.VERMILION if fired else st.GREEN
            ax.add_patch(plt.Rectangle((j - 0.46, y - 0.40), 0.92, 0.80,
                                       facecolor=face, alpha=0.85 if fired else 0.30,
                                       edgecolor="white", lw=1.2, zorder=2))
            ax.text(j, y, "HARD STOP" if fired else "silent", ha="center", va="center",
                    fontsize=st.TINY, zorder=3, fontweight="bold" if fired else "normal",
                    color="white" if fired else st.GREEN)
            # A wrong cell would be marked; none is, and the figure should show that it would be.
            if not correct:
                ax.text(j + 0.40, y + 0.30, "x", ha="right", va="top", fontsize=st.SMALL,
                        color="black", zorder=4)

    ax.set_xlim(-0.55, len(cols) - 0.45)
    ax.set_ylim(-0.55, n - 0.45)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c[0] for c in cols], fontsize=st.TINY)
    ax.set_yticks(range(n))
    ax.set_yticklabels([st.wrap(_row_label(labels, r), 34) for r in rows][::-1],
                       fontsize=st.TINY, linespacing=1.35)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(f"Safety interlock: sensitivity {m['hardstop_sensitivity']:.0%}, "
                 f"specificity {m['hardstop_specificity']:.0%}, judged gene-blind",
                 loc="left", fontweight="bold", pad=22)
    return st.save(fig, outdir, "figS2_safety_matrix")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ROOT)
    print(build(os.path.join(ROOT, "figures", "out")))
