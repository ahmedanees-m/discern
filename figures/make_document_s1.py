"""Compile Document S1: all supplemental information as a single PDF.

Cell Press requires supplemental information to arrive as one compiled PDF rather than as loose
files, aside from large data tables which may be separate spreadsheets. This assembles the contents
page, the supplemental methods, the five supplemental figures with their legends, and the smaller
supplemental tables, in that order.

Text pages are typeset with matplotlib and the figure pages are the original vector PDFs, merged
with pypdf, so the figures keep their vector quality rather than being re-rastered.

Run:  python -m figures.make_document_s1 [--figures DIR] [--tables DIR] [--out FILE]
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from pypdf import PdfReader, PdfWriter  # noqa: E402

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
DOCS = os.path.join(ROOT, "docs")

PAGE = (8.27, 11.69)          # A4 portrait, inches
MARGIN_L, MARGIN_T = 0.09, 0.94
BODY, HEAD, SUB = 8.2, 13.5, 10.0
LEADING = 0.0165

# Supplemental tables small enough to typeset inside the document. The large ones ship as Excel.
# Ordered and numbered by first citation in the main text; the four large data tables are cited
# first and ship as Excel, so the tables typeset here are S5 to S9.
INLINE_TABLES = ["tableS1_per_criterion_kappa", "tableS2_ceiling_attribution",
                 "tableS4_safety_scenarios", "tableS8_lirical_arms",
                 "tableS6_champ_chbmp_recall"]
FIGURE_ORDER = ["figS5_worked_example", "figS1_gene_term_sensitivity", "figS2_safety_matrix",
                "figS3_diagnosis_calibration", "figS4_champ_chbmp_recall"]


class Pages:
    """Accumulates text pages; each `flush` starts a new sheet."""

    def __init__(self, tmpdir):
        self.tmp = tmpdir
        self.paths = []
        self._new()

    def _new(self):
        self.fig = plt.figure(figsize=PAGE)
        self.y = MARGIN_T

    def _room(self, need):
        if self.y - need < 0.06:
            self.flush()

    def text(self, s, size=BODY, weight="normal", style="normal", color="black", indent=0.0,
             gap=0.0):
        self.y -= gap
        self._room(LEADING * (size / BODY))
        self.fig.text(MARGIN_L + indent, self.y, s, fontsize=size, fontweight=weight,
                      style=style, color=color, va="top", ha="left", wrap=False)
        self.y -= LEADING * (size / BODY) * 1.25

    def para(self, s, size=BODY, width=105, indent=0.0, gap=0.010, **kw):
        import textwrap
        self.y -= gap
        for line in textwrap.wrap(s, width=width) or [""]:
            self.text(line, size=size, indent=indent, **kw)

    def heading(self, s, size=HEAD, gap=0.030):
        self.y -= gap
        self._room(0.10)
        self.text(s, size=size, weight="bold")
        self.y -= 0.006

    def rule(self, gap=0.008):
        self.y -= gap
        self.fig.add_artist(plt.Line2D([MARGIN_L, 0.93], [self.y, self.y], color="#999999",
                                       lw=0.6, transform=self.fig.transFigure))
        self.y -= 0.010

    def flush(self):
        if self.y >= MARGIN_T - 1e-9:          # nothing written to this sheet
            plt.close(self.fig)
            self._new()
            return
        p = os.path.join(self.tmp, f"page_{len(self.paths):03d}.pdf")
        self.fig.savefig(p, format="pdf")
        plt.close(self.fig)
        self.paths.append(p)
        self._new()

    def done(self):
        self.flush()
        return self.paths


def _md_lines(path):
    """Very small markdown reader: headings, bullets, bold-lead paragraphs, plain paragraphs."""
    out = []
    for raw in open(path, encoding="utf-8").read().split("\n"):
        line = raw.rstrip()
        if line.startswith("# "):
            out.append(("h1", line[2:]))
        elif line.startswith("## "):
            out.append(("h2", line[3:]))
        elif line.startswith("### "):
            out.append(("h3", line[4:]))
        elif line.strip() == "---":
            out.append(("rule", ""))
        elif line.startswith("- ") or line.startswith("* "):
            out.append(("li", line[2:]))
        elif line.startswith("    ") and line.strip():
            out.append(("pre", line.strip()))
        elif not line.strip():
            out.append(("blank", ""))
        else:
            out.append(("p", line))
    return out


def _clean(s):
    s = re.sub(r"\*\*(.+?)\*\*", r"\1", s)
    s = re.sub(r"`(.+?)`", r"\1", s)
    s = re.sub(r"\*(.+?)\*", r"\1", s)
    s = re.sub(r"\[(.+?)\]\(.+?\)", r"\1", s)
    return s


def render_markdown(pg: Pages, path, skip_h1=False):
    buf = []

    def flush_p():
        if buf:
            pg.para(" ".join(buf))
            buf.clear()

    for kind, txt in _md_lines(path):
        if kind == "blank":
            flush_p()
        elif kind == "h1":
            flush_p()
            if not skip_h1:
                pg.heading(_clean(txt))
        elif kind == "h2":
            flush_p()
            pg.heading(_clean(txt), size=SUB, gap=0.026)
        elif kind == "h3":
            flush_p()
            pg.heading(_clean(txt), size=BODY + 0.6, gap=0.020)
        elif kind == "rule":
            flush_p()
            pg.rule()
        elif kind == "li":
            flush_p()
            pg.para("- " + _clean(txt), indent=0.012, gap=0.003)
        elif kind == "pre":
            flush_p()
            pg.text("    " + _clean(txt), size=BODY - 0.8, color="#333333")
        else:
            buf.append(_clean(txt))
    flush_p()


def render_table(pg: Pages, csv_path, title):
    with open(csv_path, encoding="utf-8") as fh:
        rows = [r for r in csv.reader(fh) if r]
    header = rows[0]
    body = [r for r in rows[1:] if any(c.strip() for c in r) and not r[0].startswith("NOTE:")]
    notes = [r[0] for r in rows[1:] if r and r[0].startswith("NOTE:")]

    pg.heading(title, size=SUB, gap=0.028)
    ncol = len(header)
    widths = [max(len(str(header[j])), *(len(str(r[j])) for r in body if len(r) > j))
              if body else len(str(header[j])) for j in range(ncol)]
    widths = [min(w, 34) for w in widths]           # cap runaway free-text columns
    total = sum(widths) or 1
    avail = 116
    widths = [max(9, int(avail * w / total)) for w in widths]

    def fmt(cells):
        return "  ".join(str(c)[:widths[j]].ljust(widths[j]) for j, c in enumerate(cells[:ncol]))

    pg.text(fmt(header), size=BODY - 1.6, weight="bold")
    pg.text("-" * min(sum(widths) + 2 * ncol, 120), size=BODY - 1.6, color="#999999")
    for r in body:
        pg.text(fmt(r + [""] * (ncol - len(r))), size=BODY - 1.6)
    for n in notes:
        pg.para(n, size=BODY - 1.4, style="italic", gap=0.008)


def build(fig_dir, tab_dir, out_path):
    tmp = tempfile.mkdtemp(prefix="docs1-")
    pg = Pages(tmp)

    # --- title and contents -----------------------------------------------------------
    pg.text("Document S1", size=18, weight="bold")
    pg.text("Supplemental Information", size=13)
    pg.y -= 0.014
    pg.para("A partitioned, calibrated framework for variant interpretation and differential "
            "diagnosis in inherited bleeding and platelet disorders", size=9.5, style="italic")
    pg.rule(gap=0.016)
    pg.heading("Contents", size=SUB)
    contents = [
        "Supplemental Methods",
        "Figure S1. A worked example: the safety interlock on a single case",
        "Figure S2. Sensitivity of the diagnosis arm to the strength of the gene term",
        "Figure S3. Treatment-divergence safety interlock, by scenario",
        "Figure S4. Diagnosis-posterior calibration on the curated benchmark",
        "Figure S5. Independent sensitivity on the CDC hemophilia mutation projects",
        "Table S5. Per-criterion agreement with expert-panel applications",
        "Table S6. Attribution of the intrinsic-evidence ceiling",
        "Table S7. Treatment-divergence scenarios",
        "Table S8. LIRICAL comparison arms",
        "Table S9. CHAMP and CHBMP recall by gene",
    ]
    for c in contents:
        pg.text(c, indent=0.012)
    pg.y -= 0.012
    pg.para("Tables S1, S2, S3 and S4 are supplied as separate Excel files, being large data "
            "tables. All supplemental figures and tables are also archived, with the scripts that "
            "generate them, in the data record.", size=BODY - 0.4, style="italic")
    pg.flush()

    # --- supplemental methods ---------------------------------------------------------
    render_markdown(pg, os.path.join(DOCS, "DISCERN_Supplemental_Methods.md"))
    pg.flush()

    # --- legends ----------------------------------------------------------------------
    pg.heading("Supplemental figure and table legends", size=HEAD)
    render_markdown(pg, os.path.join(DOCS, "DISCERN_Supplemental_Legends.md"), skip_h1=True)
    pg.flush()

    text_pages = pg.done()

    # --- merge: text, then the vector figures, then the inline tables ------------------
    writer = PdfWriter()
    for p in text_pages:
        for page in PdfReader(p).pages:
            writer.add_page(page)
    for stem in FIGURE_ORDER:
        f = os.path.join(fig_dir, f"{stem}.pdf")
        if os.path.exists(f):
            for page in PdfReader(f).pages:
                writer.add_page(page)
        else:
            print(f"  ! missing figure {stem}.pdf")

    pg2 = Pages(tmp)
    pg2.heading("Supplemental tables", size=HEAD)
    titles = {
        "tableS1_per_criterion_kappa": "Table S5. Per-criterion agreement with expert-panel applications",
        "tableS2_ceiling_attribution": "Table S6. Attribution of the intrinsic-evidence ceiling",
        "tableS4_safety_scenarios": "Table S7. Treatment-divergence scenarios",
        "tableS6_champ_chbmp_recall": "Table S9. CHAMP and CHBMP recall by gene",
        "tableS8_lirical_arms": "Table S8. LIRICAL comparison arms",
    }
    for stem in INLINE_TABLES:
        c = os.path.join(tab_dir, f"{stem}.csv")
        if os.path.exists(c):
            render_table(pg2, c, titles[stem])
            pg2.flush()
        else:
            print(f"  ! missing table {stem}.csv")
    for p in pg2.done():
        for page in PdfReader(p).pages:
            writer.add_page(page)

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "wb") as fh:
        writer.write(fh)
    return out_path, len(writer.pages)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figures", default=os.path.join(HERE, "out"))
    ap.add_argument("--tables", default=os.path.join(HERE, "tables"))
    ap.add_argument("--out", default=os.path.join(HERE, "Document_S1_Supplemental_Information.pdf"))
    args = ap.parse_args()
    path, n = build(args.figures, args.tables, args.out)
    print(f"wrote {path} ({n} pages, {os.path.getsize(path) / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
