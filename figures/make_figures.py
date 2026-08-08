"""Render every manuscript figure.

Each figure lives in its own module under ``figures/scripts/`` so it can be edited and re-rendered
on its own; this is only the driver. Geometry, palette and resolution come from
``figures/scripts/style.py``, which encodes the journal's requirements once.

Supplementary figures are numbered by first citation in the main text, which is why the module
names and the output stems do not run in step: the worked example is cited first and is Figure S1.

Run:  python -m figures.make_figures [--out DIR]
"""
from __future__ import annotations

import argparse
import os

from figures.scripts import (
    fig1_framework,
    fig2_discrimination_calibration,
    fig3_per_criterion_kappa,
    fig4_intrinsic_ceiling,
    fig5_clinvar_circularity,
    fig6_diagnosis_baselines,
    figS1_gene_term_sweep,
    figS2_safety_matrix,
    figS3_diagnosis_calibration,
    figS4_champ_chbmp,
    figS5_worked_example,
    graphical_abstract,
)

HERE = os.path.dirname(__file__)

FIGURES = [
    fig1_framework, fig2_discrimination_calibration, fig3_per_criterion_kappa,
    fig4_intrinsic_ceiling, fig5_clinvar_circularity, fig6_diagnosis_baselines,
    figS1_gene_term_sweep, figS2_safety_matrix, figS3_diagnosis_calibration,
    figS4_champ_chbmp, figS5_worked_example, graphical_abstract,
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "out"))
    args = ap.parse_args()
    for mod in FIGURES:
        paths = mod.build(args.out)
        print(f"  {mod.__name__.split('.')[-1]:34} -> {os.path.basename(paths[0])}")
    print(f"\nwrote {len(FIGURES)} figures to {args.out}")


if __name__ == "__main__":
    main()
