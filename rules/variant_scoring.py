"""Novel-variant intrinsic scoring with the variant-dependent strength trees (v3.1 Track A2).

`score_variant(gene, hgvs, ann)` assembles the variant-INTRINSIC ACMG codes from raw inputs
against the gene's CSpec thresholds (the loader's per-gene af, PM2 strength, and `computational`
PP3/BP4 cut-offs), applies the variant-dependent PVS1 (Abou Tayoun et al. 2018 decision tree)
and PS4 (case-control OR / proband-ratio) strength trees, sums Tavtigian points, and returns the
band. Predictors (gnomAD AF, REVEL, Pangolin/SpliceAI, AlphaMissense) are INJECTABLE via the
Annotations object; missing predictors degrade gracefully to reduced confidence.

Routed codes (PP4/PS3/PM3/PP1) are deliberately NOT assembled here - they belong to the other
DISCERN factors (the per-code partition). The computational cut-offs live in each spec's
`computational` block and are shared with adapters/insilico.py and adapters/splice.py (one source).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from core.schemas import Strength
from rules.acmg_codes import code_points
from rules.point_engine import BANDS, Classification
from rules.vcep.loader import get_spec

# Fallback computational cut-offs when a spec omits the block (ACMG/SVI defaults).
_COMPUTATIONAL_DEFAULT = {"pp3_revel": 0.644, "pp3_splice": 0.5, "bp4_revel": 0.290, "bp4_splice": 0.1}

_STRENGTH_WORD = {Strength.PVS: "VeryStrong", Strength.PS: "Strong", Strength.PM: "Moderate",
                  Strength.PP: "Supporting"}

# PVS1 decision-tree consequence classes (Abou Tayoun et al., Hum Mutat 2018).
_NMD_TRUNCATING = {"nonsense", "stop_gained", "frameshift", "frameshift_deletion", "frameshift_insertion"}
_SPLICE = {"canonical_splice", "splice_donor", "splice_acceptor", "splicing"}
_START = {"start_lost", "start_loss", "initiator_codon", "startloss"}
_DELETION = {"deletion", "inframe_deletion", "exon_deletion", "nonframeshift_deletion"}
_LOF_ALL = _NMD_TRUNCATING | _SPLICE | _START | _DELETION


@dataclass
class Annotations:
    """Raw per-variant inputs; any field may be None (unknown) -> graceful degradation."""
    af: float | None = None                  # gnomAD popmax/grpmax filtering allele frequency
    revel: float | None = None
    splice: float | None = None              # Pangolin (preferred) or SpliceAI score
    alphamissense: float | None = None       # secondary feature only (not a primary code)
    consequence: str = ""
    # PVS1 decision-tree inputs (Tayoun 2018):
    nmd_predicted: bool | None = None        # truncation predicted to trigger NMD
    in_functional_domain: bool | None = None  # hits a clinically established functional domain
    removes_gt10pct: bool | None = None      # truncation removes >10% of the protein
    exon_in_biorelevant_transcript: bool | None = None  # False -> downgrade to Supporting
    # PS4 inputs:
    proband_count: int | None = None
    expected_proband: float | None = None
    odds_ratio: float | None = None          # case-control OR
    or_ci_low: float | None = None           # lower bound of the OR 95% CI


@dataclass
class ScoredVariant:
    gene: str
    hgvs: str
    spec_id: str
    covered: bool
    codes: list[str]
    points: float
    classification: Classification
    confidence: str
    drivers: dict = field(default_factory=dict)


def pvs1_strength(ann: Annotations) -> Strength | None:
    """ClinGen SVI PVS1 decision tree (Abou Tayoun et al. 2018), reduced to the available inputs.

    NMD-triggering truncation -> VeryStrong; NMD-escaping but removing a critical domain / >10%
    of the protein -> Strong, else Moderate; canonical-splice mirrors truncation; initiation-codon
    and in-frame deletions -> Strong (critical domain) else Moderate. A variant not in a
    biologically-relevant transcript is downgraded to Supporting.
    """
    c = ann.consequence
    if c not in _LOF_ALL:
        return None
    if ann.exon_in_biorelevant_transcript is False:
        return Strength.PP                       # not in a biologically-relevant transcript
    critical = bool(ann.removes_gt10pct or ann.in_functional_domain)
    if c in _NMD_TRUNCATING or c in _SPLICE:
        if ann.nmd_predicted is True:
            return Strength.PVS                  # undergoes NMD -> VeryStrong
        return Strength.PS if critical else Strength.PM   # escapes NMD
    if c in _START:
        return Strength.PM                       # initiation-codon variant
    if c in _DELETION:
        return Strength.PS if critical else Strength.PM
    return None


def ps4_strength(ann: Annotations) -> Strength | None:
    """PS4 by case-control odds ratio (CI lower bound > 1) when available, else proband-ratio /
    proband-count proxy (ClinGen SVI / CFD-VCEP ratio approach)."""
    if ann.odds_ratio is not None and ann.or_ci_low is not None:
        if ann.or_ci_low <= 1.0:
            return None                          # not statistically significant
        if ann.odds_ratio >= 5.0:
            return Strength.PS
        if ann.odds_ratio >= 3.0:
            return Strength.PM
        if ann.odds_ratio >= 2.0:
            return Strength.PP
        return None
    if ann.proband_count is None:
        return None
    if ann.expected_proband:
        ratio = ann.proband_count / max(ann.expected_proband, 1e-9)
        return (Strength.PS if ratio >= 10 else Strength.PM if ratio >= 4
                else Strength.PP if ratio >= 2 else None)
    n = ann.proband_count
    return Strength.PS if n >= 15 else Strength.PM if n >= 5 else Strength.PP if n >= 2 else None


def _code(base: str, strength: Strength) -> str:
    word = _STRENGTH_WORD.get(strength)
    return f"{base}_{word}" if word else base


def _classify(pts: float) -> Classification:
    for thr, cls in BANDS:
        if pts >= thr:
            return cls
    return Classification.B


def score_variant(gene: str, hgvs: str = "", ann: Annotations | None = None) -> ScoredVariant:
    ann = ann or Annotations()
    spec = get_spec(gene)
    codes: list[str] = []
    drivers: dict = {}
    missing: list[str] = []

    # --- frequency (BA1/BS1/PM2) against the gene-specific CSpec thresholds ---
    if ann.af is not None:
        ba1, bs1, pm2 = spec.af_threshold("ba1"), spec.af_threshold("bs1"), spec.af_threshold("pm2")
        if ann.af >= ba1:
            codes.append("BA1")
            drivers["BA1"] = f"AF {ann.af:g} >= {ba1:g}"
        elif ann.af >= bs1:
            codes.append("BS1")
            drivers["BS1"] = f"AF {ann.af:g} >= {bs1:g}"
        elif ann.af < pm2:
            pm2_strength = spec.strength_for("PM2")
            codes.append(_code("PM2", pm2_strength))
            drivers["PM2"] = f"AF {ann.af:g} < {pm2:g} ({pm2_strength.name})"
    else:
        missing.append("gnomAD AF")

    # --- computational PP3 / BP4 from the spec's own cut-offs (shared with the in-silico adapter) ---
    comp = spec.computational or _COMPUTATIONAL_DEFAULT
    pp3r, pp3s = comp.get("pp3_revel", 0.644), comp.get("pp3_splice", 0.5)
    bp4r, bp4s = comp.get("bp4_revel", 0.290), comp.get("bp4_splice", 0.1)
    pp3r_mod = comp.get("pp3_revel_moderate")
    if ann.revel is not None or ann.splice is not None:
        if pp3r_mod is not None and ann.revel is not None and ann.revel >= pp3r_mod:
            codes.append("PP3_Moderate")
            drivers["PP3"] = f"REVEL={ann.revel} >= {pp3r_mod} (Moderate)"
        elif (ann.revel is not None and ann.revel >= pp3r) or (ann.splice is not None and ann.splice >= pp3s):
            codes.append("PP3_Supporting")
            drivers["PP3"] = f"REVEL={ann.revel}, splice={ann.splice}"
        elif (ann.revel is not None and ann.revel <= bp4r) and (ann.splice is None or ann.splice <= bp4s):
            codes.append("BP4_Supporting")
            drivers["BP4"] = f"REVEL={ann.revel}, splice={ann.splice}"
    else:
        missing.append("REVEL/Pangolin")

    # --- variant-dependent strength trees ---
    pvs = pvs1_strength(ann)
    if pvs:
        codes.append(_code("PVS1", pvs))
        drivers["PVS1"] = pvs.name
    ps4 = ps4_strength(ann)
    if ps4:
        codes.append(_code("PS4", ps4))
        drivers["PS4"] = ps4.name

    # --- sum Tavtigian points -> band ---
    pts = 0.0
    ba1_hit = False
    for c in codes:
        p, is_ba1 = code_points(c)
        if is_ba1:
            ba1_hit = True
        pts += p
    cls = Classification.B if ba1_hit else _classify(pts)
    confidence = ("full" if spec.covered and not missing
                  else "reduced: " + (("no VCEP spec; " if not spec.covered else "")
                                      + ("missing " + ", ".join(missing) if missing else "")).strip("; "))
    return ScoredVariant(gene=gene, hgvs=hgvs, spec_id=spec.spec_id, covered=spec.covered,
                         codes=codes, points=round(pts, 2), classification=cls,
                         confidence=confidence or "full", drivers=drivers)
