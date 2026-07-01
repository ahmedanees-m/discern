# DISCERN - Pre-Registration Protocol v1 (DRAFT for OSF)

**Status:** DRAFT. To be time-stamped on the Open Science Framework (osf.io) **before any
cohort data is analysed** (Gate G12). No paired-cohort data has been examined at the time of
registration. **Author:** Anees Ahmed Mahaboob Ali (`ahmedanees-m`). **Date drafted:** 2026-06-13.

## 1. Background & the claim under test
DISCERN computes one coupled posterior `P(D,V|E)` over disease (D) and variant (V). The novel,
**unproven** element is the COUPLING: a calibrated phenotype likelihood drives PP4 (the
disease->variant link) and thereby variant reclassification. The per-code partition (no
double-counting) and the variant-intrinsic layer are already validated on open data and are NOT
the subject of this registration. This protocol fixes, in advance, how the coupling will be
tested so the result is credible whether positive or negative.

## 2. Primary hypothesis and the explicit falsification condition
- **H6 (primary).** On paired phenotype+variant data, the coupled model reclassifies input VUS
  toward the ClinVar/VCEP 2-star+ truth label **more accurately than the variant-intrinsic-only
  model**, with the circularity guard active (each ACMG code enters exactly one factor).
- **Primary endpoint.** Net reclassification improvement of coupled vs intrinsic-only:
  (correctly reclassified VUS) - (incorrectly reclassified VUS), per patient, truth = ClinVar
  2-star+ / VCEP. Pre-specified superiority margin: coupled - intrinsic-only > 0 with a 95% CI
  excluding 0 (paired bootstrap, 10,000 resamples).
- **FALSIFICATION (reported as the headline negative result).** If the 95% CI includes or lies
  below 0, the coupling does NOT improve reclassification on real data; this is reported as the
  primary finding. No post-hoc redefinition of the endpoint is permitted.

## 3. Secondary hypotheses
- **H7a (misdiagnosis rescue).** On documented mislabelled cases (e.g. inherited thrombocytopenia
  labelled ITP; 2B labelled PT-VWD), DISCERN flags the treatment-changing competitor from
  PRE-CORRECTION evidence (corrected label hidden; Gate G6). Endpoint: flag sensitivity/specificity.
- **H7b (safety).** Confident-and-wrong rate on the confusable pairs is the headline safety
  metric (must be reported even if other metrics are favourable; Gate G4).
- **H5 (coverage, cohort half).** Mondrian split-conformal empirical coverage on real labels is
  within +/-tolerance of (1-alpha) per class.

## 4. Datasets (verified; access gated)
- **BRIDGE-BPD** EGA `EGAS00001001172` ("Sequencing of heritable Bleeding and Platelet
  Disorders"; Downes 2019, PMID 31064749) - NIHR BioResource DAC.
- **80-patient ITP cohort** (Joshi/Cooper, Blood Adv 2025;9(7):1497, PMID 39808791) - DAC.
- **a controlled-access paired cohort** - controlled access, aggregate-only, no patient data in public artifacts (Gate G7).
- Truth labels: ClinVar `variant_summary` review-status 2-star+; VCEP classifications.

## 5. Pre-specified analysis
- Circularity guard ON: PP4/PS3/PM3/PP1 owned by their single factor; verified by the
  invariance unit test + the genome-wide partition (100% coverage) before any run.
- Abstention computed BEFORE accuracy is claimed (G4); abstention threshold frozen from a
  held-out calibration split, not tuned on the test set.
- Primary: paired bootstrap NRI (coupled vs intrinsic-only). Secondary: per-class conformal
  coverage; flag sensitivity/specificity with Wilson CIs; reliability diagram + ECE + Brier.
- No subgroup or endpoint added after data inspection; any exploratory analysis labelled as such.

## 6. Pre-specified figures & tables
T1 cohort description; T2 primary NRI (coupled vs intrinsic-only, 95% CI); T3 misdiagnosis-rescue
sensitivity/specificity; F1 reliability diagram; F2 conformal coverage per class; F3 confident-
and-wrong rate. The negative-result version of every figure/table is pre-drawn.

## 7. Reader study (companion)
Pre-registered vignette bank (`eval/cases/reader_vignettes.yaml`, 14 citation-only vignettes: 6
treatment-divergent safety pairs, 7 deciding-test, 1 acquired-mimic distractor); non-specialists
read with vs without the DISCERN card. Co-primary endpoints = diagnostic accuracy and
time-to-decision; secondary = harmful-management rate (headline safety, Gate G4), correct deciding
test, confidence calibration; analysed only once after lock. Harness + tests built
(`eval/reader_study.py`, `tests/test_reader_study.py`); the DISCERN-side ceiling on the bank is
disease 85% / harmful-management flagged 100% / deciding test 62%.

## 8. Negative-result & stopping commitments
The result is published regardless of direction. No optional stopping; the analysis runs once on
the locked dataset after registration. Gate G13: no artifact will claim "the coupling works /
improves reclassification" until this protocol's primary endpoint is reported on real paired data.

*Dry run:* `eval/synthetic_coupling_harness.py` executes this pipeline on SIMULATED data to
confirm the analysis code and the circularity guard before real data (synthetic = sanity, not evidence).

## 9. Public-data Coupling Proof-of-Concept (pre-specified precursor; `eval/coupling_poc.py`)
A circularity-safe public-data precursor on the Phenopacket Store bleeding subset, pre-specified here
before the run so its result counts as a pre-registered analysis, not a fished one. It de-risks but
does NOT replace the confirmatory cohort study (Gate G13 unchanged).
- **Three disjoint streams:** Sequence (variant-intrinsic band) / Clinical-phenotype (coupling input)
  / Truth (independent functional finding or published causal status). A measurement is assigned to
  exactly one stream; the diagnostic functional finding is withheld from the coupling input.
- **Working set:** cases whose Sequence band is VUS (the only cases the coupling can act on).
- **Primary endpoint:** lift = upgrade_rate(matched) - upgrade_rate(mismatched), where matched uses the
  case's real clinical phenotype and mismatched uses a sibling disease's profile (a gene the variant is
  not in). **Secondary:** the continuous coupling PP4-equivalent P(path+LP), matched vs mismatched.
- **Falsification (stated before running):** lift ~ 0 (matched and mismatched upgrade at the same rate)
  means the coupling carries no disease-specific signal and H6 is not supported on this set. Reported
  regardless of direction, including the working-set n and all exclusions.
- **Result (2026-06-17, reported per this commitment):** public working set n=2; binary lift 0.0;
  secondary continuous lift +0.126 (matched 0.208 vs mismatched 0.082, directionally consistent). The
  public corpus is too thin to power the binary endpoint - which is the finding that motivates the
  cohort study. Full writeup: `DISCERN_Coupling_PoC_Results.md`.
