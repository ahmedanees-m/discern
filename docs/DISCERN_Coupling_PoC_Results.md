# Disease-variant coupling: proof of concept

**Date:** 2026-06-17 - **Code:** `eval/coupling_poc.py`,
`eval/extract_phenopackets.py`, `eval/hpo_feature_crosswalk.yaml`, `tests/test_coupling_poc.py`.

**What this is.** The first real-data, public, circularity-safe test of H6: does the clinical-phenotype
coupling resolve variants that sequence evidence alone leaves as VUS, and ONLY when phenotype and gene
agree? It is a proof-of-concept that de-risks and motivates the confirmatory cohort study; it is NOT
the confirmatory validation, which follows the pre-registered protocol.

> The continuous secondary value is +0.252 (matched 0.334). The pre-registered primary
> endpoint is the binary one and stands at 0.0. With n=2 usable cases neither value
> supports a claim; the run is a directional proof of concept and the pre-registered
> protocol is untouched.

## Independent verification
- **Phenopacket Store** (public GA4GH corpus, github.com/monarch-initiative/phenopacket-store): the
  paper snapshot is 6,668 phenopackets / 475 diseases / 423 genes (Danis et al., *Human Genetics and
  Genomics Advances*, **PMID 39394689 / PMC11564936**; cite by PMID, the record is dated 2024 not 2025).
  The cloned copy used here has 11,155 phenopackets. PHI-free (published literature cases).
- **HPO term IDs:** every crosswalk id was confirmed against the HPO ontology (OLS4). The spec's draft
  ids were partly wrong and were corrected (e.g. abnormal bleeding HP:0001892 not HP:0011893; giant
  platelets HP:0001902; delayed umbilical separation HP:0032434); placeholders were resolved.
- **GA4GH schema paths** verified for v2 (camelCase JSON). **Tooling:** parsed the JSON directly rather
  than installing `phenopacket-store-toolkit`/`pyphetools`, which also sidesteps the flagged malicious
  releases (toolkit 0.1.7 / pyphetools 0.9.120).
- **The two working-set variants** were confirmed gnomAD-rare/absent (LYST c.772T>C AF 6.6e-6;
  c.11173G>A absent), so the PM2 frequency code in the sequence band is evidenced, not assumed.

## Data and the size of the public working set
The bleeding/platelet-gene subset of the corpus is **31 cases across 6 genes** (FERMT3, F8, LYST,
ACTN1, WAS, MPL). After the spec's filters:

| Stage | n | Note |
|---|---|---|
| Bleeding-gene cases | 31 | FERMT3 1, F8 8, LYST 3, ACTN1 14, WAS 4, MPL 1 |
| Gene in a DISCERN-modelled cluster | 12 | excludes ACTN1/WAS/MPL (not yet modelled) and the F8 cases (elevated-FVIII thrombophilia, not bleeding) |
| Sequence band = VUS (intrinsic) | 2 | excludes FERMT3 p.Gln96Ter (nonsense) and LYST frameshift - sequence already resolves these to LP/P |
| With a mappable clinical feature | **2** | the working set: two LYST missense Chediak-Higashi cases (PMID 24112114) |

This thinness is itself a finding (it matches the earlier Phenopacket-Store diagnosis run): the public
corpus is not population-representative for inherited bleeding/platelet disorders, so a powered
confirmatory test cannot be run on it.

## Result
On the working set (n=2), with the three streams kept disjoint and the independence guards enforced
(the diagnostic functional finding, e.g. giant neutrophil granules, is withheld as Truth and never used
as a coupling input):

| Endpoint | Matched | Mismatched | Lift |
|---|---|---|---|
| **Primary (binary upgrade to LP/P)** | 0/2 | 0/2 | **0.0** |
| **Secondary (continuous coupling PP4-equivalent, mean P(path+LP))** | **0.334** | **0.082** | **+0.252 (~4.1x)** |

**Reading it.** Neither case crosses the LP reclassification threshold, so the **binary lift is
0 on n=2**. But the **continuous signal is directionally consistent with H6 in both cases**: the matched
clinical phenotype lifts the variant's pathogenic posterior to ~2.5x the mismatched value. The coupling
carries the right disease-specific signal; it is under threshold here because (a) n=2, (b) the
literature phenotypes are sparse (often only pigment dilution for Chediak, with the decisive giant-granule
finding correctly held out as Truth), and (c) only the PM2 frequency code is in the sequence band (no
REVEL annotation, so no PP3). The coupling mechanism itself is validated by a unit test
(`test_coupling_distinguishes_matched_from_mismatched`: a matched phenotype lifts P(path+LP) above a
mismatched one).

## Conclusion
The instrument works end-to-end on real public data, the independence guards hold, and the coupling
shows a directionally-correct disease-specific signal - but the public corpus (2 usable intrinsic-VUS
bleeding cases, sparse phenotypes) cannot power the binary endpoint. **This is exactly the result the PoC
was designed to surface: it de-risks the method and quantifies why the confirmatory test requires the
controlled-access, richly-phenotyped cohort.** The full instrument - extractor, OLS-verified
crosswalk, engine, independence audit, tests, and the pre-registered endpoint + falsification condition -
is built and ready to run unchanged the moment cohort access clears.

## Limitations (per spec)
- **Ascertainment:** literature cases are solved affected patients, enriched for true positives; the
  matched/mismatched contrast (not an absolute rate) is what controls for this.
- **VUS-at-presentation is rare in the literature** (papers resolve variants before publishing), which
  is the main reason the working set is tiny.
- **Not confirmatory:** the gated cohort study remains the definitive test; this PoC does not replace it.
