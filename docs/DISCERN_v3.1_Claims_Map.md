# DISCERN v3.1 - Claims Map (claim -> evidence -> status)


> **Phase R amendment (pre-submission re-analysis).** Several values below are superseded. Read
> [DISCERN_PhaseR_Results_v1.md](DISCERN_PhaseR_Results_v1.md) first; it records what changed and
> why, including the results that cost the paper a claim. In short: the "only tool emitting a
> calibrated probability" claim is **retired** (put through identical folds, REVEL reaches ECE
> 0.043 and AlphaMissense 0.029 against DISCERN's 0.017, with overlapping intervals); the AUROC gap
> to REVEL is **not significant** (DeLong p=0.29); the curated diagnosis benchmark now reads Top-1
> **100%** after Phase R fixed a missing `P(G|D)` term, but must always be quoted against its
> **93%** gene-only baseline, which it does not significantly beat; the monotone risk-coverage claim
> is **withdrawn** as retired (accuracy is now 100% at full coverage, leaving no headroom); and
> LIRICAL was **run** rather than deferred.


**Purpose (Track D1):** lock every public claim to its evidence and its status, and split
the **defensible-now methods+variant paper** from the **cohort-gated coupling/clinical paper**.
Scope language is audited against Execution Plan Section A. Gate G13: no "the coupling works"
claim until paired-cohort data reports the pre-registered endpoint.

## Paper 1 - Methods + Variant engine (DEFENSIBLE NOW; needs no cohorts)

| Claim | Evidence | Status |
|---|---|---|
| The per-code partition prevents the evidence double-counting that **phenotype-aware/coupled classifiers and ClinVar-label-consuming ACMG pipelines** commit | Genome-wide: **100% partition coverage** on 12,240 variants / 170 genes; for **33.2%** (Wilson 95% CI 32.4-34.1) the non-genetic evidence already inside the VCEP bottom-line (PP4/PS3/PP1/PM3) is band-determining, so re-adding it over-classifies them; + the variant-marginal invariance unit test. **PP4 (phenotype) is the largest routed bucket (1,443 pts)** - i.e. the dominant double-count is exactly the evidence a coupling model re-adds | **DONE (open data)** |
| Gene-specific freq thresholds reproduce VCEP freq codes on real gnomAD AFs | gnomAD per-variant cross-check: **97.8%** concordance on 629 variants (curator-cited gnomAD AFs) | **DONE** |
| ACMG combining-rule fidelity | **93.0% exact / 100% within-one-bin** on 2,653 bleeding-panel eRepo variants (arithmetic given experts' codes) | **DONE** |
| Intrinsic-only band concords with the ClinVar community label | **62.4% exact / 92.8% within-one-bin** on 10,780 matched; intrinsic-only is a *designed lower bound* (omits routed PP4/PS3/PP1/PM3) - disagreements characterised, not errors | **DONE (characterised)** |
| Variant probabilities are **calibratable** (the deliverable is calibration, NOT discrimination) | Isotonic **ECE 0.008 / Brier 0.0073** (vs uncalibrated 0.201/0.060) on 7,521 ClinVar-labelled variants. *AUC 0.999 is on P/B extremes with VUS dropped - discrimination is trivially easy there; it is reported only to confirm the isotonic map preserved ranking (monotonic transform), not as a discrimination claim.* | **DONE (H5 variant half)** |
| Gene-specific CSpec criteria extracted | BA1/BS1/PM2 + PM2_Supporting for GT/F8/F9/VWF/GP1BA (CSpec GN071/GN079/GN081 + eRepo) | **DONE** |
| Novel-variant scoring vs InterVar (H4, **full-DB, 2026-06-17**) | Literal InterVar re-run with its **complete default database set** (1000g2015aug, esp6500, avsnp147, dbnsfp42a, dbscsnv11, clinvar_20210501, rmsk, gnomAD-exome + full intervardb/OMIM - nothing dropped) on the same 1,015 ClinVar P/B variants. **Missense subset (n=364, the meaningful VUS axis): DISCERN-full AUC 0.944 DECISIVELY beats literal InterVar 0.811 (+0.133) and matches REVEL 0.942** (which DISCERN ingests as PP3). **All P/B variants (dominated by easy null/splice - the "trivially easy extremes"): DISCERN 0.882 (with frequency) / 0.912 (without) ~ InterVar 0.887 - comparable.** | **DONE (open data)** - `eval/intervar_full_eval.py` |
| ~~Overall H4 AUC edge over InterVar~~ (**CORRECTED**) | The earlier reduced-DB result ("DISCERN 0.912 > InterVar 0.874, +0.038 overall") **overstated the overall edge**: it used a handicapped InterVar (dropped 1000g/avsnp/mim2gene) AND scored DISCERN without frequency. The fair full-DB run shows the overall AUCs converge; the real, large DISCERN advantage is **on missense**, consistent with the "all-variant P/B AUC is trivially easy" caveat. | **CORRECTED - supersedes the reduced-DB claim** |
| Independent (non-ClinGen) sensitivity on LOF disease alleles | **CHAMP/CHBMP (CDC F8/F9 catalogs, 5,437 disease alleles): 91.2% LP/P recall on the 2,130 null variants by CONSEQUENCE ALONE** (no gnomAD/predictor); F8 97.7%, F9 65.8% (F9 gap fully explained by its large terminal exon / NMD-escape rules). | **DONE (open data)** - `eval/champ_chbmp_benchmark.py` |
| Novel-missense recovery is bounded by routed evidence (motivates coupling) | **Missense arm (REVEL-finalized):** on novel F8/F9 missense (cDNA not-in ClinVar), variant+ClinVar engine recovers to LP only the PS1andPM2andPP3 set - **F8 16/1307=1.2%, F9 11/511=2.2%, combined 27/1818=1.5%** (REVEL via ANNOVAR hg19_revel on VM; 28/28 gnomAD-absent confirmed). PM5 (26-47%) and no-hit (50-72%) stay VUS by point arithmetic. Intrinsic to ACMG (not DISCERN-specific); ~98.5% of novel missense need PS3/PP1/**PP4 coupling**. | **DONE (open data)** - `eval/champ_chbmp_missense_arm.py` |
| Full confusable-cluster coverage C1-C10 | 10 clusters; every discrimination LR PMID-sourced (CI-guarded) **and source-verified 2026-06-17** (all ~95 values verified or plausible; 10 misattributed citations corrected); per-cluster safety/contraindication map | **DONE (B)** - LR *calibration* (vs a labeled cohort) is Paper 2 |
| Management-aware safety incl. leading-call hard-stop | interlock fix + per-cluster contraindications (DDAVP/2B, splenectomy/inherited-IT, rFXIII-A2/F13B, platelets/Quebec, HSCT) + regression tests | **DONE (B)** |
| Selective/conformal prediction (machinery) | Mondrian split-conformal implemented (`jointdx/conformal.py`) + mechanics test (synthetic-sampling guarantee) | **IMPLEMENTED (Paper 1)** - the *empirical coverage guarantee* (held-out coverage = nominal 1-alpha) is **NOT claimed in Paper 1**; scoped to **Paper 2** (needs a labeled diagnosis cohort; the curated n=42 is too small for per-class calibration) |
| Variant arm vs the current ACMG/predictor tool set (Track 1) | On the missense P/B axis (n=342, GeneBe-annotated, ClinVar-blinded): DISCERN AUROC 0.935 tracks REVEL 0.948 (its PP3 input) and clears InterVar 0.811 / AlphaMissense 0.921. **RETIRED by Phase R R2** (was: "the only tool emitting a calibrated probability"). Under an identical fold and isotonic protocol the comparators calibrate too - REVEL ECE 0.043, AlphaMissense 0.029, DISCERN 0.017, intervals overlapping. GeneBe and InterVar remain class-only and uncalibratable. The differentiator is the calibrated, auditable **classification system** (band + criterion trail + partition + abstention + safety hard-stop), not calibration alone. GeneBe's apparent AUROC 1.0 is ClinVar reproduction (class matches ClinVar direction 337/337), demonstrating empirically why ClinVar is not a fair truth surface for ClinVar-consuming tools. | **DONE (open data)** - `bench/track1_*`, `DISCERN_Benchmark_Results_v1.md` (Gate G-T1) |
| Trustworthiness layer as named contributions (Track 3) | Calibration (variant ECE 0.017; diagnosis ECE 0.141 on n=42 with **0 confidently-wrong**); safety hard-stop **sensitivity 100% AND specificity 100%** on 5 treatment-divergent scenarios; risk-coverage **withdrawn by Phase R R4**: with the diagnosis arm now at 100% at full coverage there is no headroom for abstention to demonstrate, so the curve is uninformative on n=42. | **DONE (open data)** - `bench/track3_*`, `DISCERN_Benchmark_Results_v1.md` (Gate G-T3) |

## Paper 2 - Coupling + Clinical (COHORT-GATED; pre-registered)

| Claim | Evidence path | Status |
|---|---|---|
| The coupling improves VUS reclassification over intrinsic-only | Pre-registered endpoint H6 (+ explicit falsification), on a controlled-access paired phenotype-genotype cohort | **NOT CLAIMED** until paired data (G12 reg, G13 gate) |
| Public coupling proof-of-concept (precursor; de-risks, not confirmatory) | `eval/coupling_poc.py` on the Phenopacket Store bleeding subset: circularity-safe matched-vs-mismatched lift. Public working set n=2; binary lift 0.0; continuous lift +0.126 (matched 0.208 vs mismatched 0.082, directionally consistent). Endpoint pre-specified in the OSF protocol before running. | **DONE (open data); reported** - public corpus too thin to power the endpoint, motivating the cohort study; does not claim the coupling works (G13 unchanged) |
| Misdiagnosis rescue (label-hidden case-control) | `eval/misdx_rescue.py` on the ITP cohort / BRIDGE-BPD | gated (DAC) |
| Diagnosis accuracy at scale | curated cases **Top-1 100% / Top-3 100% / abstention 7% (n=42)** after the Phase R `P(G|D)` fix (81% before it). **Not quotable alone:** a phenotype-blind gene lookup scores **93%** on the same cases, delta +7% (95% CI 0 to +17), McNemar p=0.25, and the subset where the gene does not settle the answer is **n=3**. LIRICAL restricted to the same cluster reaches Recall@1 57% on the 23 HPO-runnable cases. Cohorts carry the headline. | **partial (A); gated (B/C)** - the arm rests on abstention, zero confident errors and the safety interlock, not on Top-1 |
| Per-patient VUS-reclassification rate vs 3-star truth | needs paired phenotype+variant | gated (B/C) |
| Reader study (usefulness) | pre-registered vignette protocol + **built/tested instrument** (`eval/reader_study.py`, 14-vignette citation-only bank, `tests/test_reader_study.py`); DISCERN-side ceiling disease 85% / harmful-management flagged 100% / deciding test 62% | **own (pre-reg)** - instrument DONE; randomized reader arms are user-run after OSF lock |

## Scope-language rules (audited)
"Partition / arithmetic fidelity / calibration / cluster coverage" = **proven on open data**.
"Coupling" = **pending or negative**, never "validated", until Paper 2's paired-data result.
"Coverage" (cluster breadth) is never conflated with "accuracy". Every synthetic result is
labelled sanity, never a headline. No patient data appears in any public artifact (G7).
