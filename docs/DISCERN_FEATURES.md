# DISCERN - Full Feature Map (start to end)

**DISCERN = Diagnostic Inference from Shared-mechanism Coupling of Evidence in Rare Nosology.**

*Every feature: what it tells you, what you feed in, what you get back, how it is computed, the data behind it, and its validation status.*

**Status legend.** Validated (open data) = benchmarked against public reference data, number reported. Cohort-gated = pre-registered, needs a paired phenotype-genotype cohort to claim (Paper 2). Deferred = setup done or scoped but not run, with the reason stated. Not predicted = a measured clinical endpoint that the system never guesses. All probabilities are 0 to 1.

---

## In one paragraph (the simple version)

Inherited bleeding and platelet disorders are a set of rare diseases that look almost identical at the bedside (they present with bruising, nosebleeds, heavy periods, low platelets) but need very different treatment, and their genetic tests keep coming back "uncertain". DISCERN is a decision-support engine that takes what a clinician already has (the affected gene, a few clinical and laboratory findings, and optionally a specific variant) and returns four things, each computed by rule and never invented: a ranked differential across the look-alike diseases, a calibrated classification of the variant with the exact ACMG criteria applied, a safety alert when the likely diagnosis makes a planned treatment dangerous, and the single next test that would resolve the case. It covers 10 confusable disease clusters, 31 diseases, 29 causal genes on the clinical side, and a 170-gene bleeding/platelet panel on the variant side. Its one hard rule is no fabrication: every number is a tool-computed value with a stated scope, and when the evidence does not support a call the engine abstains and tells you which test would settle it.

## The problem it solves

- **Look-alike diseases, divergent treatment.** Type 2B von Willebrand disease and platelet-type VWD present the same way, but desmopressin helps one and harms the other. Bernard-Soulier syndrome is routinely mistaken for immune thrombocytopenia, leading to unnecessary splenectomy. Inherited thrombocytopenia with leukaemia risk (RUNX1) can be mistaken for ITP, with catastrophic consequences if a relative carrying the same variant is used as a transplant donor. Getting the differential right is treatment-changing.
- **The variant-of-uncertain-significance (VUS) problem.** Most missense variants in these genes come back as VUS. On the eRepo bleeding set, 974 of 2,239 expert-classified single-nucleotide variants (about 43 percent) are Uncertain Significance. A VUS stalls the diagnosis. DISCERN scores the variant and, where sequence evidence alone cannot cross the threshold, names the functional or phenotype evidence that would.
- **Double-counting and over-calling.** Many tools re-use the same evidence in more than one place (for example, letting the phenotype both suggest the disease and boost the variant), which manufactures false confidence. DISCERN enforces a strict partition: each piece of evidence is owned by exactly one factor.

## What you feed in (input) and what you get back (output)

**Input (all optional except the gene):**
- Gene (from the covered panel), grouped by disease cluster.
- Clinical and laboratory features, each marked present, absent (a pertinent negative), or not assessed.
- A variant: consequence plus REVEL/allele-frequency, or GRCh38 coordinates for an automatic lookup.
- A planned management action (for example, desmopressin) to check the safety interlock.
- No patient identifiers, ever (Gate G7).

**Output:**
- **Differential diagnosis:** the diseases in the relevant cluster ranked by posterior probability, with the leading call and a decided-or-abstains verdict.
- **Variant classification:** the ACMG band (Pathogenic to Benign), a calibrated probability, and the criteria applied, each shown under the single factor that owns it.
- **Safety:** any treatment-divergence hard-stop, shown first when present.
- **Recommended next test:** the one assay that best resolves the case (and, when the variant is a VUS, would supply its missing evidence).
- **Provenance:** the VCEP specifications, the likelihood ratios, and the primary references used.

---

## STAGE 1 - INTAKE (turn a case into structured evidence)

| Feature | What it tells you | How it is computed | Data / source | Output and status |
|---|---|---|---|---|
| Gene to cluster routing | which look-alike disease group the case belongs to | deterministic map from gene to one or more confusable clusters (`diseases/ontology.py`) | the curated cluster knowledge base | cluster id(s) - Validated (open data) |
| Feature intake (tri-state) | records each clinical/lab finding as present, absent, or unassessed | present is evidence for, absent is an explicit pertinent negative, unassessed is omitted (never assumed) | the per-cluster discriminating-feature list (63 features) | typed evidence set - Validated (open data) |
| Scope guard | whether the case is inside the modelled panel | genes outside the 10 clusters return a plain out-of-scope message, not a degraded guess | the covered gene list | in-scope or explicit refusal - Validated (open data) |

---

## STAGE 2 - VARIANT CLASSIFICATION (the ACMG band, calibrated, no double-counting)

| Feature | What it tells you | How it is computed | Data / source | Output and status |
|---|---|---|---|---|
| ACMG point classification | Pathogenic / Likely pathogenic / VUS / Likely benign / Benign | ClinGen-style Tavtigian point system over the applied codes (`rules/variant_scoring.py`); P at >=10 points, LP 6-9, VUS 0-5 | ACMG/AMP framework (Richards 2015); Tavtigian 2018, 2020 | class + points - Validated (open data) |
| PVS1 decision tree | the correct strength for a loss-of-function variant | Abou Tayoun 2018 decision tree (NMD, biologically-relevant transcript, critical exon) | the published PVS1 tree | PVS1 at graded strength - Validated (open data) |
| Gene-specific frequency codes | BA1 / BS1 / PM2 at the right thresholds for this gene | per-VCEP frequency thresholds; cross-checked on real gnomAD allele frequencies | ClinGen VCEP CSpecs (GT / F8 / F9 / VWF / GP1BA) + gnomAD v2.1.1/v4 | frequency code - Validated: 97.8 percent concordance on 629 curator-cited AFs |
| In-silico code (PP3/BP4) | computational support for or against pathogenicity | REVEL at ClinGen-calibrated thresholds (Pejaver 2022); AlphaMissense available | REVEL, AlphaMissense | PP3/BP4 - Validated (open data) |
| Calibrated probability | a pathogenicity probability that means what it says | Tavtigian point-to-probability posterior; and an isotonic recalibration for the reliability claim | 7,521 ClinVar-labelled variants | 0 to 1 - Validated: isotonic ECE 0.008 / Brier 0.0073 (from 0.201 / 0.060 uncalibrated) |
| No-double-counting partition | proof that no evidence is used twice | each ACMG code is assigned to exactly one factor (variant-intrinsic, functional, phenotype-coupling, segregation, phasing); a unit test enforces variant-marginal invariance | the partition rule set | per-code owning factor - Validated: 100 percent coverage on 12,240 variants / 170 genes; 33.2 percent (95 percent CI 32.4-34.1) would be over-classified without it |
| ACMG combining-rule fidelity | that the arithmetic matches expert panels | replays the experts' own applied codes and checks the resulting band | 2,653 bleeding-panel eRepo variants | Validated: 93.0 percent exact / 100 percent within-one-bin |

---

## STAGE 3 - DIFFERENTIAL DIAGNOSIS (rank the look-alike diseases)

| Feature | What it tells you | How it is computed | Data / source | Output and status |
|---|---|---|---|---|
| Cluster posterior | probability of each disease in the cluster | a factor-graph joint over disease and variant states; phenotype likelihood ratios per feature, marginalised to a per-disease posterior (`jointdx/`) | 10 clusters, 31 diseases; ~95 PMID-sourced feature likelihood ratios | ranked posteriors summing to ~1 - Validated (open data) |
| Discrimination likelihood ratios | how strongly each feature separates the diseases | `P(feature | disease)` per feature, each with sample size and a source PMID; a CI guard requires a PMID on every entry | primary literature, independently source-verified | per-feature LR - Validated: all ~95 values verified or plausible; 10 misattributed citations found and corrected |
| Curated diagnosis accuracy | how often the leading call is right | run over hand-curated, citation-only published cases (no identifiers) | 42 cases across all 10 clusters, every PMID NCBI-verified | Validated: Top-1 81 percent / Top-3 100 percent / abstention 10 percent; every non-Top-1 is a same-cluster confusable held in the Top-3 |
| Independent LOF sensitivity | that the engine catches disease alleles it never trained on | classify the null (loss-of-function) variants in the CDC F8/F9 catalogs by consequence alone | CHAMP + CHBMP, 5,437 disease alleles (2,130 null) | Validated: 91.2 percent LP/P recall by consequence alone (F8 97.7, F9 65.8, the F9 gap explained by its terminal-exon NMD rules) |

---

## STAGE 4 - DISEASE-VARIANT COUPLING (the novel core)

| Feature | What it tells you | How it is computed | Data / source | Output and status |
|---|---|---|---|---|
| PP4 coupling | lets a confirmed disease help resolve an uncertain variant | a calibrated phenotype likelihood drives the disease-to-variant link (PP4), tilting the variant posterior on-gene; each ACMG code still enters exactly one factor (the partition holds) | the joint model | reclassification signal - Cohort-gated (Paper 2), pre-registered on OSF |
| Public proof-of-concept | a first, circularity-safe test of the coupling on open data | three disjoint streams (sequence band / clinical phenotype / independent truth); matched-versus-mismatched lift; the deciding functional finding withheld from the coupling input | Phenopacket Store bleeding subset (PMID 39394689) | Reported: public working set n=2; binary lift 0.0; continuous lift +0.126 (directionally consistent). The public corpus is too thin to power the endpoint, which is what motivates the cohort study |
| Confirmatory endpoint (H6) | whether the coupling improves VUS reclassification versus intrinsic-only | pre-registered net reclassification improvement on a paired cohort, with an explicit null and the independence audit as a hard gate | BRIDGE-BPD (EGA EGAS00001001172) / an IRB cohort | Cohort-gated - Not claimed until the endpoint is reported on real paired data (Gate G13) |

---

## STAGE 5 - SAFETY INTERLOCK (treatment-divergence hard-stops)

| Feature | What it tells you | How it is computed | Data / source | Output and status |
|---|---|---|---|---|
| Hard-stop on contraindication | when a planned treatment is dangerous given the likely diagnosis | if the planned treatment is a contraindication of any non-excluded disease with non-trivial probability, a high-severity hard-stop fires (including on the leading call) | per-cluster contraindication map (DDAVP/2B, splenectomy/BSS, related-donor-transplant/RUNX1, platelet-transfusion/Quebec, recombinant-FVIII/2N) | hard-stop alert - Validated: sensitivity 100 percent (fires when it must) AND specificity 100 percent (silent when the treatment is harmless), on 5 treatment-divergent scenarios |
| Management-divergence flag | when a competitor diagnosis would change management even without a planned action | competitor-only divergence flag weighted by probability and severity | the per-cluster treatment map | competitor flag - Validated (open data) |

---

## STAGE 6 - NEXT-BEST TEST (value of information)

| Feature | What it tells you | How it is computed | Data / source | Output and status |
|---|---|---|---|---|
| Deciding-test recommendation | the single assay that best resolves the case | value-of-information over the cluster's candidate observations, each with its outcome likelihood ratios (`nextobs/`) | 30 curated next-observation assays with outcome LRs | recommended test - Validated (open data) |
| Resolves the VUS too | when the variant is uncertain, the same test that names the disease also supplies the variant's missing functional evidence | the deciding observation is linked to the variant factor it would inform | the observation-to-factor links | linked recommendation - Validated (open data) |

---

## STAGE 7 - TRUSTWORTHINESS (calibration, abstention, selective prediction)

| Feature | What it tells you | How it is computed | Data / source | Output and status |
|---|---|---|---|---|
| Variant calibration | that the probability is reliable, not just discriminative | isotonic recalibration; expected calibration error and Brier score reported | 7,521 ClinVar / eRepo variants | Validated: ECE 0.008 (ClinVar set), 0.017 (eRepo-primary), and holds on a time-split panel |
| Diagnosis calibration + confidently-wrong | that high-confidence diagnoses are not wrong | reliability of the leading-call confidence; count of confident-and-wrong calls | 42 curated cases | Validated: 0 confidently-wrong calls at the 0.8 threshold (diagnosis ECE 0.141, limited by n=42) |
| Abstention / risk-coverage | that abstaining raises accuracy on the cases kept | accuracy on the retained cases as a function of coverage | 42 curated cases | Validated: accuracy rises 81 percent at full coverage to 100 percent on the most-confident half (monotone) |
| Conformal coverage guarantee | a per-class held-out coverage guarantee | Mondrian split-conformal machinery implemented and unit-tested | needs a labelled diagnosis cohort | Implemented; empirical per-class coverage is Cohort-gated (n=42 is too small) |

---

## STAGE 8 - HEAD-TO-HEAD POSITIONING (versus current tools)

| Feature | What it tells you | How it is computed | Data / source | Output and status |
|---|---|---|---|---|
| Variant arm versus the tool set | where DISCERN stands against current ACMG tools and predictors | AUROC / calibration on the missense axis versus GeneBe, REVEL, AlphaMissense, InterVar | 342 missense P/B (ClinVar) and 425 (eRepo expert-panel) | Validated: DISCERN AUROC 0.935 (ClinVar) / 0.939 (eRepo) tracks REVEL; DISCERN is the only tool emitting a calibrated probability (ECE 0.017 vs REVEL 0.09 / AlphaMissense 0.16; GeneBe and InterVar are class-only) |
| ClinVar-circularity finding | why a ClinVar-consuming tool cannot be graded on ClinVar | GeneBe reproduces the ClinVar label (AUROC 1.0, class matches direction on all resolved calls, holds on the time-split) | the same sets | Validated finding: motivates the eRepo-primary and time-split surfaces |
| Per-code kappa versus experts | the partition made visible against an expert-panel surface | Cohen kappa of DISCERN's applied codes versus eRepo's applied codes | 2,239 eRepo variants | Validated: PVS1 kappa 0.81, PP3 kappa 0.93; PS3/PS4/PM1/PM5 applied = 0 (the evidence-stream codes DISCERN does not derive from sequence alone) |
| BIAS-2015 (fair classifier) | the strongest current fair-classifier comparator, ClinVar-blinded | full annotation pipeline stood up (Nirvana + eRepo annotated + PP5/BP6 zeroing implemented) | the eRepo set | Deferred: the tool's own preprocessing is broken in its released tags; published eRepo numbers cited instead (pathogenic sensitivity 73.99 vs InterVar 64.31) |
| Diagnosis arm versus LIRICAL/Exomiser | an external diagnosis baseline | within-cluster head-to-head attempted | public data | Deferred: unsupportable on public data (curated cases lack HPO depth; the public phenopacket corpus lacks confusable-pair density). Reported as a motivating negative; the within-cluster contest needs the cohort (Paper 2) |
| Reader study (does it help) | whether a non-specialist diagnoses better with DISCERN | a pre-registered vignette instrument (14 citation-only vignettes) plus the DISCERN-side ceiling | the vignette bank | Instrument built and tested; DISCERN-side ceiling disease 85 percent / harmful-management flagged 100 percent / deciding test 62 percent. The randomized reader arms are the user-run step |

---

## The one rule across all of it

Every number above is a tool-computed value with a stated scope. Validated axes are benchmarked on public data with the number reported (negatives included). The coupling is pre-registered and never claimed as clinically validated until its paired-data endpoint is reported (Gate G13). Anything outside scope - a patient's actual bleeding severity, an in-vivo response, a clinical outcome - is a known-unknown that the system does not predict. The product is a traceable answer, or an explicit abstention that names the test which would resolve it, never a fabricated one.
