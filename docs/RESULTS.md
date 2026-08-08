# Results

Every number DISCERN reports, with the dataset it was measured on and the harness that produces
it. Each harness writes a JSON file; the test suite reads those files and fails if a value quoted
here or in the README drifts from its source.

All evaluation uses public data. No patient-level data appears in this repository or in any
artifact derived from it.

## Variant classification

Surface: variants in the ClinGen Evidence Repository carrying an expert-panel classification,
restricted to the bleeding and platelet gene set. The primary panel is the full set; the
time-split panel holds out records approved after the split date, on which no tool could have
observed the expert classification.

| Measure | Surface | n | Result |
|---|---|---|---|
| AUROC, DISCERN | eRepo-primary, missense | 425 | 0.939 (95% CI 0.906 to 0.967) |
| AUROC, REVEL | eRepo-primary, missense | 425 | 0.954 (95% CI 0.931 to 0.972) |
| AUROC, AlphaMissense | eRepo-primary, missense | 424 | 0.932 (95% CI 0.906 to 0.953) |
| AUROC, DISCERN | time-split, missense | 383 | 0.927 (95% CI 0.882 to 0.961) |
| ECE, DISCERN | eRepo-primary | 425 | 0.017 (95% CI 0.011 to 0.047) |
| ECE, REVEL calibrated | eRepo-primary | 425 | 0.043 (95% CI 0.030 to 0.073) |
| ECE, AlphaMissense calibrated | eRepo-primary | 424 | 0.029 (95% CI 0.024 to 0.065) |
| ACMG combining-rule fidelity | eRepo bleeding panel | 2653 | 93.0% exact, 100% within one bin |
| Frequency threshold agreement | curator-cited gnomAD frequencies | 629 | 97.8% |

Comparators are calibrated under the identical out-of-fold isotonic protocol on shared folds, so
the calibration comparison is like for like. The paired AUROC difference between DISCERN and REVEL
is not significant (DeLong p = 0.29). DISCERN emits a calibrated probability as delivered, whereas
the comparators require a post-hoc fit against labels; that is a difference in what the tool
provides, not in discrimination.

Harnesses: `bench/phase_r_variant.py`, `bench/track1b_erepo_headtohead.py`,
`eval/erepo_reconstruction.py`.

## Evidence partition

Every ACMG/AMP criterion is assigned to exactly one factor of the joint model. Coverage is
measured genome-wide rather than on the bleeding panel alone.

| Measure | Dataset | Result |
|---|---|---|
| Partition coverage | eRepo, 12240 variants across 170 genes | 38802 applied criteria, all routed, none unrecognised |
| Bottom-line reuse | same surface | 33.2% of variants band-determined by non-genetic evidence (95% CI 32.4 to 34.1) |

Exclusivity is enforced by construction and asserted by a test that supplies the joint model a
variant carrying PP4, PS3, PP1 and PM3 and requires the variant marginal to be unchanged to
within 1e-12.

Harness: `eval/erepo_genomewide.py`. Partition: `rules/vcep/partition.py`.

## The ceiling on intrinsic evidence

Restricted to evidence a sequence-only pipeline can supply, no expert-classified missense variant
in the panel reaches a pathogenic band.

| Measure | n | Result |
|---|---|---|
| Missense variants reaching Likely Pathogenic on intrinsic evidence | 425 | 0 |
| Highest total reached | 425 | 3 points, against a 6-point threshold |
| Recovered by restoring variant-intrinsic criteria with no input here | 316 pathogenic | 27 |
| Recovered by restoring criteria the partition routes to other factors | 316 pathogenic | 52 |
| Recovered by restoring both | 316 pathogenic | 89 |

Neither cause accounts for the ceiling on its own. PM1, PM5, PS1 and PS4 remain variant-intrinsic
under the partition and have no input in this pipeline, which distinguishes an attribution
to the partition from an attribution to missing input.

Harness: `bench/phase_r_variant.py`.

## Agreement with expert criterion applications

Cohen's kappa is computed on the 1265 variants carrying a pathogenic or benign expert
classification, which is the subset on which both sides have a comparable decision. An unmodified
code is resolved to its criterion's default strength before comparison, so PP3 and PP3_Supporting
are treated as the same assertion.

Agreement is high where the input is available and the rule is mechanical (PVS1 0.98, PP3 0.95)
and low where the panels apply a criterion this pipeline has no input for. Coverage and the
direction of each disagreement are reported per criterion rather than as a single summary.

Harness: `bench/track1b_erepo_headtohead.py`.

## Independent sensitivity

The CDC hemophilia mutation projects supply an F8 and F9 truth set independent of the expert-panel
surface. Recall is scored on consequence alone, with no frequency or predictor input.

| Gene | Result |
|---|---|
| F8 | 97.7% |
| F9 | 65.8% |
| Overall | 91.2% on 2130 null variants of 5437 disease alleles |

The F9 figure reflects that gene's large terminal exon, where premature termination codons escape
nonsense-mediated decay and are held at uncertain significance rather than called pathogenic.

Harness: `eval/champ_chbmp_benchmark.py`. Detail: `DISCERN_CHAMP_CHBMP_Benchmark_v1.md`.

## Differential diagnosis

Benchmark: 42 published cases spanning all ten discrimination clusters, each represented as a
PubMed identifier, the causal gene and the extracted clinical and laboratory findings. No article
text is reproduced.

| Arm | n | Result |
|---|---|---|
| DISCERN, Top-1 | 42 | 100% |
| Phenotype-blind gene lookup, Top-1 | 42 | 93% |
| Prior only, Top-1 | 42 | 21% |
| Uniform random floor, Top-1 | 42 | 32% |

The causal gene is supplied as input, so the phenotype-blind gene lookup is the baseline that
matters: the difference against it is not significant (exact McNemar p = 0.25). The result is
in-sample, because these cases were used while developing the joint model; the same benchmark
scores 81% under the configuration that preceded the disease-posterior correction. It is reported
as a characterization of the implementation, not as a held-out comparison.

Against LIRICAL 2.4.1 on the 23 cases carrying any phenotype term, two comparisons support
different claims and are not combined:

| Comparison | Evidence | DISCERN | LIRICAL | Paired test |
|---|---|---|---|---|
| Reasoning on identical evidence | the 13 HPO-expressible findings, gene withheld from both | 74% | 57% | McNemar p = 0.29 |
| Encoding the evidence that decides these cases | all 48 findings | 91% | 57% | McNemar p = 0.02 |

The first does not establish an advantage in inference at this sample size. The second is
significant but is a claim about what the model can represent: 35 of the 48 discriminating
features are laboratory assay results that the Human Phenotype Ontology does not encode, so a
phenotype-driven ranker is denied the discriminating channel rather than handicapped equally.

Exomiser was not run. It requires a VCF, and seeding one per case with zygosity matched to each
report's stated inheritance introduces a correctness risk larger than the comparison's marginal
value, given that LIRICAL had already established the phenotype-channel result on the same inputs.

Harnesses: `eval/curated_case_benchmark.py`, `eval/gene_only_baseline.py`, `eval/lirical_arm.py`,
`eval/phenotype_tool_comparison.py`.

## Treatment safety

The interlock fires when a planned management is contraindicated by any disease the evidence has
not excluded, including the leading call. Adjudication uses a gene-blind posterior: sequencing one
gene does not exclude a disease of another gene, which may never have been sequenced.

| Measure | n | Result |
|---|---|---|
| Sensitivity, fires on contraindicated management | 5 scenarios | 100% |
| Specificity, silent on harmless management | 5 scenarios | 100% |

Behaviour is invariant across the full range of the gene-term likelihood ratio.

Harness: `bench/track3_trustworthiness.py`. Implementation: `safety/interlock.py`.

## Abstention and calibration of the diagnosis posterior

On the curated benchmark the diagnosis arm makes no confidently wrong call above a 0.8 confidence
threshold (0 of 42), with mean predicted confidence 0.882 against observed accuracy 1.000 and an
expected calibration error of 0.118. The arm sits at ceiling at full coverage, so this
characterizes calibration rather than demonstrating a coverage-accuracy trade-off; a monotone
risk-coverage result is not claimed and would require a benchmark that is not at ceiling.

Harness: `bench/track3_trustworthiness.py`.

## Disease-variant coupling

The coupling that motivates the architecture is implemented and unit-tested. Its evaluation on
real paired phenotype-genotype data is pre-registered with an explicit falsification condition and
is not claimed here. Public corpora do not contain enough cases carrying both a classified variant
and the laboratory findings that separate the confusable diseases; a proof-of-concept over the
GA4GH Phenopacket Store returns a binary reclassification lift of 0.0, which is the expected
result when the discriminating channel is absent from the corpus rather than evidence against the
mechanism.

Harness: `eval/coupling_poc.py`. Detail: `DISCERN_Coupling_PoC_Results.md`. Protocol:
`DISCERN_OSF_PreRegistration_v1.md`.

## Scope of claims

The following are measured and reported above: discrimination and calibration on the expert-panel
surface, partition coverage and exclusivity, the intrinsic-evidence ceiling and its attribution,
independent sensitivity on the CDC catalogs, differential diagnosis against phenotype-blind and
external baselines, and the treatment-safety interlock.

The following are not claimed. DISCERN is not the only tool that can emit a calibrated
probability: under identical folds and protocol, REVEL and AlphaMissense calibrate to overlapping
intervals. The AUROC advantage over REVEL is not established. A monotone risk-coverage
relationship is not claimed. An advantage in reasoning over a phenotype-driven ranker on identical
evidence is not established at this sample size. The disease-variant coupling endpoint is not
claimed until its pre-registered evaluation is reported.
