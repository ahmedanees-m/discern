# Supplemental Methods

Detail supporting the Material and methods section of the main text. Nothing here is required to
follow the argument; everything here is required to reproduce it. Every procedure below is
implemented in the archived code, and the file that implements it is named at the end of each
subsection.

---

## Variant Curation Expert Panel specifications as implemented

Gene-specific rule sets are encoded as versioned YAML, one file per specification, loaded at run
time rather than compiled into the scoring logic. Four elements are captured per specification: the
frequency thresholds that determine BA1, BS1 and PM2; any modification to a criterion's default
strength; the computational thresholds that determine PP3 and BP4; and the identifier of the
specification itself, so a classification can be traced to the rule set that produced it.

**Frequency criteria.** Thresholds were taken from the ClinGen Criteria Specification Registry
where a specification exists, and cross-checked against the allele frequencies curators cite in the
corresponding Evidence Repository records. This cross-check corrected several values during
development and is reported in the main text as 97.8% concordance on 629 variants for which a
curator-cited frequency was available.

**Strength modification.** Where a panel modifies a criterion's default strength, the modification
is encoded explicitly. The most consequential is `PM2` at Supporting rather than the ACMG default of
Moderate, which several panels specify and which materially changes the point total for rare
missense variants, PM2 being the most frequently applied criterion in this domain.

**Computational thresholds.** Each specification carries its own REVEL and SpliceAI cut-offs for
PP3 and BP4 where the panel defines them. Where a panel does not, the ClinGen-calibrated bands are
used as a fallback. The RUNX1 specification uses the Myeloid Malignancy panel's version 2 values
(PP3 at REVEL >= 0.88, BP4 at REVEL <= 0.50) rather than generic defaults.

**Residual simplifications, stated.** Two elements of the published specifications are not fully
encoded. PVS1 is implemented as the Abou Tayoun decision tree rather than as each panel's
gene-specific elaboration of it, and PS4 is implemented as a generic odds-ratio and proband-count
tree rather than each panel's semi-quantitative counting rules. Both are documented simplifications
rather than omissions; neither affects the results reported here, since PS4 is never applied for
want of case-control input and PVS1 does not apply to the missense surface that carries the
analysis.

*Implemented in `rules/vcep/loader.py` and `rules/vcep/specs/*.yaml`.*

## The evidence partition

Each ACMG/AMP criterion is mapped to exactly one owning factor by a static table. Criterion strings
are normalised before lookup, so a strength-modified form such as `PM2_Supporting` resolves to the
same owner as `PM2`. The mapping is total over the criterion vocabulary encountered in the Evidence
Repository: the genome-wide analysis reports 100.0000% coverage over 38,802 applied criteria with
no unrecognised code.

Exclusivity is enforced by construction rather than by convention, and verified by a test that
supplies the joint model a variant with the criteria PP4, PS3, PP1 and PM3 added explicitly and
asserts the variant marginal is unchanged to within 1e-12. Because those criteria are owned by
other factors, re-adding them to the genetic stream must have no effect; the test fails if the
partition is ever weakened.

*Implemented in `rules/vcep/partition.py`; verified by `tests/test_discern_joint.py`.*

## The gene term and its two bounds

The joint factorises as

    P(D, V | E)  proportional to  P(E_pheno | D) . P(E_geno | V) . P(E_func | D, V)
                                  . P(V | G, D) . P(G | D) . P(D)

`P(V | G, D)` normalises over the five variant states, so it cancels when the joint is marginalised
to a disease posterior. Without the separate `P(G | D)` term the gene therefore carries no
information about disease identity, which is the defect described in the main text. `P(G | D)` is
implemented as an on-gene/off-gene likelihood pair, giving a likelihood ratio of 4.

Two bounds constrain that value, and both are architectural rather than fitted:

1. **It must not exceed a cluster's deciding observation.** The sharpest discriminator in the
   knowledge base is the ristocetin-induced platelet aggregation mixing study, at a likelihood ratio
   of approximately 11.5. If the gene term exceeded that, no laboratory result could overturn the
   gene, and the value-of-information layer could never predict a change of leading diagnosis. A
   ratio of 4 is the largest tested value at which the mixing study still overturns the gene.
2. **It must not veto a treatment-divergent competitor.** Finding a variant in one gene does not
   exclude a disease of another gene, which may never have been sequenced or whose variant may be
   benign. The safety interlock is therefore evaluated on a gene-blind posterior, obtained by
   omitting `P(G | D)`, while the leading diagnosis reported to the user remains the gene-aware one.

A sweep across the full range of the parameter (Figure S1) shows hard-stop sensitivity and
specificity invariant at 100% throughout, the curated diagnosis result flat above the committed
value, and the deciding-assay bound failing at ratios of 9 and above.

*Implemented in `jointdx/factorgraph.py`; swept by `bench/phase_r_gene_term_sensitivity.py`.*

## Calibration and its out-of-sample estimation

Variant scores are calibrated by isotonic regression. Estimation is out-of-fold throughout: a
5-fold stratified split with shuffling and a fixed seed is constructed once per surface; the
isotonic map is fit on the four training folds and applied to the held-out fold; expected
calibration error and Brier score are accumulated on held-out predictions only. No variant is ever
scored by a calibrator that has seen its label.

Three properties are asserted in code at run time rather than checked afterwards: that no fold
trains and evaluates on the same index, that no index appears in two evaluation folds, and that
every index is held out exactly once. Any violation raises rather than returning a number.

The fold assignment itself is published (`benchmarks/variant_arm/calibration_folds.json` in the data
record), listing for every variant which fold held it out, together with the full training and
evaluation index sets. A reader can verify the three properties directly rather than taking the
protocol on trust.

Expected calibration error uses ten equal-width bins over the unit interval, weighting each bin's
absolute gap between mean predicted probability and observed frequency by that bin's occupancy. The
final bin includes probability 1.0.

**Comparator calibration.** REVEL and AlphaMissense were put through the identical fold assignment
and the identical isotonic procedure, and are reported both raw and calibrated. This is what makes
the calibration comparison like-for-like, and it is what caused the calibration-uniqueness claim of
an earlier draft to be retired.

*Implemented in `bench/phase_r_variant.py`; the calibration primitive is `core/stats.ece`.*

## Statistical procedures

**Confidence intervals on proportions.** Wilson score intervals. Preferred to the normal
approximation because the reported rates sit far from one half at large sample size, where Wald
intervals are unreliable.

**Confidence intervals on metrics.** Percentile bootstrap, 1,000 resamples, fixed seed. Resamples
in which a class is absent are discarded rather than scored. Intervals are the 2.5th and 97.5th
percentiles of the resampling distribution.

**Comparing two areas under the curve.** DeLong's test in the fast form of Sun and Xu, computed on
the variants both tools score, with ties handled by mid-ranks. Reported alongside a paired bootstrap
on the difference, which resamples variants and rescores both tools on the same resample, so the
interval reflects the pairing.

**Comparing two classifiers on the same cases.** Exact McNemar, using the binomial test on the
discordant pairs. Chosen over the chi-squared approximation because the discordant counts here are
small.

**Rank metrics.** Recall at k and mean reciprocal rank from one-based hit positions, with a
non-hit contributing zero to the reciprocal rank.

Every one of these is implemented once, in a module with no dependency beyond numpy and scipy, and
is covered by tests against closed forms and hand-computed values. They were extracted from the
analysis harnesses precisely so that the arithmetic behind the reported intervals could be
exercised by continuous integration rather than assumed.

*Implemented in `core/stats.py`; verified by `tests/test_stats.py`.*

## The ClinVar-blinded protocol, per tool

Several comparators consume ClinVar through the PP5 and BP6 criteria, which would let them reproduce
rather than predict the label on any ClinVar-derived truth surface. The protocol therefore differs
by tool and is stated for each.

- **DISCERN** applies no ClinVar-derived criterion. PS1 and PM5, which require a same-residue
  ClinVar lookup, are implemented but deliberately not wired in for these runs. This is the reason
  they appear zero times in the per-criterion comparison, and it is a protocol choice rather than a
  gap in the engine.
- **REVEL and AlphaMissense** are precomputed scores and consume no classification.
- **GeneBe** does apply PP5 and BP6. Rather than exclude it, it is retained as an exhibit: it
  attains a perfect area under the curve on the eRepo surface and continues to do so under the
  time-split, which demonstrates the circularity concretely.
- **InterVar and BIAS-2015** are positioned by published metrics obtained on different variant sets,
  and are labelled as such wherever they appear.

The time-split panel additionally restricts grading to variants whose expert classification
post-dates the ClinVar snapshot bundled with the comparator tools (2021-05-01), so that no tool
could have memorised the label irrespective of protocol.

## The curated case benchmark and the HPO crosswalk

Cases were curated from published reports, each represented as the causal gene together with the
clinical and laboratory findings the report describes, and each carrying the PubMed identifier of
its source. No article text is reproduced anywhere in the deposit; the benchmark ships as
identifiers, extracted terms, and the expected diagnosis.

**Representability in the Human Phenotype Ontology** was determined by inverting the committed
crosswalk, which maps HPO terms to the engine's internal finding vocabulary, and asking which
findings used anywhere in the benchmark have at least one HPO term. Thirteen of the forty-eight
distinct findings do. The thirty-five that do not are laboratory assay results with no HPO
representation, including the ristocetin mixing study, flow cytometry for CD42 and alphaIIbbeta3,
multimer patterns, light transmission aggregometry, and the prothrombinase assay. Nineteen of the
forty-two cases carry no HPO term at all and cannot be ranked by a phenotype-driven tool.

*Implemented in `eval/phenotype_tool_comparison.py`; verified by `tests/test_diagnosis_baselines.py`.*

## The phenotype-blind baselines

Three reference points, none of which sees a phenotype:

- **Gene lookup.** Ranks the cluster by placing the diseases that list the case's gene first, each
  block ordered by prior. A static table with no likelihood ratios and no findings. Verified to
  produce an identical ranking when every finding is stripped from the case.
- **Prior only.** Ranks the cluster by prior and ignores the gene entirely.
- **Uniform random within cluster.** The floor, reported as the analytic expectation together with a
  simulated interval.

Cases are stratified into three mutually exclusive groups: those whose gene maps to exactly one
disease within its cluster, those carrying no causal gene, and the informative subset whose gene is
present but maps to more than one disease, so that something other than the gene must break the tie.
The three groups partition the benchmark exactly.

*Implemented in `eval/gene_only_baseline.py`.*

## LIRICAL run configuration

LIRICAL version 2.4.1 was run in phenotype-only mode inside an `eclipse-temurin:17-jre` container.
Nothing was installed on the host.

Per case, the observed HPO terms were supplied via `--observed-phenotypes` and, where the case
records an explicitly absent finding with an HPO representation, the corresponding terms were
supplied via `--negated-phenotypes`, so that LIRICAL receives the same pertinent negatives the
engine uses. Genome build was left unset, which selects phenotype-only operation and requires no
VCF. Only the twenty-three cases carrying at least one HPO term were run; the remaining nineteen
cannot be ranked by a phenotype-driven tool and are reported as such rather than scored as failures.

**A practical note for reproduction.** LIRICAL's own downloader retrieves `mim2gene_medgen` over
FTP, which many institutional networks block. The same file is served over HTTPS from the identical
NCBI path, and substituting it allows the remaining resources to download normally.

**Scoring.** LIRICAL ranks diseases genome-wide by OMIM and Orphanet identifier, whereas the engine
ranks within a cluster of three to eight. Scoring is therefore performed at gene and
disease-family resolution: a hit counts if the top-k contains any disease identifier the Human
Phenotype Ontology annotates to the case's causal gene. The gene-to-disease crosswalk is derived
from the committed HPO annotation file rather than asserted, so every identifier traces to a public
source. This is the generous reading for LIRICAL, since it cannot be penalised for differences in
OMIM granularity.

A second, like-for-like arm collapses LIRICAL's genome-wide ranking onto the same cluster members
the engine ranks, so both tools order the same three to eight diseases. That arm is the one the main
text quotes.

**Exomiser** was not run. It requires a VCF, and seeding one per case with zygosity matched to each
report's stated inheritance introduces a correctness risk judged larger than the comparison's
marginal value once LIRICAL had established the phenotype-channel result on the same inputs. This is
recorded as a scope decision rather than presented as a handicapped run.

*Implemented in `eval/lirical_arm.py`.*

## Reproducibility

All reported values are produced by scripts archived with the code record and are locked by
regression tests that read the committed outputs, so a change to any benchmark that moves a reported
number fails the build rather than propagating silently. A further test greps the documentation for
claims retired during pre-submission analysis and fails if any is asserted again without its
retraction.

Third-party databases are referenced by versioned manifest giving source URL, version or snapshot
date, and checksums of the exact files used, rather than redistributed. Software and database
versions are listed in Table S9.
