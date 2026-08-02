# DISCERN - Phase R pre-submission results

Phase R exists to find the answers to five reviewer objections before a reviewer does. Every result
below is reported as it came out. Two of them cost the paper a claim, one of them cost it a
headline number, and one exposed a modelling gap that had to be fixed. All four are recorded here
with the corrected framing, because a stated weakness is survivable and a discovered one is not.

All numbers regenerate from a clean clone at the current commit. The harnesses are
`bench/phase_r_variant.py`, `bench/phase_r_gene_term_sensitivity.py`, `eval/gene_only_baseline.py`,
`eval/phenotype_tool_comparison.py`, `eval/lirical_arm.py` and `eval/erepo_genomewide.py`, and the
headline values are locked by `tests/test_phase_r.py`.

---

## Disposition table

| ID | Check | Gate | Result | Favourable | Manuscript action |
|---|---|---|---|---|---|
| R1 | Out-of-sample calibration | G-R1 | Already out-of-fold; now with CIs and a leakage assertion | Yes | State the protocol in Methods; add CIs |
| R2 | Calibrated comparators | G-R2 | Post-hoc calibration closes most of the gap; CIs overlap | **No** | **Retire "only calibrated tool"**; reframe on the classification system |
| R3 | AUROC CIs and paired test | G-R3 | Gap to REVEL is not significant (p=0.29); added-value analysis added | Mixed | "Tracks REVEL" now defensible; report the intrinsic-only ceiling |
| R4 | Gene-only baseline | G-R4 | Baseline 93%; exposed a missing model term, now fixed | **No** (as a claim) | Diagnosis arm cannot carry an accuracy claim; report the confound |
| R5 | LIRICAL on the curated cases | G-R5 | Run; DISCERN 100% vs 57% within-cluster on n=23 | Yes, weakly | Report with the HPO-coverage caveat |
| R6 | Reproduce the ~33% figure | G-R6 | Reproduces exactly with CI | Yes | Keep, with an explicit denominator |
| R7 | Partition novelty framing | - | n/a | n/a | Intro and Discussion edit |

---

## R1 - calibration is estimated out-of-sample

The audit found the protocol was already correct: isotonic regression is fit inside each training
fold of a 5-fold stratified split and the metric is accumulated only on held-out folds. Phase R
makes that explicit rather than implicit, asserts fold disjointness in code, and attaches bootstrap
confidence intervals (1,000 resamples, percentile method).

| Surface | n | ECE | 95% CI | Brier |
|---|---|---|---|---|
| eRepo-primary | 425 | 0.0168 | 0.0105 - 0.0466 | 0.0635 |
| time-split (approved after 2021-05-01) | 383 | 0.0186 | 0.0117 - 0.0485 | 0.0646 |

**Disposition: favourable.** The headline calibration number survives, and it is now reported as a
generalisation estimate with an interval rather than a point.

## R2 - the comparators, calibrated on the identical folds

REVEL and AlphaMissense were put through the same isotonic protocol and the same folds.

| Tool (eRepo-primary, n=425) | ECE | 95% CI | Brier |
|---|---|---|---|
| DISCERN, isotonic out-of-fold | 0.0168 | 0.0105 - 0.0466 | 0.0635 |
| AlphaMissense, isotonic out-of-fold | 0.0288 | 0.0239 - 0.0652 | 0.0939 |
| REVEL, isotonic out-of-fold | 0.0429 | 0.0298 - 0.0732 | 0.0732 |
| REVEL, raw score | 0.0919 | 0.0791 - 0.1215 | 0.0810 |
| AlphaMissense, raw score | 0.1674 | 0.1345 - 0.1996 | 0.1383 |

The time-split surface behaves the same way (DISCERN 0.0186, AlphaMissense 0.0310, REVEL 0.0339).

**Disposition: unfavourable, and the claim must change.** DISCERN is best on every point estimate
and on Brier, but the confidence intervals overlap the calibrated comparators. A raw metapredictor
score *can* be post-hoc calibrated, and once it is, most of the advantage disappears.

**Retired claim.** "DISCERN is the only tool emitting a calibrated probability."

**Replacement claim.** Calibration is not a property REVEL or AlphaMissense ship: obtaining it
required fitting an isotonic layer on labelled expert data, which is a modelling step the tools do
not provide and which itself consumes the labels a user is trying to predict. DISCERN is calibrated
as delivered. More importantly, a calibrated score is still only a score - it cannot emit an ACMG
classification, an auditable criterion trail, a partition guarantee, an abstention decision, or a
treatment-safety hard stop. The differentiator is the calibrated, auditable *classification
system*, not calibration alone.

## R3 - the AUROC gap, and what the pipeline adds

| Surface | DISCERN | REVEL | AlphaMissense |
|---|---|---|---|
| eRepo-primary (n=425) | 0.9386 (0.9064-0.9665) | 0.9536 (0.9315-0.9717) | 0.9315 (0.9058-0.9531) |
| time-split (n=383) | 0.9272 (0.8825-0.9606) | 0.9511 (0.9240-0.9728) | 0.9292 (0.9012-0.9547) |

Paired on the same variants:

| Surface | delta (DISCERN - REVEL) | DeLong z | p | paired bootstrap 95% CI |
|---|---|---|---|---|
| eRepo-primary | -0.0150 | -1.064 | 0.287 | -0.0435 to +0.0118 |
| time-split | -0.0239 | -1.425 | 0.154 | -0.0598 to +0.0076 |

**The gap is not statistically significant on either surface**, and the interval spans zero. REVEL
is DISCERN's own PP3 input, so "tracks REVEL" was always the right verb; it is now supported rather
than asserted.

### The intrinsic-only ceiling - the most consequential finding in Phase R

Measuring what the pipeline adds turned up something the manuscript did not previously state. On
the eRepo missense surface, using intrinsic sequence evidence alone, **DISCERN reaches a maximum of
3.0 ACMG points against the 6 required for Likely Pathogenic, and assigns a pathogenic band to zero
of 425 variants.**

| Decision rule (eRepo-primary, n=425) | coverage | accuracy on resolved | pathogenic recall | benign recall |
|---|---|---|---|---|
| DISCERN ACMG class | 18.8% | 0.938 | **0.000** | 0.688 |
| DISCERN calibrated probability, coverage-matched | 16.2% | 0.913 | 0.098 | 0.294 |
| REVEL at ClinGen thresholds | 80.7% | 0.956 | 0.851 | 0.541 |
| REVEL, coverage-matched to DISCERN | 18.6% | 0.975 | 0.123 | 0.349 |

Every band DISCERN resolves on this surface is benign-side. This is not a defect: PM2 plus PP3
cannot reach 6 points, and the criteria that carry missense pathogenicity - PS3 functional, PS4
case-control, PM1 hotspot, PM5, PP4 phenotype - are owned by the functional, disease and coupling
factors and are not derivable from sequence. It is the partition behaving exactly as specified.

It is also the clearest argument the paper has for the coupling. The evidence the intrinsic surface
is missing is precisely the evidence a disease model supplies. On the 830 missense variants the
expert panel itself left uncertain, DISCERN assigns a band to 188 and a ClinGen-threshold REVEL rule
to 518 - but neither number is an accuracy, because eRepo VUS carry no independent truth. That is
the gap the pre-registered coupling endpoint exists to close.

**Disposition: mixed.** The discrimination claim is now defensible with a CI. A REVEL threshold
rule makes more calls than the ACMG framework permits, and at matched coverage its calls are at
least as accurate. Both facts are stated.

## R4 - the gene-only baseline, and the modelling gap it exposed

A naive gene-to-most-common-disease lookup, with no phenotype and no likelihood ratios, scores
**Top-1 93%** on the same 42 cases. That alone disqualifies the pooled Top-1 as evidence for the
phenotype engine.

Building the baseline also exposed a real gap. The joint factorised as
`P(E_pheno|D) x P(E_geno|V) x P(E_func|D,V) x P(V|D) x P(D)`, and `P(V|D)` normalises over the five
variant states - so the gene the variant was sequenced in cancelled out on marginalisation. A
variant in F8 argued no more for haemophilia A than for haemophilia B. The factorisation was
missing `P(G|D)`, and eight of the 42 cases failed for exactly this reason.

Adding the missing term (`jointdx/factorgraph.gene_loglik`) completes the factorisation. Each
stream still enters once: `P(V|G,D)` is the variant's state given its location, `P(G|D)` is the
location itself, and no other factor claims the gene. The partition invariance test is unaffected.

Two constraints govern the term, both architectural rather than fitted:

- **It must not veto a treatment-divergent competitor.** Finding a variant in one gene does not
  rule out a disease of another gene, which may never have been sequenced. The safety interlock is
  therefore evaluated on a gene-blind posterior against the engine's own leading call, so a
  platelet-type von Willebrand case still hard-stops desmopressin for type 2B.
- **It must stay weaker than a cluster's deciding observation**, or no laboratory result could ever
  overturn the gene and the value-of-information layer would be decorative.

The sweep in `bench/phase_r_gene_term_sensitivity.py` reports the whole grid:

| gene likelihood ratio | curated Top-1 | shared-gene | unique-gene | hard-stop sens/spec | deciding assay still wins |
|---|---|---|---|---|---|
| 1.0 (gene inert, prior behaviour) | 81% | 100% | 78% | 1.00 / 1.00 | yes |
| 1.5 | 88% | 100% | 86% | 1.00 / 1.00 | yes |
| 2.33 | 90% | 100% | 89% | 1.00 / 1.00 | yes |
| **4.0 (committed)** | **100%** | **100%** | **100%** | **1.00 / 1.00** | **yes** |
| 9.0 | 100% | 100% | 100% | 1.00 / 1.00 | no |
| 19.0 | 100% | 100% | 100% | 1.00 / 1.00 | no |

The committed value is the largest at which the deciding assay still overturns the gene. Safety is
invariant across the entire grid. The result is flat above 4.0, so it does not sit on a knife edge.

With the term in place:

| Stratum | n | DISCERN | gene lookup | prior only | random floor | delta (95% CI) | McNemar p |
|---|---|---|---|---|---|---|---|
| pooled | 42 | 100% | 93% | 21% | 32% | +7% (0% to +17%) | 0.25 |
| unique-gene | 36 | 100% | 100% | 17% | 32% | 0% | 1.00 |
| shared-gene | 3 | 100% | 33% | 33% | 33% | +67% (0% to +100%) | 0.50 |
| no gene given | 3 | 100% | 67% | 67% | 31% | +33% (0% to +100%) | 1.00 |

**Disposition: unfavourable as a claim, and the reframe stands.** Top-1 100% on 42 cases is not
evidence the phenotype engine works, because a static lookup table gets 93% and the difference is
not significant. The informative subset - cases where the gene does not settle the answer - is
**n=3**. That is the honest finding, and it is the same structural point as the coupling negative:
public curated cases supply the gene and almost no phenotype depth, so there is nothing left for a
phenotype engine to do. The diagnosis arm is presented as a scoped component resting on abstention,
zero confident errors, and the safety interlock, not on Top-1 accuracy.

Downstream effects of the fix, all reported: curated Top-1 81% to 100%, diagnosis ECE 0.1412 to
0.1180, abstention 10% to 7%, confidently-wrong count unchanged at 0, safety unchanged at 100%
sensitivity and 100% specificity. The risk-coverage curve is now uninformative - accuracy is 100%
at full coverage, so there is no headroom for abstention to demonstrate - and the previous claim
that "accuracy rises to 100% on the most-confident half" is withdrawn.

## R5 - LIRICAL on the curated cases

Run rather than declined. LIRICAL 2.4.1, phenotype-only mode, in a container on the VM, with both
observed and negated HPO terms supplied so it sees the same pertinent negatives DISCERN uses.

The plan assumed sparse HPO handicaps every tool equally. It does not, and the measurement says why:
**only 13 of the benchmark's 48 discriminating features have any HPO representation.** The 35 that
do not are the laboratory assays that actually separate these diseases - RIPA and its mixing study,
flow cytometry for CD42 and alphaIIbbeta3, multimer patterns, light transmission aggregometry, the
prothrombinase assay. Nineteen of 42 cases therefore carry no HPO term at all and LIRICAL cannot
rank them; median terms per case is 1, mean 1.5.

On the 23 cases it can accept:

| Arm | n | Recall@1 | Recall@3 | Recall@5 | MRR |
|---|---|---|---|---|---|
| LIRICAL, genome-wide | 23 | 13% | 26% | 35% | 0.215 |
| LIRICAL, restricted to the cluster | 23 | 57% | 96% | 96% | 0.754 |
| DISCERN, same 23 cases | 23 | 100% | 100% | 100% | 1.000 |

**Only the restricted row is a fair comparison** - it is the one where both tools order the same
three to eight diseases. The genome-wide row is not: LIRICAL ranks roughly 8,600 diseases without
being given the gene, while DISCERN ranks within one cluster and is given it.

**Disposition: favourable but weakly.** DISCERN is ahead on the like-for-like arm (23/23 against
13/23), and the reason LIRICAL trails is not that it is a weak ranker but that the discriminating
channel in this domain has no HPO encoding. That is a fact about the domain, and it is the same
finding as R4 and as the coupling negative, arriving from a third direction.

## R6 - the ~33% inflation figure reproduces

Re-run at the current commit over the full ClinGen eRepo (`eval/erepo_genomewide.py`):

- 12,240 variants across 170 genes, 38,802 applied evidence codes
- partition vocabulary coverage 100.0000%, zero uncovered codes (Gate G9)
- variants carrying a non-genetic owned code: 6,357 (51.9%)
- **inflation prevented: 4,068 variants, 33.2%, 95% CI 32.4% - 34.1%** (Wilson)

**Denominator, stated precisely:** 33.2% of *all* eRepo variants carrying at least one applied code
would be placed in a higher (more pathogenic) ACMG band by a naive scorer that sums every applied
criterion into a single factor, compared with the band produced by the variant-intrinsic factor
alone. Among the at-risk subset - the 6,357 variants that actually carry a non-genetic owned code -
the rate is 64.0%.

**Disposition: favourable.** The figure is kept, with the denominator and mechanism stated, and is
locked by a test.

## R7 - positioning the partition against prior art

Framing only, no computation. The Introduction and Discussion now separate three things that a
reviewer citing ClinGen SVI guidance would otherwise collapse:

1. That double-counting is undesirable is known and long-stated. That is not the contribution.
2. The contribution is the **formal partition** that makes exclusive ownership structural rather
   than advisory inside a coupled disease-variant model, where the temptation to re-add PP4 or PS3
   on top of an already-absorbing bottom line is created by the coupling itself.
3. It is **empirically demonstrated against expert-panel classifications** (per-criterion kappa,
   with PS3/PS4/PM1/PM5 applied zero times on the intrinsic surface) and **quantified** (R6).

This is distinct from the PP5/ClinVar circularity problem, which BIAS-2015 itself flags and which
DISCERN also demonstrates through the GeneBe exhibit. The mechanisms are different: circularity is
re-reading a conclusion, the partition prevents re-using the evidence behind one.

---

## What Phase R changed in the repository

- `jointdx/factorgraph.py` - added `P(G|D)`; separated the incidental-variant constant from the
  gene term; `joint(..., gene_evidence=False)` for the safety view
- `jointdx/orchestrate.py` - safety evaluated gene-blind against the engine's own leading call
- `safety/interlock.py` - `flags(..., lead_id=...)` so flags always report against the shown call
- `bench/phase_r_variant.py`, `bench/phase_r_gene_term_sensitivity.py` - R1, R2, R3, and the sweep
- `eval/gene_only_baseline.py` - R4
- `eval/phenotype_tool_comparison.py`, `eval/lirical_arm.py` - R5
- `tests/test_phase_r.py` - 10 guards locking every headline above; suite is 178 tests

## What still needs the author

- OSF time-stamp for the pre-registration (Gate G12)
- A tagged release to mint the Zenodo code DOI, then the supporting-data DOI
- The randomized reader-study arms
- Paired-cohort access for the coupling endpoint H6 (Gate G13)
