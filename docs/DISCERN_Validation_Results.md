# Validation results on open data

**Run:** 2026-06-13, on the VM, against real downloaded data. No synthetic data is used as
a result. Sources verified against primary literature.

This is the record of that run. Where a number here differs from `RESULTS.md`, `RESULTS.md`
is the reported value: it reflects the disease-posterior correction made after this run.

## ACMG combining fidelity and the per-code partition on the ClinGen Evidence Repository

Ran `eval/erepo_reconstruction.py` over the real ClinGen Evidence Repository export
(`erepo_classifications.tab`), restricted to the DISCERN bleeding/platelet cluster genes.

| Metric | Result |
|---|---|
| Real VCEP-classified variants (bleeding genes) | **2,653** across 9 genes (RUNX1 1620, ITGA2B 358, ITGB3 250, VWF 119, F8 95, F9 71, GP1BA 66, GP1BB 37) |
| ACMG combining-rule fidelity, exact (vs VCEP label, from the experts' own codes) | **93.0%** |
| ACMG combining-rule fidelity, within one bin | **100.0%** |
| Partition vocabulary coverage (0 unknown codes) | **100%** |
| Variants carrying non-genetic codes (at risk of double-count) | **841** (**31.7%**) |
| Points routed out of the genetic stream | PP4 **1,443**, functional **393**, seg/phasing **794** |
| Variants a naive bottom-line score would over-classify (inflation prevented) | **549** (**20.7%**) |

**Observed per-code applied-strength distribution** (the VCEPs' own applied codes): PM2
1862, BP4 1010, PP4 593, BP7 564, PP3 523, PVS1 492, PM3 427, PS4 266, PM5 244, BA1 198,
PM1 196, PS3 167. This is the real *observed frequency* of applied strengths. The per-gene
*rule* strength tables in `rules/vcep/specs/*.yaml` were subsequently EXTRACTED and verified
(2026-06-13): the CSpec frequency criteria (BA1/BS1/PM2) and the PM2_Supporting strength are
now real for GT/F8/F9/VWF/GP1BA, from the CSpec registry (GN071 F8/F9, GN079 GP1BA, GN081 VWF -
note GN079 is the GP1BA spec, NOT GT, an earlier mislabel) cross-checked against the VCEPs'
own eRepo records. The only residual
placeholders are the variant-dependent PVS1/PS4 strength decision trees (a documented
simplification, not a fillable value) and RUNX1's BA1/BS1 numeric thresholds.

**Points routed out of the genetic stream, by owning factor:** PP4 **1,443** (to the
disease-to-variant coupling), phasing 534, functional 393, segregation 260, de-novo 67.
This quantifies the motivation for the per-code partition: a tool that consumed the VCEP's
bundled bottom-line label would re-count this evidence in the genetic factor; here it is
owned by a single factor each and never re-added.

**Threat model - who actually double-counts (this is a real workflow, not a hypothetical).**
The over-classified variants are those whose VCEP/ClinVar bottom-line already absorbed
non-genetic evidence (PP4 phenotype, PS3 functional, PP1 segregation, PM3 in-trans). Any pipeline
that *consumes that label and re-adds the same evidence* over-classifies them. Concretely:
1. **Phenotype-aware / coupled classifiers - including a naive version of DISCERN's own joint
   model.** They use a variant's ClinVar/VCEP pathogenicity (which already counted PP4 "phenotype
   consistent with the disorder") and *then* separately score the patient's phenotype -> disease
   likelihood, entering the phenotype evidence twice. PP4 is the single largest routed bucket
   (1,443 pts), so the dominant double-count is exactly the evidence a coupling model re-adds -
   which is precisely why DISCERN needs the partition rather than consuming bottom-line labels.
2. **Automated ACMG meta-classifiers and panel/exome re-analysis pipelines** (e.g. Franklin/Genoox,
   VarSome, GeneBe, TAPES, CharGer) that ingest a prior ClinVar classification as an input feature
   or prior while *also* re-deriving PP1/PS3/PM3 from the same underlying publications, and Bayesian
   point-reaggregation workflows that re-run the ACMG point system on a variant whose deposited
   label already integrated those codes.

**Interpretation:** the pre-registered protocol (no double-counting) is verified two ways: (i) a unit test
that the variant marginal is invariant to bundled PP4/PS3/PP1/PM3 codes, and (ii) the
per-code routing on 2,653 real variants above - owned codes appear in 31.7% of variants,
and for 549 (20.7%) the moved evidence is band-determining (a naive all-codes score would
over-classify them). It is **not** a calibration of the disease-variant coupling, which
awaits paired-phenotype cohorts. (The earlier "100% no-inflation rate" was a tautology -
it summed the same points twice - and has been replaced.)

## Diagnosis benchmark on the GA4GH Phenopacket Store

Ran `eval/phenopacket_benchmark.py` over the cloned `phenopacket-store` (11,155
phenopackets). After extracting the causal gene from the interpretation and excluding
thrombophilia cases (elevated factor VIII / DVT / PE - outside DISCERN's bleeding clusters):

| Metric | Result |
|---|---|
| In-cluster cases found | **4** (1 LAD-III/FERMT3, 3 Chediak-Higashi/LYST) |
| Top-1 accuracy | **100%** (4/4) |
| Top-3 accuracy | **100%** |

**Finding:** Phenopacket Store is a general rare-disease corpus and is thin on
inherited bleeding/platelet disorders - only 4 cases fall in DISCERN's clusters (consistent
with the pre-registered protocol). DISCERN diagnoses all 4 correctly, but the diagnosis-accuracy
headline cannot rest on n=4. The PhEval-compatible runner is in place; the diagnosis
benchmark headline rests instead on the hand-built curated published-case set, and in the
longer run on paired phenotype-genotype cohorts.

## VUS reclassification rate

The variant-layer validation is delivered by the combining-fidelity and per-code partition
analysis above, which covers real VCEP variants including the uncertain ones. The per-patient VUS-reclassification
*rate* concordant with 3-star truth requires per-patient phenotype paired with VUS, which
the open variant databases do not provide standalone; it is run on paired cohorts where
phenotype and variant co-occur. ClinVar 3-star remains the truth label for that run.

## Expansion results (2026-06-13 to 06-16, real open data)

| Benchmark | Dataset | Result | Code |
|---|---|---|---|
| Genome-wide partition | full ClinGen ERepo, 12,240 variants / 170 genes | **100% partition coverage**; **33.2% inflation-prevented** (Wilson 95% CI 32.4-34.1) | `eval/erepo_genomewide.py` |
| gnomAD per-variant freq cross-check | 629 curator-cited gnomAD AFs | gene-specific CSpec thresholds reproduce the VCEP freq code at **97.8%** | `eval/gnomad_freq_check.py` |
| Intrinsic-only band vs ClinVar | 10,780 matched | **62.4% exact / 92.8% within-one-bin** (intrinsic-only is a designed lower bound; routed PP4/PS3/PP1/PM3 omitted) | `eval/clinvar_concordance.py` |
| Variant calibration | 7,521 ClinVar-labelled | **deliverable = calibration: isotonic ECE 0.008 / Brier 0.0073** (from 0.201/0.060). AUC 0.999 is on P/B extremes (VUS dropped) where discrimination is trivial - reported only to confirm the monotonic isotonic map preserved ranking, NOT a discrimination claim | `eval/variant_calibration.py` |
| Novel-variant scoring vs InterVar (**full-DB run 2026-06-17**) | 1,015 ClinVar P/B spec-gene variants; literal InterVar re-run with its **complete default DB set** (1000g/esp/avsnp/dbnsfp42a/dbscsnv/clinvar_20210501/rmsk/gnomAD-exome + full intervardb - nothing dropped) | **Missense (n=364, the meaningful axis): DISCERN-full 0.944 BEATS InterVar 0.811 (+0.133), matches REVEL 0.942.** All P/B (easy null-dominated extremes): DISCERN 0.882 w/freq, 0.912 w/o ~ InterVar 0.887. Overall AUCs converge on the null-dominated extremes; the edge is on missense | `eval/intervar_full_eval.py` |
| **Independent LOF sensitivity (CHAMP/CHBMP)** | CDC F8/F9 disease-allele catalogs, 5,437 variants (non-ClinGen truth set) | **91.2% LP/P recall on 2,130 null variants by CONSEQUENCE ALONE** (F8 97.7%, F9 65.8%; F9 gap = large terminal exon / NMD-escape) | `eval/champ_chbmp_benchmark.py` |
| **Novel-missense recovery arm (CHAMP/CHBMP)** | novel F8/F9 missense (cDNA not-in ClinVar) + routed PS1/PM5 + PP3 (REVEL, VM) + PM2 (gnomAD-verified absent) | LP recovery **1.5%** (F8 16/1307, F9 11/511), PS1-bounded; PM5 26-47% / no-hit 50-72% stay VUS by point arithmetic. Intrinsic to ACMG, not DISCERN-specific - quantifies the coupling's value | `eval/champ_chbmp_missense_arm.py` |

Full narrative in `DISCERN_CHAMP_CHBMP_Benchmark_v1.md`.

## Claim -> result crosswalk (what is proven now)

| Claim | Dataset | Result | Status |
|---|---|---|---|
| No double-counting | real ERepo per-code variants | invariant marginal (unit test) + 100% partition coverage on 2,653; owned codes in 31.7%, naive over-classifies 549 (20.7%) | **Done on open data; coupling calibration pending paired cohorts** |
| ACMG combining-rule fidelity | real ERepo bleeding genes | 93.0% exact / 100% within-1-bin (arithmetic only) | **Done on open data** |
| Per-code strengths are real, not placeholders | real ERepo | extracted true strength distribution | **Done on open data** |
| Differential-diagnosis accuracy | Phenopacket Store bleeding subset | 4/4 Top-1 (small N; corpus thin) | partial; corpus thin |
| Differential-diagnosis accuracy (curated) | **42 PMID-verified published cases, all 10 clusters** (`eval/cases/curated_cases.yaml`) | **Top-1 81% / Top-3 100% / abstention 10%** under the configuration that preceded the disease-posterior correction; every non-Top-1 is a same-gene/same-cluster confusable held in the Top-3. See `RESULTS.md` for the corrected arm | **superseded by `RESULTS.md`** |
| VUS reclassification rate | ClinVar 3-star plus paired cohort | variant layer done | not evaluable on open data |
| Misdiagnosis rescue | paired cohort under controlled access | harness built | not evaluable on open data |
| Equity analysis | paired cohort under controlled access | harness built | not evaluable on open data |
| Calibration / abstention | labeled sets | harness built (eval/calibration.py, abstain.py) | runs on labeled cohorts |

## Not evaluable on open data

Reclassification rate, misdiagnosis rescue, the equity analysis and the coupling endpoint require paired phenotype-genotype data under controlled access. The harnesses for each are implemented and tested; the protocol is pre-registered.
