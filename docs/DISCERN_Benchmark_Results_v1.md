# DISCERN - Benchmark Results v1 (Tracks 1 and 3)

**DISCERN = Diagnostic Inference from Shared-mechanism Coupling of Evidence in Rare Nosology.**

**Date:** 2026-06-19 - **Author:** Anees Ahmed
Mahaboob Ali (`ahmedanees-m`). **Scope:** Paper 1 positioning on data in hand. The coupling proof and
prospective clinical validation remain Paper 2 (cohort-gated) and are out of scope here.

**What ran this round (open data, committed + reproducible):** Track 1 (variant arm vs the current
ACMG/predictor tool set) and Track 3 (the trustworthiness layer). Tracks 2, 4 and 5 are scoped with
their requirements in the last section. All reference DOIs/PMIDs in the plan were re-verified against
Crossref/NCBI on 2026-06-19 (all six resolve correctly, including the 2026 Bioinformatics ACMG-tools
evaluation and DeepRare).

---

## TRACK 1 - Variant arm head-to-head (missense / VUS-adjacent axis)

**Set:** the h4set variant panel (GRCh38), missense subset n=1144; **P/B truth (ClinVar): n=342**
(241 pathogenic/likely-pathogenic, 101 benign/likely-benign); 802 ClinVar-VUS held out of the
discrimination metric. **Annotation source:** GeneBe public API (snapshot 2026-06-19, cached at
`bench/data/genebe_h4set.jsonl`), which supplies gnomAD AF + REVEL + AlphaMissense + consequence, so
DISCERN is scored locally (`rules.variant_scoring.score_variant`) on the identical variants with NO
ClinVar-derived codes (the ClinVar-blinded protocol).

| Tool | AUROC | AUPRC | sens@90%spec | Calibrated probability? |
|---|---|---|---|---|
| **DISCERN** (intrinsic, af-on) | **0.935** | 0.954 | 0.909 | **Yes - isotonic ECE 0.017** |
| DISCERN (frequency-blind) | 0.934 | 0.953 | 0.842 | yes |
| REVEL (the PP3 input) | 0.948 | 0.976 | 0.888 | no (raw score; ECE 0.091) |
| AlphaMissense | 0.921 | 0.967 | 0.772 | no (raw score; ECE 0.161) |
| GeneBe (acmg_score) | 1.000* | 1.000* | 1.000* | no (class + points only) |
| InterVar (full-DB, prior) | 0.811 | - | - | no (class only) |

`*` GeneBe's perfect score is **ClinVar circularity, not skill** (see below) - it is not a fair number.

**ACMG concordance (categorical; abstention is a deliberate non-call, not an error):**
- **DISCERN:** abstention 79% on missense; of the calls it does make it resolves only to the benign
  side (73 calls, 94.5% correct) and makes **zero pathogenic calls on intrinsic-only evidence** - the
  designed conservative lower bound (a missense needs PP4/PS3/PP1 to reach LP, which is exactly the
  coupling motivation). It does not over-call.
- **GeneBe:** abstention 1%, accuracy-on-resolved 1.000, MCC 1.000 - because it reproduces ClinVar.

**The ClinVar-circularity finding (plan section 0, guard #2 - made empirical):**
GeneBe's `acmg_score` is perfectly separable on this ClinVar P/B set (benign max -1, pathogenic min
+4) and its class matches the ClinVar direction on **337/337 resolved calls**. The named ClinVar codes
(PP5/BP6) appear on only 12% of calls, yet the AUROC stays 1.000 even on the PP5/BP6-name-blind subset
(n=302) - so the dependence is systemic, not just the named codes. **A ClinVar-consuming classifier
therefore cannot be fairly graded on a ClinVar-derived truth set.** DISCERN, REVEL and AlphaMissense
apply no ClinVar-derived evidence, so their numbers are not inflated this way. This is itself a
publishable methodological result and motivates the eRepo-primary / time-split run (next section).

**Reading + Gate G-T1.** On raw discrimination DISCERN (0.935) tracks REVEL (0.948) - expected,
since DISCERN ingests REVEL as its PP3 signal - and clears the legacy InterVar anchor (0.811) and
AlphaMissense (0.921). DISCERN does NOT claim a discrimination win over its own metapredictor input.
The differentiator, stated as the gate requires, is that **DISCERN is the only tool in the set that
emits a calibrated pathogenicity probability** (ECE 0.017 vs REVEL 0.091 / AlphaMissense 0.161;
GeneBe and InterVar emit a class only), and the only one whose conservatism (abstain rather than
over-call) is principled. **Gate G-T1: met via the explicit calibration differentiator + the stated
surface (discrimination ~ REVEL, by construction).**

Artifacts: `bench/track1_variant_headtohead.csv`, `bench/track1_metrics.json`,
`bench/track1_variant_headtohead.py`, `bench/genebe_client.py`.

---

## TRACK 1' - eRepo-primary + time-split (removes the ClinVar-circularity caveat) [Phase A2 Part A]

The Track-1 ClinVar surface let GeneBe score a circular AUROC 1.0. This re-run grades on the
**FDA-recognized ClinGen eRepo expert-panel** classifications instead, on two surfaces: **eRepo-primary**
(all cluster-gene P/B) and a **time-split** panel (eRepo variants expert-approved after 2021-05-01, so
InterVar's `clinvar_20210501` bundle could not have memorised the label). Truth = eRepo assertion;
set = 2,239 eRepo SNVs in the cluster genes (GRCh38 coords taken from the eRepo HGVS column - no new
download). DISCERN scored locally from the GeneBe annotation (`bench/genebe_erepo.py`); the whole
re-run is disk-free.

| Tool | eRepo-primary AUROC (missense, n=425) | time-split AUROC (n=383) | ECE (eRepo / time-split) |
|---|---|---|---|
| **DISCERN** | **0.939** | **0.927** | **0.017 / 0.019** |
| REVEL | 0.954 | 0.951 | 0.092 / 0.091 |
| AlphaMissense | 0.931 | 0.929 | 0.167 / 0.180 |
| GeneBe (exhibit) | 1.000 (circular) | 1.000 (circular) | n/a (class-only) |

**The calibration differentiator holds on the independent surface and on the time-split** (DISCERN ECE
0.017-0.019 vs REVEL ~0.09 / AlphaMissense ~0.17), and DISCERN's discrimination tracks REVEL exactly as
before. GeneBe stays at 1.000 even on the time-split - it queries live ClinVar, so it cannot be
time-split and is retained only as the circularity exhibit (a ClinVar-consuming tool cannot be fairly
graded on a ClinVar-derived truth). **Gate G-T1' met** for the DISCERN/predictor comparison; the
ClinVar-circularity caveat on the Track-1 claim is superseded by this time-split panel.

**Per-ACMG-code Cohen's kappa vs eRepo (the partition made visible on an expert surface):** DISCERN
agrees with the experts on the codes it derives intrinsically - **PVS1 kappa 0.81, PP3 kappa 0.93**, BA1
0.64 - and applies **none** of the evidence-stream codes it is not given the inputs for (PS3 functional,
PS4 case-control, PM1 hotspot, PM5 = 0 applied), which is the no-double-counting partition: DISCERN does
not re-derive evidence it cannot independently compute. (Divergences: DISCERN applies
PM2_Supporting broadly - 479 vs eRepo's 8 - and BP4 less often than the experts; both are application-
style differences, noted not hidden.)

**InterVar on the eRepo set - disk boundary (deferred).** A *fair* InterVar needs dbNSFP42a
for its PP3/BP4 in-silico evidence; the ANNOVAR hg38 dbNSFP42a database is ~250 GB uncompressed, which
exceeds the VM's ~100 GB free space. Running InterVar without it would be the *handicapped* configuration
the earlier H4 correction explicitly flagged as overstating DISCERN's edge, so it is not run that way.
The legacy InterVar anchor is therefore the prior **full-DB ClinVar-set missense AUROC 0.811**
(labeled as the ClinVar set, per integration item 4); a fair eRepo InterVar run is deferred until disk
allows. The eRepo-primary headline (calibration edge, circularity removal) does not depend on InterVar.

**BIAS-2015 with PP5/BP6 = 0 - setup completed; classifier deferred.** The full
annotation pipeline was stood up on the VM: .NET 6.0.36 + Nirvana 3.18.1 installed, the GRCh38 data
sources (~53 GB) downloaded, the 2,238-variant eRepo VCF (`bench/erepo_to_vcf.py`) annotated to
`bench/data/erepo_nirvana.json.gz`, and the PP5/BP6 ClinVar-blinding implemented (a one-line early
return in `get_bp6`/`get_pp5`). **BIAS-2015 itself could not be run** because its auxiliary-data
preprocessing is broken in the released tags (v2.1.1 and v3.0.0): `src/preprocessing/__init__.py` is
missing while `preprocessing.py` imports bare function names from the package, the
`generate_submitter_counts` import does not match its module (`generate_clinvar_submitter_counts`), and
even patched it requires un-shipped ClinVar `variant_summary`/`submission_summary` + UniProt/AVADA inputs
and a multi-step run. Per the run-spec's explicit fallback (BIAS-2015 must not block the eRepo re-run),
it is deferred, and the **published BIAS-2015 eRepo numbers are cited as the external classifier
reference**: on the FDA-recognized eRepo set BIAS-2015 reports pathogenic sensitivity 73.99% vs InterVar
64.31% and benign sensitivity 80.23% vs 53.91% (Eisenhart 2025, *Genome Medicine*, DOI
10.1186/s13073-025-01581-y) - i.e. it is a class-only ACMG tool that beats InterVar but, like InterVar,
emits no calibrated probability. The 53 GB of Nirvana data was reclaimed (disk back to 102 GB);
re-running BIAS-2015 needs only a ~10-minute Nirvana re-download once its preprocessing is fixed.
Artifacts: `bench/erepo_extract.py`, `bench/genebe_erepo.py`, `bench/erepo_to_vcf.py`,
`bench/track1b_erepo_headtohead.py`, `bench/track1b_erepo_metrics.json`, `bench/data/erepo_bleeding.{tsv,vcf}`,
`bench/data/erepo_nirvana.json.gz`.

---

## Diagnosis arm (Phase A2 Part B) - disposition (not run; deferred to the cohort)

The spec's Part B (LIRICAL/Exomiser head-to-head) has a primary gate of **within-cluster** discrimination
- DISCERN's home turf. That comparison **cannot be run fairly on public data**, for two independent reasons
surfaced during setup (both documented rather than worked around):

1. **The curated case set lacks HPO depth for HPO-driven tools.** Its discriminators are DISCERN's
   discrete clinical + functional features; mapped to HPO they give a mean of **1.3 terms/case, with 20 of
   42 cases having zero** clinical HPO terms (their deciding features are lab assays - RIPA, flow CD42,
   PS-exposure - that have no HPO equivalent). Feeding LIRICAL/Exomiser this would be an unfair, starved
   input, not a real contest.
2. **The public phenopacket corpus lacks confusable-pair density.** The Phenopacket Store bleeding subset
   has rich HPO (mean 5.0 terms) and real GRCh38 variants, but its genes are ACTN1 x14 / F8 x8 / WAS x4 /
   LYST x3 / MPL x1 / FERMT3 x1 - there are no *cluster-mates* present to discriminate against (F8
   haemophilia with no other coagulation-factor cases; LYST/Chediak with no HPS/GPS cases). The
   within-cluster task is undefined on it.

This is the **same public-corpus thinness the Coupling PoC already quantified**, and it is reported here as
a **motivating negative**: the within-cluster diagnosis head-to-head is a genuine reason the confirmatory
work needs the richly-phenotyped, confusable-dense **cohort (Paper 2)**, not a gap in DISCERN. A
genome-wide Recall@k on the 31 phenopackets was considered but declined: genome-wide ranking is Exomiser's
design point, not DISCERN's (DISCERN is a specialist that abstains outside its clusters), so it would
neither test DISCERN's differentiators (within-cluster, calibration, safety) nor be a like-for-like
contest. The diagnosis arm's standing public evidence remains the curated within-cluster proof
(Top-1 81% / Top-3 100%, n=42) plus the trustworthiness results (Track 3); the external head-to-head is
explicitly a cohort deliverable.

---

## TRACK 3 - Trustworthiness layer (calibration / safety / abstention)

### 3a Calibration
- **Variant posterior:** DISCERN isotonic ECE **0.017** vs REVEL 0.091 vs AlphaMissense 0.161 (n=342).
  DISCERN is the calibrated one; the predictor scores are not gene-calibrated probabilities.
- **Diagnosis posterior (n=42 curated cases):** ECE 0.141 (modest, limited by n=42),
  accuracy 81%, and critically **0 confidently-wrong calls** (no diagnosis made with confidence >=0.8
  was incorrect) - the property that matters for a safety-critical aid.

### 3b Safety interlock (n=5 treatment-divergent scenarios: DDAVP/2B, splenectomy/BSS,
related-donor-transplant/RUNX1-FPD, platelet-transfusion/Quebec, rFVIII-monotherapy/2N)
- **Hard-stop sensitivity = 100%** (5/5): the interlock fires whenever the planned treatment is
  contraindicated by a non-excluded competitor.
- **Hard-stop specificity = 100%** (5/5): the interlock stays silent when the planned treatment is
  harmless (tranexamic acid) - it does not cry wolf. Both directions measured, not just sensitivity.

### 3c Abstention / selective prediction (risk-coverage, n=42)
- The confidence ranking is informative: accuracy rises from **81% at full coverage to 100% on the
  most-confident 50%** (monotone payoff confirmed) - the operational value of calibration.
- DISCERN's current `decide()` operating point: coverage 93%, accuracy-on-decided 79%, abstention 7%.
  Note: this threshold abstains on ~3 low-confidence-but-correct cases rather than on the
  errors, so on-decided accuracy (79%) sits just below full (81%). The signal is good (the curve is
  monotone); the threshold is conservatively tuned on n=42 and is a calibration-tuning item, not a
  signal failure.

**Gate G-T3: met** - ECE reported (variant 0.017; diagnosis 0.141 on n=42); safety sensitivity 100%
at specificity 100%; monotone risk-coverage curve demonstrated.

Artifacts: `bench/track3_metrics.json`, `bench/track3_risk_coverage.csv`, `bench/track3_trustworthiness.py`.

---

## Limitations and the next steps

- **ClinVar-derived truth.** Track 1 uses ClinVar P/B labels (the secondary surface). Because GeneBe
  reproduces ClinVar, the pre-registered primary surface is **ClinGen eRepo (expert-panel)** and/or a
  **time-split** (grade on variants classified after each tool's data snapshot). The eRepo bleeding
  set already exists in the repo (Tier A1); running every tool on it is the next Track-1 step.
- **n=42 diagnosis set** is a proof-of-concept; the diagnosis ECE and risk-coverage are reported with
  that caveat. The cohorts carry the diagnosis headline (Paper 2).
- **Track 2 (vs LIRICAL / Exomiser).** Requires the Java tools + their data bundles on the VM and
  seeding each curated case's variant into a background exome; it is the diagnosis head-to-head and is
  the next heavy build. Not started this round.
- **Track 4 (vs DeepRare / agentic LLM).** Requires a frontier-LLM endpoint to run the hard cases N>=5
  times for variance + safety-flag + hallucination contrast. Lower engineering than Track 2; deferred.
- **Track 5 (reader study to TRIPOD+AI / DECIDE-AI standard).** Documentation + the OSF time-stamp
  (registration guide already written); deferred to the reader-study run.

## Reproducibility
```
python -m bench.genebe_client                 # fill the GeneBe cache (resumable; snapshot dated)
python -m bench.track1_variant_headtohead     # -> bench/track1_*.{csv,json}
python -m bench.track3_trustworthiness        # -> bench/track3_*.{json,csv}
python -m pytest tests/test_bench.py -q        # regression guards on the headline numbers
```
All comparator data is public; no patient data is used (Gate G7).
