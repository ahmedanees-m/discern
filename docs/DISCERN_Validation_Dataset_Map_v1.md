# Validation dataset map

**Date:** 2026-06-15 - **Companion to:** `DISCERN_Coverage_Architecture_v1.md`.
**Purpose:** the cited, accessibility-tagged catalog of every dataset that can prove DISCERN's two abilities, mapped to the diseases/genes it covers - and a clear statement of which abilities can be validated **now on public data** vs which remain **cohort-gated**.

> **Independently verified 2026-06-16:** every dataset is real; all citations/DOIs/PMIDs confirmed against primary sources. Metadata corrections folded in: GoldVariants = **2021**, DOI 10.1111/jth.15459, PMID 34355501 (+ corrigendum 36737374); EAHAD McVey 2020 = **PMID 32166871**, "free of charge / per-record web (no bulk API)" - not a copyright/research-use restriction; **CHAMP/CHBMP are publicly bulk-downloadable** (PMID 23280990 / 24498619); Bastida PMID 28983057; UK GAPP PMID 27479822; Panel-VWD = Ramanan *RPTH* 2025;9(2):102730 PMID 40242192; BRIDGE-BPD DOI 10.1182/blood.2018891192. DMS confirmed non-available for the covered genes.

---

## Requirements by capability

DISCERN does two coupled jobs; each needs a different kind of truth set:

1. **Variant / VUS classification** ("is this variant pathogenic, and can we reclassify a VUS?") -> needs **variant-level truth**: curated variants with expert classifications and/or functional/clinical labels. *Richly served by public data.*
2. **Differential diagnosis between confusable diseases** + the **coupling** ("which disease, and does phenotype sharpen the variant call?") -> needs **paired phenotype+genotype at the case level**. *Mostly cohort-gated.*

**Coverage being validated:** ~25 disease entities across **10 confusable clusters**; variant engine over **~80 genes** (ISTH TIER1-aligned).

---

## Part 1. Variant truth datasets, largely public

| Dataset | What it gives | Access | DISCERN coverage / use |
|---|---|---|---|
| **ClinGen eRepo** | VCEP expert classifications with per-code SEPIO evidence | **public** (.tab download) | Already used (A1): 12,240 variants / 170 genes; partition + combining-fidelity truth |
| **ClinVar** (`variant_summary`) | Community classifications + review-status stars (>=2*, 3* expert panel) | **public** (NCBI FTP) | Already used (A1/A3/H4): concordance + calibration + novel-variant labels |
| **gnomAD** v4.1 | Population allele frequencies | **public** (per-variant AFs to be pulled) | Already used (A1): frequency-rule (BA1/BS1/PM2) cross-check |
| **GoldVariants (ISTH SSC)** | Curated rare variants for **diagnostic-grade BTPD genes**, multi-center, bulk-transferred to ClinVar | **public / open-access** (REDCap resource) | **NEW - directly on-target.** Expands variant truth across *exactly* DISCERN's panel genes. **Folding status (2026-06-17): COVERED-BY-INCLUSION** - GoldVariants is bulk-deposited in ClinVar, so its variants are already inside DISCERN's ClinVar-derived validations (H3 n=10,780, H4 n=1,015, calibration n=7,521). A *dedicated named-subset* re-slice cannot be pulled from ClinVar's `[Submitter]` field (the 821 variants are attributed to the 30 individual submitting centers, not one queryable org - verified 8 name forms -> 0); it would need the ISTH GinTH twice-yearly download + an H4-style annotation pass, and the result would *mirror* the existing ClinVar concordance. Attribution-only; not gating Paper 1. |
| **EAHAD Coagulation Factor Variant DBs** (F7, **F8**, **F9**, **VWF**) | Per-variant genotype **+ phenotype (lab + clinical) + structure/function**; >3,000 F8, >1,200 F9 variants | **web-accessible, research-use** (f8-db.eahad.org, f9-db.eahad.org) | **NEW - partially paired.** Truth for the coagulation cluster (hemophilia A/B, FXIII) + VWF; the phenotype annotation supports *limited* coupling checks on these genes |
| **CHAMP / CHBMP (CDC)** | F8 / F9 mutation projects (severity, inhibitor risk) | **public** | Independent F8/F9 truth + clinical correlates |
| **LOVD** (Leiden Open Variation Database) | Gene-specific variant records (VWF and others; EAHAD mirrors here) | **public** | Cross-reference / additional per-gene truth |
| **HGMD** | Large curated disease-variant catalog | **licensed** (academic version limited) | Optional; cite-only unless licensed |

Between eRepo, ClinVar and gnomAD (used here) and GoldVariants, EAHAD and LOVD (public, not used here), variant classification can be validated across the covered genes well beyond the present run, without controlled access.

---

## Part 2. Paired phenotype and genotype datasets, largely access-controlled

| Dataset | What it gives | Access | Use |
|---|---|---|---|
| **Phenopacket Store (GA4GH)** | Standardized phenotype+gene case records | **public** | Already used (n=4 in-cluster); expandable as more bleeding cases are deposited |
| **Curated published cases** (case reports + cohort tables) | Real cases with phenotype + variant + diagnosis | **public (manual curation)** | Already used (n=10); **expandable to dozens** from the cohort/case literature below |
| **UK GAPP** (Genotyping & Phenotyping of Platelets) | WES + deep platelet phenotyping in inherited thrombocytopenia/PFD | **published cohort** (patient-level likely controlled) | Curate published case tables now; full data via collaboration/DAC |
| **Bastida 2018 / panel-IPD cohorts** | HTS diagnostic yield in inherited platelet disorders | **published** | Curated cases (Top-1/Top-3 benchmark) |
| **Panel-based VWD cohorts** (e.g. PMC12002656) | Pre/post-sequencing subtype changes incl. **2B vs PT-VWD, 2N vs HA** | **published** | Curated cases for exactly DISCERN's VWD clusters |
| **BRIDGE-BPD / NIHR BioResource** (`EGAS00001001172`) | **2,396 index patients** paired phenotype+genome (Downes 2019; not ~13k - corrected), bleeding/platelet | **controlled (EGA DAC)** | The primary coupling/diagnosis cohort - DAC application |
| **Genomics England 100,000 Genomes** | Paired genome+phenotype incl. bleeding/platelet panels | **controlled (Research Environment)** | Secondary coupling cohort |
| **UDN - Solve-RD** | Undiagnosed/rare-disease paired data | **controlled** | Off-critical-path secondary |
| **Disease registries** (WBDR clinical hemophilia/VWD registry; VWD/ITP/Glanzmann registries). *Note: EUHASS is treatment-safety pharmacovigilance, NOT a genomic/phenotype cohort - verified.* | Phenotype +/- genotype at scale | **controlled access** | A controlled-access paired cohort is the nearest-term paired source |

Case-level paired data for diagnostic accuracy and for the coupling is largely behind data-access or ethics approval. The curated published-case benchmark can be extended from public reports; the coupling endpoint still requires the controlled-access cohorts.

---

## Part 3. Functional and orthogonal truth

- **EAHAD structure/function annotations** - per-variant structural/functional consequence for F7/F8/F9/VWF (use as orthogonal support for PS3-type evidence).
- **MaveDB / ProteinGym (deep mutational scanning)** - *sparse for bleeding genes* (DMS is concentrated on cancer/known genes). **Check per covered gene**; do not assume availability. Where a DMS map exists for a covered gene, it is strong orthogonal functional truth for VUS.
- **ClinGen VCEP functional-evidence tables** - already feed PS3/BS3 routing.

---

## Part 4. Per-cluster data coverage

| Cluster | Key genes | Variant truth (Part 1) | Diagnosis truth (Part 2) |
|---|---|---|---|
| C1 integrin (GT/LAD-III/RASGRP2) | ITGA2B, ITGB3, FERMT3, RASGRP2 | **rich** (eRepo VCEP, ClinVar, GoldVariants) | curated cases; cohort gated |
| C2 VWF-GPIb (2B/PT-VWD/2A) | VWF, GP1BA | **rich** (VWF VCEP, EAHAD-VWF, GoldVariants, ClinVar) | panel-VWD published cases; cohort gated |
| C3 macrothrombocytopenia vs ITP | GP1BA/B/9, MYH9, ACTN1, TUBB1, FLNA | **good** (ClinVar, GoldVariants) | published cohorts; gated |
| C4 RUNX1/ANKRD26/ETV6 vs ITP | RUNX1, ANKRD26, ETV6 | **rich** (MM-VCEP, ClinVar, GoldVariants) | Joshi/Cooper ITP cases; gated |
| C5 2N VWD vs mild HA | VWF, F8 | **rich** (EAHAD-F8/VWF, VCEPs, GoldVariants) | EAHAD phenotype + curated cases |
| C6 storage-pool (HPS/Chediak) | HPS1-6, AP3B1, LYST | **moderate** (ClinVar) | case reports; gated |
| C7 alpha-granule (GPS/GFI1B/ARC/Quebec) | NBEAL2, GFI1B, VPS33B, PLAU | **thin-moderate** (ultra-rare) | case reports |
| C8 FXIII (F13A1 vs F13B) | F13A1, F13B | **good** (ClinVar, EAHAD/registries) | FXIII registry/case data |
| C9 type-1/low-VWF vs mild PFD | VWF + PFD genes | **rich** (VWF) | published cohorts |
| C10 Scott (ANO6) | ANO6 | **thin** (ultra-rare) | rare case reports |

**Reading it:** common, well-studied genes (F8, F9, VWF, ITGA2B/ITGB3, GP1BA, RUNX1, MYH9) have **abundant** variant truth across multiple public databases. The **ultra-rare** entities (ANO6/Scott, PLAU/Quebec, VPS33B/ARC, GFI1B) are data-thin **everywhere** - a real limit, mitigated only by case-report curation and the reduced-confidence tagging already in the engine.

---

## Part 5. Conclusion

- **Variant classification.** Public resources (GoldVariants, EAHAD F7/F8/F9/VWF, LOVD) together with eRepo, ClinVar and gnomAD cover the common genes in scope, so the variant engine can be validated more broadly than in the present run without controlled access.
- **Differential diagnosis & the coupling - PARTIALLY now, fully only gated.** The curated published-case benchmark can grow from ~10 to dozens from public cohort/case literature (improves the diagnosis-accuracy evidence). But the **coupling's headline validation** - calibrated phenotype sharpening the variant call on paired data - still requires the controlled cohorts (BRIDGE-BPD DAC, local Glanzmann IRB). No public dataset substitutes for it.
- **Ultra-rare clusters (C7/C10)** are data-thin in every database; validated weakly and tagged reduced-confidence.

The variant half can be validated broadly on public databases. The coupling half can be strengthened with curated published cases but is fully testable only on access-controlled cohorts.

---

## Sources
- GoldVariants: ISTH SSC Genomics subcommittee; *J Thromb Haemost* (open-access REDCap resource; 821 initial variants, 30 centers/14 countries; bulk -> ClinVar).
- EAHAD Coagulation Factor Variant DBs (F7/F8/F9/VWF): McVey et al., *Haemophilia* 2020;26(2):306-313 (DOI 10.1111/hae.13947); f8-db.eahad.org, f9-db.eahad.org. CHAMP/CHBMP: CDC Hemophilia A/B Mutation Projects.
- ISTH TIER1 gene list: Megy et al., *J Thromb Haemost* 2019;17(8):1253-1260 (PMID 31179617).
- UK GAPP: Johnson et al., *Haematologica* 2016;101:1170-9. Bastida et al., *Haematologica* 2018;103:148-162.
- Panel VWD: *Clinical utility of panel-based genetic sequencing for VWD* (PMC12002656).
- BRIDGE-BPD: `EGAS00001001172` (NIHR BioResource); Downes et al., *Blood* 2019 (PMID 31064749).
- ClinVar / gnomAD / ClinGen eRepo / Phenopacket Store / LOVD / HGMD - standard public (HGMD licensed).

*Verification note: dataset contents and access terms drift - re-check GoldVariants/EAHAD/LOVD bulk-access and licensing at use time; EAHAD content is research-use and copyright-protected (per-record access, not necessarily bulk download).*

*End of DISCERN Validation Dataset Map v1*
