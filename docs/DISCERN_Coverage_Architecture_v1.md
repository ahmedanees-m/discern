# DISCERN - Coverage Architecture v1 (literature-grounded)

**Version:** 1.0.0 - **Date:** 2026-06-13
**Purpose:** the authoritative, cited answer to *which diseases, which confusable clusters, which genes* DISCERN must cover - derived from the field's classification bodies, curated gene panels, and diagnostic-pitfall literature (not from clinical intuition). This document supplies the concrete content for Track A (variant gene panel) and Track B (diagnosis clusters).

---

## 0. How coverage was determined (sources of truth)

Coverage is anchored to the bodies that define this space, so every disease/gene/cluster is traceable:

- **Gene set:** ISTH SSC **curated gene list** for bleeding, thrombotic, and platelet disorders - Megy, Downes, Simeoni et al., *J Thromb Haemost* 2019;17(8):1253-1260 (TIER1/TIER2; list maintained at isth.org, updated July 2024). Cross-checked against diagnostic panels (e.g. Mayo GNPLT, 70 platelet-disorder genes).
- **Platelet-disorder taxonomy & diagnosis:** Gresele/SSC *Diagnosis of inherited platelet function disorders* (*J Thromb Haemost* 2015); Lentaigne, Freson, Laffan et al. *Inherited platelet disorders: toward DNA-based diagnosis* (*Blood* 2016;127(23):2814-2823); recent diagnostic review (*Biomolecules* 2025;15(6):846).
- **VWD classification:** *Diagnosis of von Willebrand disease* (*Blood Adv* 2025;9(22):5870); *How I treat type 2 VWD* (*Blood* 2015;125(6):907); NHLBI/ASH-ISTH-NHF-WFH guidelines.
- **Rare bleeding disorders (coagulation):** ASH Education *Treatment of rare factor deficiencies* (2016); de Moerloose et al. (*Haemophilia* 2016, PMID 27405678); Franchini (*Ann Transl Med* 2018).
- **Inherited-thrombocytopenia <-> ITP misdiagnosis & malignancy predisposition:** *Genetic sequencing in chronic ITP* (*Blood Adv* 2025;9(7):1497); Galera et al. (*Int J Lab Hematol* 2019, doi:10.1111/ijlh.12999); Zaninetti et al. (*J Thromb Haemost* 2017, doi:10.1111/jth.13855); WHO germline-predisposition category.

**Design principle (from v3).** DISCERN's value has two faces, and "coverage" means something different for each:
1. **Variant engine ->** breadth: cover the **entire curated gene panel** (the partition + ACMG scoring is gene-agnostic, so wide coverage is nearly free and is the generalization claim).
2. **Diagnosis engine ->** depth on the **confusable clusters**: the finite, well-documented set of look-alikes where misdiagnosis causes real harm. This set is bounded - covering *all* of it is what makes the tool great, not padding with non-confusable disorders.

---

## PART 1 - The disease/gene taxonomy (what exists)

Four domains. Genes listed are the load-bearing ones; the variant engine targets the full TIER1 panel (Part 4).

### 1. Inherited platelet disorders (IPDs)

**1A. Thrombocytopenias (low count).**
- *Macrothrombocytopenia (large/giant platelets):* **MYH9**-related disorder (May-Hegglin/Sebastian/Fechtner/Epstein); **Bernard-Soulier syndrome** (GPIb-IX-V: **GP1BA/GP1BB/GP9**, biallelic; monoallelic *GP1BA/GP1BB* also); **ACTN1**; **TUBB1**; **FLNA** (X-linked); **DIAPH1**; **TPM4**.
- *Normal/small-platelet thrombocytopenia:* **WAS/X-linked thrombocytopenia** (micro-platelets); **ANKRD26**-RT; **RUNX1**/FPDMM; **ETV6**-RT; **MPL** (CAMT); **THPO**; **CYCS**; **GP1BA** (some); **TAR/RBM8A**.

**1B. Platelet function disorders (normal/near-normal count).**
- *Adhesion (GPIb-VWF axis):* Bernard-Soulier (GP1BA/GP1BB/GP9); **platelet-type VWD** (**GP1BA** gain-of-function).
- *Aggregation (aIIbb3 / inside-out signaling):* **Glanzmann thrombasthenia** (**ITGA2B/ITGB3**); **LAD-III** (**FERMT3**/kindlin-3); **RASGRP2**/CalDAG-GEFI.
- *delta-granule / storage-pool (often syndromic):* Hermansky-Pudlak (**HPS1, HPS3-6, AP3B1, DTNBP1, BLOC1S3/6**, +); Chediak-Higashi (**LYST**); Griscelli; nonsyndromic delta-SPD.
- *alpha-granule:* Gray platelet syndrome (**NBEAL2**); **GFI1B**; ARC syndrome (**VPS33B/VIPAS39**); Quebec platelet disorder (**PLAU** duplication).
- *Signaling/secretion/receptors:* **P2RY12**; thromboxane receptor **TBXA2R**; **TBXAS1**; **GNAS**; **PLA2G4A**; **FLI1**; **SRC**.
- *Procoagulant (PS exposure):* **Scott syndrome** (**ANO6**).

**1C. Germline predisposition to myeloid/lymphoid malignancy (overlaps 1A).** **RUNX1** (FPDMM -> AML/MDS), **ANKRD26** (THC2 -> MDS/AML), **ETV6** (THC5 -> ALL/MDS). Among inherited thrombocytopenias these have prevalence ~3% / 18% / 5% respectively (Galera 2019). WHO recognizes "myeloid neoplasms with germline predisposition and pre-existing platelet disorder" for exactly these three.

### 2. Von Willebrand disease (VWF; GP1BA for the PT-VWD mimic)
Quantitative: **Type 1** (~75%), **1C** (accelerated clearance), **Type 3** (near-absent). Qualitative: **2A** (loss of HMW multimers), **2B** (GPIb gain-of-function; ~40% thrombocytopenic), **2M** (reduced GPIb/collagen binding, normal multimers), **2N** (reduced FVIII binding - mimics hemophilia A). Plus the **"low VWF"** entity (~35% lack a VWF variant) and **acquired VWS** (a non-inherited mimic).

### 3. Coagulation factor deficiencies
- *Common:* **Hemophilia A** (**F8**), **Hemophilia B** (**F9**), **FXI deficiency** (**F11**).
- *Rare bleeding disorders (RBDs, ~3-5% of factor deficiencies, mostly AR):* fibrinogen (**FGA/FGB/FGG**; a-/hypo-/dysfibrinogenemia), prothrombin (**F2**), **FV** (**F5**), **combined FV+FVIII** (**LMAN1/MCFD2**), **FVII** (**F7**), **FX** (**F10**), **FXIII** (**F13A1/F13B**), **FXII** (**F12**), **combined vitamin-K-dependent factors** (**GGCX/VKORC1**).

### 4. Fibrinolytic / vascular / connective-tissue bleeding (mimics & overlap)
alpha2-antiplasmin deficiency (**SERPINF2**), PAI-1 deficiency (**SERPINE1**), hereditary hemorrhagic telangiectasia (**ENG/ACVRL1/SMAD4**), vascular Ehlers-Danlos (**COL3A1**). *(Lower priority for v1; relevant as differentials for "bleeding with normal coagulation tests.")*

---

## PART 2 - The confusable-cluster catalog (what DISCERN diagnoses)

These are the documented diagnostic **traps** - where two+ disorders look alike but diverge in management or surveillance. Each entry is a DISCERN cluster. *** = patient-safety-critical** (wrong call -> contraindicated drug, unnecessary surgery, or missed cancer).

| # | Cluster | Constituent disorders (gene) | The decisive observation | Management divergence / stakes | Source | Status |
|---|---|---|---|---|---|---|
| **C1** | Aggregation defect, normal count | Glanzmann (ITGA2B/ITGB3) - **LAD-III** (FERMT3) - RASGRP2 | Leukocytosis + recurrent infection; flow aIIbb3; PAC-1 | * LAD-III needs **HSCT** + has immunodeficiency; GT/RASGRP2 do not | Gresele 2015; OMIM | [verified] built |
| **C2** | Enhanced/abnormal RIPA | VWD **2B** (VWF) - **PT-VWD** (GP1BA) - VWD 2A | **RIPA mixing** (plasma vs platelet origin); GP1BA vs VWF sequencing | * **DDAVP contraindicated in 2B**; 2B->VWF concentrate vs PT-VWD->platelets | Blood 2015;125:907; Blood Adv 2025 | [verified] built |
| **C3** | Macrothrombocytopenia vs ITP | BSS (GP1BA/B/9) - MYH9-RD - ACTN1 - TUBB1 - FLNA - **vs ITP** | Lifelong/familial low count, **giant platelets**, MYH9 inclusions, normal MPV-for-gene | * Avoid unnecessary **splenectomy / steroids / IVIG**; giant platelets undercounted by analyzers | Lentaigne 2016; Blood Adv 2025;9:1497 | [pending] build |
| **C4** | Normal-size IT + malignancy risk vs ITP/MDS | **RUNX1** - **ANKRD26** - **ETV6** - vs ITP, vs sporadic MDS | Family history, normal platelet size, AD pattern, marrow dysmegakaryopoiesis; germline panel | ** **Missed leukemia surveillance**; avoid immunosuppression; **affected relative must not be HSCT donor** | Galera 2019; Zaninetti 2017; Blood Adv 2025;9:1497 | [pending] build (highest priority) |
| **C5** | Isolated low FVIII | VWD **2N** (VWF) - **mild/moderate hemophilia A** (F8) | **VWF:FVIII-binding assay** + genetics; FVIII:C/VWF:Ag ratio | Inheritance (AR vs X-linked) -> counseling; treatment (VWF-containing vs FVIII); female carriers | Blood Adv 2025;9:5870; Mayo VWD8B | [pending] build |
| **C6** | delta-granule / storage-pool with syndromic clues | HPS (HPS1/3-6, AP3B1) - Chediak (LYST) - nonsyndromic delta-SPD | Oculocutaneous albinism, pulmonary fibrosis/colitis (HPS); immunodeficiency/HLH (Chediak); EM dense-granule count | * Chediak->**HSCT** + HLH risk; HPS->pulmonary/colitis surveillance | Gresele 2015; Lentaigne 2016 | [pending] build |
| **C7** | alpha-granule deficiency | Gray platelet (NBEAL2) - GFI1B - ARC (VPS33B/VIPAS39) - Quebec (PLAU) | Gray platelets on smear, EM alpha-granule absence; myelofibrosis (GPS); urokinase (Quebec) | Myelofibrosis surveillance (GPS); Quebec needs **antifibrinolytics, not platelets** | Lentaigne 2016 | [pending] build |
| **C8** | FXIII deficiency (subunit matters) | **F13A1** vs **F13B** | FXIII activity (normal routine coags), antigen, subunit/genetics; delayed umbilical/ICH, recurrent miscarriage | * **Recombinant FXIII-A2 works for F13A1 but not F13B**; lifelong prophylaxis; high ICH risk | de Moerloose 2016; FXIII review (PMC12825037) | [pending] build |
| **C9** | Mild bleeding, normal/borderline tests | **Type 1 VWD / low VWF** (VWF) - mild platelet function defect - normal | Repeat VWF panel, LTA, bleeding-assessment-tool score | Over-/under-diagnosis; "low VWF" is not classic VWD; avoid over-treatment | Blood Adv 2025;9:5870; Gresele 2015 | [pending] build |
| **C10** | Procoagulant defect missed by routine tests | **Scott syndrome** (ANO6) | PS-exposure / prothrombinase assay (normal LTA & aggregation) | Easily missed (standard platelet tests normal); transfusion planning | Lentaigne 2016 | [pending] optional |

**Coverage logic.** C1-C2 are built. C3-C9 are the documented, harm-bearing traps that complete the diagnosis engine for this domain; **C4 is first** (the leukemia-predisposition group has the gravest miss). C10 is optional (rare, but a clean demonstration of "tests normal yet bleeding"). Beyond these, additional disorders are *catalog entries for the variant engine* but not separate diagnosis clusters unless a genuine confusion is documented - that is the line that keeps the tool great rather than bloated.

---

## PART 3 - Build prioritization (by harm, not by effort)

Ordered by **clinical stakes x misdiagnosis frequency x whether gene-specific ACMG rules exist** (VCEP coverage strengthens the coupled variant call). Per the directive, this is *not* gated by effort - it is the order that maximizes the tool's clinical value.

1. **C4 - RUNX1/ANKRD26/ETV6 vs ITP/MDS.** Gravest miss (cancer surveillance + donor selection); high real-world misdiagnosis rate; **RUNX1 VCEP** exists. *Build first.*
2. **C3 - Macrothrombocytopenia vs ITP.** Very common misdiagnosis -> unnecessary splenectomy; **Platelet Disorders VCEP** covers BSS/GP1BA. 
3. **C5 - 2N vs mild hemophilia A.** Clean decisive test (VWF:FVIIIB); **VWD VCEP (VWF)** + **Coagulation VCEP (F8)** both exist - strong coupled-variant story across two genes.
4. **C8 - FXIII (F13A1 vs F13B).** Genotype changes the *drug*; dramatic, well-bounded.
5. **C6 - Storage-pool/HPS/Chediak.** Syndromic clues route to HSCT/surveillance; high-value but more genes.
6. **C7 - alpha-granule (gray platelet etc.).** Myelofibrosis surveillance; moderate frequency.
7. **C9 - mild VWD / low VWF.** Highest *volume*, lowest per-case drama; valuable for calibration/abstention demonstration.
8. **C10 - Scott syndrome.** Optional capstone.

Each new cluster is **YAML + sourced LRs only** (no engine change) - the generalization claim, and the reason the order above costs build-time, not redesign.

---

## PART 4 - Variant-engine gene panel (breadth target)

The variant side targets the **full ISTH TIER1 curated panel** (Megy 2019, maintained list) across all four domains, not just cluster genes - this is what the genome-wide partition validation (v3 Track A1) and novel-variant scoring (A2) run on. Concretely the panel spans, at minimum:

- **Platelet (~70 genes):** ITGA2B, ITGB3, FERMT3, RASGRP2, ITGB2, GP1BA, GP1BB, GP9, MYH9, ACTN1, TUBB1, FLNA, DIAPH1, TPM4, WAS, ANKRD26, RUNX1, ETV6, MPL, THPO, CYCS, RBM8A, NBEAL2, GFI1B, VPS33B, VIPAS39, PLAU, HPS1, HPS3, HPS4, HPS5, HPS6, AP3B1, DTNBP1, BLOC1S3, BLOC1S6, LYST, P2RY12, TBXA2R, TBXAS1, GNAS, PLA2G4A, FLI1, SRC, ANO6, ... (full TIER1 list).
- **VWD:** VWF (+ GP1BA for PT-VWD).
- **Coagulation:** F8, F9, F11, F2, F5, F7, F10, F12, F13A1, F13B, FGA, FGB, FGG, LMAN1, MCFD2, GGCX, VKORC1.
- **Vascular/fibrinolytic (catalog):** SERPINF2, SERPINE1, ENG, ACVRL1, SMAD4, COL3A1.

**Gene-specific ACMG rules (ClinGen VCEPs) that anchor the coupled variant call** - already partly extracted (v2): Platelet Disorders VCEP (**ITGA2B/ITGB3**, **GP1BA**/BSS), Coagulation Factor Deficiency VCEP (**F8, F9**), von Willebrand Disease VCEP (**VWF**), RUNX1 VCEP (**RUNX1**). Genes outside a VCEP fall back to general ACMG with reduced confidence (Gate G2) - and the genome-wide validation (A1) reports exactly how many panel variants that affects.

---

## PART 5 - How this updates the v3 plan

- **Track A (variant):** panel = Part 4 (full TIER1), not the 38-gene bleeding subset; A1's genome-wide run and A2's novel-variant scoring use it. VCEP anchors per Part 4.
- **Track B1 (clusters):** replace the sketched candidate list with the **C1-C10 catalog (Part 2)**, built in the **Part 3 priority order (C4 first)**. Each cluster's LRs and pertinent negatives are curated to the same `(freq, n_cases, PMID)` bar as the VCEP specs; each cluster's decisive observation and management divergence are encoded for the next-observation (EIG) and safety-interlock modules.
- **Safety interlock (B3):** the divergence/contraindication map extends to every new cluster - at minimum: **DDAVP!->2B** (C2, fix the leading-call bug), **avoid splenectomy/immunosuppression** when an inherited thrombocytopenia is plausible (C3/C4), **recombinant FXIII-A2!->F13B** (C8), **antifibrinolytics-not-platelets** for Quebec (C7), and **HSCT-relevant** flags for LAD-III/Chediak (C1/C6). Every flag sourced.
- **Benchmark (B4):** the curated published-case set is populated **per cluster** from the case literature cited above (e.g. ANKRD26-as-MDS, ANKRD26/RUNX1-as-ITP), giving real, cited cases across C3-C8.

---

## Sources (load-bearing)

1. Megy K, Downes K, Simeoni I, et al. Curated disease-causing genes for bleeding, thrombotic, and platelet disorders: communication from the SSC of the ISTH. *J Thromb Haemost.* 2019;17(8):1253-1260. (TIER1 list, isth.org, updated Jul 2024.)
2. Gresele P; SSC. Diagnosis of inherited platelet function disorders: guidance from the SSC of the ISTH. *J Thromb Haemost.* 2015.
3. Lentaigne C, Freson K, Laffan MA, et al. Inherited platelet disorders: toward DNA-based diagnosis. *Blood.* 2016;127(23):2814-2823.
4. Diagnosis of von Willebrand disease. *Blood Adv.* 2025;9(22):5870.
5. How I treat type 2 variant forms of von Willebrand disease. *Blood.* 2015;125(6):907.
6. Treatment of rare factor deficiencies in 2016. *ASH Hematology Educ Program.* 2016.
7. de Moerloose P, Schved JF, Nugent D. Rare coagulation disorders: fibrinogen, FVII and FXIII. *Haemophilia.* 2016. PMID 27405678.
8. Galera P, et al. Inherited thrombocytopenia and platelet disorders with germline predisposition to myeloid neoplasia. *Int J Lab Hematol.* 2019. doi:10.1111/ijlh.12999.
9. Zaninetti C, et al. Inherited thrombocytopenia caused by ANKRD26 mutations misdiagnosed and treated as MDS. *J Thromb Haemost.* 2017. doi:10.1111/jth.13855.
10. The role of genetic sequencing in the diagnostic workup for chronic ITP. *Blood Adv.* 2025;9(7):1497.
11. Factor XIII deficiency: biology, testing, treatment. *(PMC12825037.)*

*Verification note: the ISTH TIER1 list is the canonical, version-controlled gene set and should be pulled fresh at build time (it updates). Gene-disease and cluster claims above are traceable to the sources cited; new clusters undergo the same source-verification pass as the VCEP specs.*

*End of DISCERN Coverage Architecture v1.0.0*
