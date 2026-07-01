# DISCERN - CHAMP/CHBMP Independent Sensitivity Benchmark (v1)

**Date:** 2026-06-16 - **Code:** `eval/champ_chbmp_benchmark.py` - **Track:** A2 extension (variant engine).
**Datasets (verified, public, PHI-free):** CDC **CHAMP** (Hemophilia A Mutation Project, F8) and
**CHBMP** (Hemophilia B Mutation Project, F9), 2022 lists, downloaded 2024-11-07 from
`cdc.gov/hemophilia/mutation-project`. Citations: Payne 2013 (PMID 23280990), Li 2013 (PMID 24498619).
Verified against primary sources.

## What this tests
CHAMP/CHBMP catalog variants **reported in patients with haemophilia** - an independent, curated set
of disease-causing F8/F9 alleles (with ISTH-criteria severity + inhibitor history), assembled by a
source entirely separate from ClinGen eRepo / ClinVar. Truth label for every catalogued allele =
**pathogenic**, so this is a **sensitivity (recall)** benchmark of DISCERN's gene-specific
(CFD-VCEP F8/F9) scoring, broken out by consequence.

**Design choice (conservative):** the NULL subset is scored on **consequence alone** - the PVS1
decision tree (Abou Tayoun 2018), with **no gnomAD frequency and no predictor input**. PM2 is
deliberately withheld; because every failing call below is capped at PVS1_Strong/Moderate, adding
PM2_Supporting (these alleles are absent from gnomAD) never flips a call - so 91.2% is a floor.

## Catalog composition
| Variant type | F8 (n=4,038) | F9 (n=1,399) |
|---|---|---|
| Missense | 1,803 | 755 |
| Frameshift | 1,052 | 244 |
| Nonsense | 460 | 124 |
| Splice site change | 351 (185 canonical +/-1,2) | 122 (65 canonical +/-1,2) |
| Large structural (>50 bp) | 236 | 50 |
| Small structural (in-frame) | 86 | 36 |
| Synonymous / UTR / promoter | 50 | 68 |

## Result 1 - NULL-subset recall (LP/P), consequence only

| Gene | Consequence | n | LP/P | Recall |
|---|---|---|---|---|
| F8 | canonical splice (+/-1,2) | 185 | 185 | **100.0%** |
| F8 | frameshift | 1,052 | 1,024 | **97.3%** |
| F8 | nonsense | 460 | 449 | **97.6%** |
| F9 | canonical splice (+/-1,2) | 65 | 65 | **100.0%** |
| F9 | frameshift | 244 | 148 | 60.7% |
| F9 | nonsense | 124 | 72 | 58.1% |
| **F8 all-null** | | **1,697** | **1,658** | **97.7%** |
| **F9 all-null** | | **433** | **285** | **65.8%** |
| **Overall** | | **2,130** | **1,943** | **91.2%** |

**The F9 gap is a real, explainable finding, not a miss.** F9's terminal **exon 8 is unusually large**
(it encodes the entire serine-protease catalytic domain). Per ClinGen SVI rules, PTCs in the last
exon escape NMD and receive **PVS1_Strong/Moderate**, not Very Strong - 4 or 2 Tavtigian points =
**VUS on consequence alone**. The engine is *correctly cautious*: a C-terminal truncation's impact
depends on how much of the catalytic domain is lost, which needs functional/segregation evidence
(which the catalog has, via patient phenotype - i.e. the coupling layer). All 187 non-pathogenic
calls are last-exon NMD-escape truncations; none are mis-scored against the rules.

Non-canonical "splice site change" entries (166 F8 / 57 F9 - e.g. +3, exonic ESE, deep-intronic)
are **excluded** from the PVS1 subset by design (they need a SpliceAI/Pangolin call, not PVS1).

## Result 2 - MISSENSE intrinsic-only ceiling (why coupling is needed)
With the CFD-VCEP F8/F9 spec, a missense variant tops out at **PP3_Supporting + PM2_Supporting =
2 points = VUS**, even with a maximal predictor (REVEL 0.95) and absence from gnomAD. So all
2,558 catalogued missense (1,803 F8 + 755 F9) are **VUS on intrinsic+predictor evidence alone** -
they require routed PS3 (functional) / PP1 (segregation) / **PP4 (the disease<->variant coupling)**
to reach LP. This is the designed limitation that motivates Paper 2, demonstrated on real disease
alleles: the variant engine confidently resolves LOF/null variants, but missense resolution is
coupling-gated. (The H4 benchmark already showed that *with* REVEL PP3 + ClinVar PS1/PM5,
DISCERN's gene-specific scoring beats literal InterVar - AUC 0.912 vs 0.874 - on the ClinVar
missense set.)

## Result 3 - MISSENSE arm: routed PS1/PM5 + PP3 (empirical; `champ_chbmp_missense_arm.py`)
The natural extension: score each catalogued missense as DISCERN would a **novel** variant, adding
the routed ClinVar codes (**PS1** = same AA change pathogenic in ClinVar; **PM5** = different change,
same residue) from `adapters/clinvar.py` to PP3 (REVEL) + PM2. Anti-circularity: only variants whose
**exact cDNA is NOT in ClinVar** (truly novel) are credited - so PS1/PM5 always come from a
*different* variant. (Local ClinVar `variant_summary_20260503`; F8 has 281 residues / 375 pathogenic
AA changes in the index, F9 136 / 202.)

| Gene | catalogued missense | novel (cDNA not-in ClinVar) | PS1 | PM5 | neither |
|---|---|---|---|---|---|
| F8 | 1,803 | 1,307 | 16 (**1.2%**) | 346 (26.5%) | 945 (72.3%) |
| F9 | 755 | 511 | 13 (**2.5%**) | 242 (47.4%) | 256 (50.1%) |

**Point arithmetic bounds the recall exactly.** A novel missense reaches LP (>=6 pts) only as
**PS1(4)+PM2(1)+PP3(1)=6**; PM5(2)+PM2(1)+PP3(1)=4 and no-hit=2 both stay **VUS**. So the **PS1 rate
is the hard ceiling on novel-missense LP recovery: 1.2% (F8) / 2.5% (F9)** by the variant+ClinVar
engine alone - PM5 (26-47%) and no-hit (50-72%) cannot reach LP regardless of predictor.

**Finalized with the VM REVEL pass (2026-06-16).** REVEL for the 29 PS1 cases was obtained on the VM
(ANNOVAR `hg19_revel`, `table_annovar` refGene+revel; ANNOVAR's recomputed AAChange validated the
strand/coordinate of every variant), and gnomAD absence was confirmed via the gnomAD v2.1.1 API
(**28/28 SNVs gnomAD-absent -> PM2 fires for all**). Result:

| Gene | PS1 cases | REVEL>=0.6 (PP3) | gnomAD-absent (PM2) | **LP/P recall** |
|---|---|---|---|---|
| F8 | 16 | 16/16 (0.647-0.992) | 16/16 | **16/1307 = 1.2%** |
| F9 | 13 | 11/12 SNV (p.Ala279Thr=0.456 fails; 1 delins not SNV-scored) | 12/12 | **11/511 = 2.2%** |
| **Both** | 29 | 27 | 28/28 | **27/1818 = 1.5%** |

So the exact figure (1.5% combined) lands right at the PS1 ceiling, REVEL trimming F9 from 2.5% to
2.2% (one sub-cutoff substitution + one delins REVEL cannot score). Two CHAMP records (F8 c.74A>C,
F9 c.836C>T) carry cDNA/protein annotation inconsistencies - ANNOVAR's RefSeq consequence (Y25S,
A279V) differs from the catalog's label (Y25C, A279T); flagged, <=2-case effect.

**This is not a DISCERN weakness - it is intrinsic to ACMG missense interpretation** (InterVar and
any code-based classifier face the identical ceiling; REVEL/AlphaMissense *rank* but do not
*classify*). It quantifies, on real disease alleles, exactly how much is left to the functional
(PS3) / segregation (PP1) / **disease-coupling (PP4)** layers: ~97.5-98.8% of novel missense disease
alleles are unreachable by sequence+frequency+ClinVar alone. **Lever:** if the CFD-VCEP permits
**PM5_Strong** at residues with multiple distinct pathogenic substitutions, the PM5 fraction
(26-47%) becomes partly LP-reachable - DISCERN currently applies PM5 at the conservative Moderate.

## Framing
- **Strengthens:** an *independent, separately-curated* truth set (CDC, not ClinGen/ClinVar)
  confirms DISCERN classifies LOF/null F8/F9 disease alleles as LP/P at **91-98%** on consequence
  alone, with the one systematic gap (F9 terminal exon) fully explained by NMD-escape rules.
- **Does not claim:** missense sensitivity from intrinsic evidence (ceiling = VUS by design), nor
  any coupling result (still gated). No specificity arm here (the catalog has no benign controls);
  specificity is covered by the eRepo/ClinVar calibration benchmarks (Brier 0.0073, ECE 0.008).
- **No PHI.** Public CDC variant catalog; the raw .xlsx are not committed (download URL + date recorded).
