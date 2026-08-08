# Supplemental figure and table legends

**Figure S1. Sensitivity of the diagnosis arm to the strength of the gene term.**
Curated Top-1 accuracy, hard-stop sensitivity and hard-stop specificity across the full range of the
P(G|D) likelihood ratio, from an inert gene term (ratio 1) to a ratio of 19. Safety behavior is
invariant throughout. The diagnosis result is flat above the committed value of 4, which is the
largest ratio at which a cluster's deciding observation can still overturn the gene; shaded regions
mark ratios at which it cannot. Related to Figure 4 and the Material and methods.

**Figure S2. Treatment-divergence safety interlock, by scenario.**
Each of five treatment-divergent scenarios against the two behaviors the interlock must exhibit:
firing when the planned management is contraindicated for a disease that has not been excluded, and
remaining silent when the planned management is harmless. Sensitivity and specificity are both
100%. Adjudication is performed on a gene-blind posterior against the engine's own leading call, so
that a variant in one gene cannot silently retire a treatment-divergent disease of another. Related
to the trustworthiness results.

**Figure S3. Diagnosis-posterior calibration on the curated benchmark.**
Mean predicted confidence against observed accuracy on the curated case set (n = 42), with the
expected calibration error. The shaded band marks the region above the 0.8 confidence threshold, in
which no incorrect call was made. Reported as a proof-of-concept characterization: with the
diagnosis arm at ceiling on this benchmark the calibration figure reflects mild underconfidence
rather than a discriminating result, and the abstention layer's operational value requires a
benchmark that is not at ceiling. Related to the withdrawn risk-coverage claim.

**Figure S4. Independent sensitivity on the CDC hemophilia mutation projects.**
Likely-pathogenic or pathogenic recall on null variants in the CHAMP and CHBMP catalogs (2,130 null
variants of 5,437 disease alleles), scored by consequence alone with no frequency or predictor
input. The lower F9 figure reflects that gene's large terminal exon, where premature termination
codons escape nonsense-mediated decay and are correctly held at uncertain significance rather than
called pathogenic. Related to the independent sensitivity result.

**Figure S5. A worked example: the safety interlock on a single case.**
(A) The evidence supplied: a GP1BA missense variant absent from gnomAD, enhanced low-dose
ristocetin-induced platelet aggregation, a mixing study of platelet origin, plasma origin
explicitly excluded, thrombocytopenia, and desmopressin as the planned therapy.
(B) The criterion trail. PM2 is applied and owned by the variant factor. PP4 and PS3 are owned by
the disease and functional factors and are not available from sequence; PM5 and PS4 are
variant-intrinsic but have no input under this protocol. The variant therefore remains of uncertain
significance.
(C) The ranked differential, in which platelet-type von Willebrand disease leads at 99.3%.
(D) The decision output. Desmopressin is contraindicated in type 2B von Willebrand disease, a VWF
disease, whereas the variant here is in GP1BA. Because safety is adjudicated gene-blind, type 2B
retains 0.9% in the safety view and the hard stop fires, together with the sequencing test that
would resolve the question. This is the framework's clinical argument on one case: a low-probability
but treatment-divergent competitor is not silently retired by the gene. Related to Figure 1.

---

**Table S1. Per-criterion agreement with expert-panel applications.**
Every ACMG/AMP criterion encountered on the Evidence Repository surface, with the number of times
the expert panels applied it, the number of times DISCERN applied it, Cohen's kappa where both
applied it at least once, the number of variants on which both applied it, and the fraction of
those on which both applied it at the same ClinGen strength. Computed on the 1,265 variants
carrying a pathogenic or benign expert classification, since agreement requires a graded variant on
both sides; each kappa reconciles with the applied counts in its own row at that denominator. An
unmodified code is resolved to its criterion's default strength before the strength comparison.
Criteria applied zero times carry the reason, in the same three categories Table S2 uses:
a protocol choice, a data-availability limit, or an engine scope gap.

**Table S2. Attribution of the intrinsic-evidence ceiling.**
Variants reaching Likely Pathogenic under each restoration condition: as scored; restoring the
variant-intrinsic criteria this pipeline cannot annotate; restoring the criteria the partition routes
to other factors; and restoring both. Followed by the reason each unapplied criterion has no input,
classified as a protocol choice, a data-availability limit, or an engine scope gap. Counts use
default criterion strengths and are therefore lower bounds.

**Table S3. Discrimination likelihood ratios for the ten confusable-disease clusters.**
Every likelihood ratio in the disease-discrimination model, given as the probability of the finding
under each disease, with the sample size and the PubMed identifier of the primary source. Continuous
integration fails if any entry lacks a source. Supplied as a separate Excel file.

**Table S4. Treatment-divergence scenarios.**
Each safety scenario with the disease pair, the planned management, the expected behavior, and the
observed behavior, for both the contraindicated and the harmless plan.

**Table S5. The curated published-case benchmark.**
All 42 cases with PubMed identifier, causal gene, cluster, expected diagnosis, extracted Human
Phenotype Ontology terms, stratum, and the leading call and Top-1 outcome for DISCERN and for the
phenotype-blind gene lookup. No article text is reproduced. Supplied as a separate Excel file.

**Table S6. CHAMP and CHBMP recall by gene.**
Likely-pathogenic or pathogenic recall on null variants for F8, F9, and the combined catalogs, with
the explanation of the F9 result.

**Table S7. Third-party data sources.**
Each source with the files used, the access URL, the version or snapshot date, and its role. None is
redistributed; checksums for the exact files used appear in the deposit manifest. Supplied as a
separate Excel file.

**Table S8. LIRICAL comparison arms.**
Every arm of the external comparison with its sample size, recall at 1, 3 and 5, and mean reciprocal
rank, followed by the paired tests against LIRICAL restricted to the same cluster. The
hpo_representable_only arm supports claims about inference on identical evidence; the phenotype_only
arm supports the architectural claim about what the model can represent. The post-correction arm is
in-sample and is not a head-to-head result.

**Table S9. Software, tool, and database versions.**
Every component used to produce the reported results, with its version and role. DISCERN's own
commit hash is the authoritative version identifier. Supplied as a separate Excel file.
