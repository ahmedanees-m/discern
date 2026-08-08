# Reproducibility checklist

Every reported number is reproducible from a clean clone and the named open datasets. No
third-party database is redistributed; `data/manifest.json` records each source with its URL,
version and checksum so the pinned inputs can be recovered.

## Environment

```bash
git clone https://github.com/ahmedanees-m/discern && cd discern
conda env create -f environment.yml        # or: pip install -e ".[dev]"
make test                                  # ruff and pytest
```

Python 3.11 or later. The test suite runs without the heavy bioinformatics dependencies; the
harnesses that read large third-party files require them.

## Datasets

Pull each source fresh at run time and record the version, since upstream releases drift.

| Dataset | Source | Version used here |
|---|---|---|
| ClinGen Evidence Repository | erepo.clinicalgenome.org, `.tab` export | 2026-05 export, 170 genes, 12240 variants |
| ClinVar | NCBI FTP, `tab_delimited/variant_summary.txt.gz` | 2026-05-04 |
| gnomAD | gnomad.broadinstitute.org | v4.1 callset, v4.1.1 annotations |
| REVEL | dbNSFP release | pinned at run time |
| AlphaMissense | Zenodo record 8360242, CC BY 4.0 | pinned at run time |
| Human Phenotype Ontology | hpo.jax.org | release recorded in the manifest |
| CDC CHAMP and CHBMP | cdc.gov hemophilia mutation projects | 2022 variant lists |
| GA4GH Phenopacket Store | monarch-initiative/phenopacket-store | release recorded in the manifest |

## Reproducing the reported numbers

Each command writes the JSON file that the tests and `RESULTS.md` read.

| Result | Command |
|---|---|
| ACMG combining-rule fidelity | `python -m eval.erepo_reconstruction <erepo.tab>` |
| Genome-wide partition coverage and reuse | `python -m eval.erepo_genomewide <erepo.tab>` |
| Per-gene frequency threshold triangulation | `python -m eval.erepo_thresholds <erepo.tab>` |
| CDC CHAMP and CHBMP recall | `python -m eval.champ_chbmp_benchmark` |
| Variant arm, discrimination and out-of-fold calibration | `python -m bench.phase_r_variant` |
| Per-criterion agreement and the circularity exhibit | `python -m bench.track1b_erepo_headtohead` |
| Gene-term sensitivity sweep | `python -m bench.phase_r_gene_term_sensitivity` |
| Safety interlock and abstention | `python -m bench.track3_trustworthiness` |
| Curated-case diagnosis | `python -m eval.curated_case_benchmark` |
| Phenotype-blind baselines and strata | `python -m eval.gene_only_baseline` |
| HPO representability of the discriminating features | `python -m eval.phenotype_tool_comparison` |
| External comparison against LIRICAL | `python -m eval.lirical_arm` |
| Coupling proof of concept | `python -m eval.coupling_poc` |

## Verification

```bash
make ci        # lint and tests against the tracked file set, as continuous integration sees it
```

Several tests are regression guards rather than unit tests. They fail if a value quoted in the
documentation drifts from its harness output, if a discrimination likelihood ratio loses its
source PMID, if an ACMG criterion changes factor, or if a claim the results do not support
appears in a document.

## Not reproducible from open data

The following require controlled-access data and are not part of any reported result. Their
evaluation is pre-registered in `DISCERN_OSF_PreRegistration_v1.md` and will be reported
regardless of outcome.

- The disease-variant coupling endpoint on real paired phenotype-genotype data.
- Cohort diagnosis, misdiagnosis rescue, and per-patient VUS rate.
- The randomized reader study.

## Release

A Zenodo DOI is minted on each tagged GitHub release. No patient-level data appears in any
artifact. Container images for reproducible deployment are built by the release workflow.
