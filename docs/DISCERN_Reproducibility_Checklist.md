# DISCERN v3.1 - Reproducibility Checklist (Track D2)

Every public number is reproducible from a clean clone + the named open datasets. Heavy runs
execute on the VM inside Docker (Gate G8: no host installs; secrets from env, never in git).

## Environment
```
git clone https://github.com/ahmedanees-m/discern && cd discern
conda env create -f environment.yml        # or: pip install -e ".[dev]"
make test                                   # ruff + pytest (135 tests, CI-green)
```

## Datasets (open; pull fresh at run time - versions drift)
| Dataset | Source | Version used |
|---|---|---|
| ClinGen eRepo | erepo.clinicalgenome.org (download .tab) | 2026-05 export, 170 genes / 12,240 variants |
| ClinVar variant_summary | NCBI FTP `tab_delimited/variant_summary.txt.gz` | 2026-05-04 |
| gnomAD | gnomad.broadinstitute.org | v4.1 callset / v4.1.1 annotations (per-variant AFs needed for the freq cross-check - not yet on the VM) |
| Predictors | REVEL (dbNSFP), Pangolin (tkzeng v1.0.1 / Invitae v1.4.x), AlphaMissense (CC BY 4.0) | pin at run time |

## Reproduce the public numbers
| Result | Command | Expected |
|---|---|---|
| Genome-wide partition (H1/H2) | `python3 -m eval.erepo_genomewide <erepo.tab>` | 100% coverage; 33.2% inflation-prevented (CI 32.4-34.1) |
| ACMG combining fidelity | `python3 -m eval.erepo_reconstruction <erepo.tab>` | 93.0% exact / 100% within-1-bin |
| ClinVar concordance (H3) | `python3 -m eval.clinvar_concordance <erepo.tab> <variant_summary.gz>` | 62.4% exact / 92.8% within-1-bin |
| Variant calibration (H5) | `python3 -m eval.variant_calibration <erepo.tab> <variant_summary.gz>` | isotonic ECE 0.008 / Brier 0.0073 |
| CSpec threshold triangulation | `python3 -m eval.erepo_thresholds <erepo.tab>` | per-gene BA1/BS1/PM2 + PM2_Supporting |
| Curated case diagnosis (B4) | `python3 -m eval.curated_case_benchmark` | Top-1 80% / Top-3 100% (n=10) |
| Coupling sanity (synthetic) | `python3 -m eval.synthetic_coupling_harness` | runs; circularity guard True (NOT a validation) |
| CI guards | `pytest tests/test_erepo_eval.py tests/test_v31_clusters.py` | partition coverage 1.0; every LR PMID-sourced |

## Data-gated / external (cannot be reproduced from open data)
- gnomAD per-variant frequency cross-check (A1) - needs per-variant AFs pulled to the VM.
- H4 novel-variant AUC vs InterVar - needs real REVEL/Pangolin/AlphaMissense per variant.
- The coupling validation, cohort diagnosis, misdiagnosis rescue, reader study, per-patient VUS
  rate - **DAC/IRB gated** (BRIDGE-BPD EGAS00001001172; ITP cohort; a controlled-access paired cohort).
  Pre-registered (`DISCERN_OSF_PreRegistration_v1.md`) before analysis; reported regardless of outcome.

## Release
Zenodo DOI minted on the first GitHub Release; no patient data in any artifact (G7); Docker
images for reproducible deployment (`deploy/compose.vm.yml`).
