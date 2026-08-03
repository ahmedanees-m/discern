"""Assemble the manuscript package and the Zenodo supporting-data deposit.

Collects from three places - this repository, the author's Drive (manuscript and planning
documents), and the VM (the LIRICAL run's raw per-case outputs) - and writes a self-describing
tree with a README, a file map, a checksum manifest, and licences.

The redistribution rule is enforced here, not left to judgement: third-party databases are never
copied. What ships is the identifiers, the derived outputs, the configurations, the scripts, and a
manifest giving each source's URL, version and checksum so anyone can re-pull and re-run. Two hard
exclusions: no patient-level data, and no copyrighted case-report text.

Run:  python -m deploy.build_submission_package [--out DIR] [--skip-vm]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DRIVE = os.path.join(os.path.dirname(ROOT), "v3.1", "Current status_DISCERN")

VERSION = "v0.1.0"
DEPOSIT = f"discern-paper1-supplement-{VERSION}"


def _sh(*cmd):
    return subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT).stdout.strip()


def md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def copy(src, dst_dir, name=None):
    if not os.path.exists(src):
        return None
    os.makedirs(dst_dir, exist_ok=True)
    dst = os.path.join(dst_dir, name or os.path.basename(src))
    if os.path.isdir(src):
        shutil.copytree(src, dst, dirs_exist_ok=True)
    else:
        shutil.copy2(src, dst)
    return dst


def fetch_vm_artifacts(dest):
    """Pull the LIRICAL raw per-case outputs. Paramiko only; no rclone, no host installs."""
    try:
        import paramiko
    except ImportError:
        print("  ! paramiko unavailable - skipping VM collection")
        return None
    host = os.environ.get("VM_HOST", "10.30.158.35")
    user = os.environ.get("VM_USER", "anees_22phd0670")
    pw = os.environ.get("VM_PASSWORD")
    if not pw:
        print("  ! VM_PASSWORD not set - skipping VM collection")
        return None
    remote = f"/home/{user}/lirical_run/lirical_raw_outputs.tar.gz"
    os.makedirs(dest, exist_ok=True)
    local = os.path.join(dest, "lirical_raw_outputs.tar.gz")
    # Retry: SFTP over this link occasionally short-reads, and one optional artefact failing must
    # not abort the whole package build. Size is verified against the remote stat, not assumed.
    for attempt in (1, 2, 3):
        t = None
        try:
            # Constructing the Transport already opens the socket, so it must be inside the guard:
            # an unreachable VM is a warning, never a reason to abandon the package.
            t = paramiko.Transport((host, 22))
            t.connect(username=user, password=pw)
            sftp = paramiko.SFTPClient.from_transport(t)
            expect = sftp.stat(remote).st_size
            # Chunked read rather than sftp.get(): the prefetch path short-reads reproducibly on
            # this link, and a truncated archive is worse than none.
            with sftp.open(remote, "rb") as rf, open(local, "wb") as wf:
                rf.prefetch(expect)
                while True:
                    block = rf.read(32768)
                    if not block:
                        break
                    wf.write(block)
            got = os.path.getsize(local)
            if got == expect:
                print(f"  + VM: {os.path.basename(local)} ({got / 1e6:.1f} MB)")
                return local
            print(f"  ! VM transfer short ({got} != {expect}), attempt {attempt}")
        except Exception as exc:                                     # noqa: BLE001
            print(f"  ! VM fetch failed on attempt {attempt}: {exc}")
        finally:
            if t is not None:
                t.close()
    print("  ! VM artefact not collected; the deposit is complete without it (the derived LIRICAL "
          "rankings are already in benchmarks/diagnosis_arm/lirical_run/lirical_ranked.tsv)")
    if os.path.exists(local):
        os.remove(local)
    return None


# ------------------------------------------------------------------------------------------
DATA_MANIFESTS = {
    "erepo": {
        "source": "ClinGen Evidence Repository", "url": "https://erepo.clinicalgenome.org/evrepo/",
        "files": ["data/raw/erepo/erepo_classifications.tab"],
        "role": "expert-panel truth surface (primary)", "redistributed": False,
        "license_note": "ClinGen data are freely available; not redistributed here to keep the "
                        "deposit to derived outputs and identifiers."},
    "clinvar": {
        "source": "NCBI ClinVar", "url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/",
        "files": ["data/raw/clinvar/clinvar_20260503.vcf.gz",
                  "data/raw/clinvar/variant_summary.txt.gz"],
        "role": "supporting surface; deliberately blinded in the primary protocol",
        "redistributed": False, "license_note": "public domain; excluded for size."},
    "gnomad": {
        "source": "gnomAD", "url": "https://gnomad.broadinstitute.org/downloads",
        "files": ["data/raw/gnomad/gnomad.v4.0.constraint_metrics.tsv"],
        "role": "population allele frequency", "redistributed": False,
        "license_note": "gnomAD data are released without restriction; excluded for size."},
    "hpo": {
        "source": "Human Phenotype Ontology", "url": "https://hpo.jax.org/app/data/ontology",
        "files": ["data/raw/hpo/hp.json", "data/raw/hpo/phenotype.hpoa",
                  "data/raw/hpo/genes_to_phenotype.txt"],
        "role": "HPO terms and the gene-to-disease crosswalk used by the LIRICAL arm",
        "redistributed": False, "license_note": "HPO custom licence permits redistribution; "
                                                "excluded for size."},
    "revel_alphamissense": {
        "source": "REVEL and AlphaMissense", "url": "https://sites.google.com/site/revelgenomics/ ; "
                                                    "https://zenodo.org/records/8360242",
        "files": [], "role": "missense predictors consumed as PP3/BP4 and as comparators",
        "redistributed": False,
        "license_note": "AlphaMissense is CC BY 4.0; REVEL is free for academic use. Scores reach "
                        "this work through the committed GeneBe annotation cache."},
    "champ_chbmp": {
        "source": "CDC CHAMP and CHBMP", "url": "https://www.cdc.gov/hemophilia/php/mutation-project/",
        "files": [], "role": "independent F8/F9 disease-allele truth", "redistributed": False,
        "license_note": "US CDC public resource; the deposit carries the derived recall table only."},
    "phenopacket_store": {
        "source": "GA4GH Phenopacket Store",
        "url": "https://github.com/monarch-initiative/phenopacket-store",
        "files": ["eval/data/bleeding_subset.jsonl"], "role": "coupling proof-of-concept corpus",
        "redistributed": True,
        "license_note": "the extracted bleeding subset is small and derived; included."},
    "lirical": {
        "source": "LIRICAL", "url": "https://github.com/TheJacksonLaboratory/LIRICAL",
        "version": "2.4.1", "files": [], "role": "external phenotype-driven comparator",
        "redistributed": False,
        "license_note": "run in a container on the VM; its data bundle is fetched by the tool. Note "
                        "that its downloader uses FTP for mim2gene_medgen, which many networks "
                        "block; fetch that file over HTTPS from the same NCBI path if so."},
}


def build(out_root, skip_vm=False):
    pkg = os.path.abspath(out_root)
    dep = os.path.join(pkg, "zenodo_deposit", DEPOSIT)
    for d in (pkg, dep):
        os.makedirs(d, exist_ok=True)

    sha = _sh("git", "rev-parse", "HEAD")
    short = _sh("git", "rev-parse", "--short", "HEAD")
    dirty = bool(_sh("git", "status", "--porcelain"))
    print(f"repository at {short}{' (DIRTY)' if dirty else ''}")

    # ---- manuscript layer -------------------------------------------------------------
    print("manuscript layer")
    for fn in ("DISCERN_Paper1_Manuscript_v2.md", "DISCERN_Paper1_Submission_Package.md",
               "DISCERN_PreSubmission_Analysis_Plan_PhaseR.md", "DISCERN_Zenodo_Deposit_Manifest.md"):
        if copy(os.path.join(DRIVE, fn), pkg):
            print(f"  + {fn}")

    # ---- figures and tables, regenerated then placed in both layers --------------------
    print("figures and tables")
    from figures.make_figures import main as figmain
    from figures.make_tables import main as tabmain
    for fn, target in ((figmain, "figures"), (tabmain, "tables")):
        argv = sys.argv
        sys.argv = ["x", "--out", os.path.join(pkg, target)]
        try:
            fn()
        finally:
            sys.argv = argv
    for sub in ("figures", "tables"):
        target = os.path.join(dep, sub)
        if os.path.isdir(target):
            shutil.rmtree(target)
        shutil.copytree(os.path.join(pkg, sub), target)

    # ---- deposit: configs --------------------------------------------------------------
    print("deposit: configs")
    cfg = os.path.join(dep, "configs")
    copy(os.path.join(ROOT, "rules", "vcep", "specs"), os.path.join(cfg, "vcep_specs"))
    copy(os.path.join(ROOT, "diseases", "clusters"), os.path.join(cfg, "cluster_definitions"))
    from jointdx import factorgraph
    from rules.point_engine import BANDS
    from rules.vcep.partition import FACTOR_OF
    with open(os.path.join(cfg, "partition_map.json"), "w", encoding="utf-8") as fh:
        json.dump({"description": "ACMG criterion -> owning factor. Exclusive by construction; a "
                                  "unit test asserts the variant marginal is invariant to re-adding "
                                  "codes owned elsewhere.",
                   "code_to_factor": FACTOR_OF}, fh, indent=2)
    with open(os.path.join(cfg, "locked_config.json"), "w", encoding="utf-8") as fh:
        json.dump({"description": "the frozen configuration behind every reported number",
                   "commit": sha, "version": VERSION,
                   "gene_term": {"on_gene": factorgraph.ON_GENE, "off_gene": factorgraph.OFF_GENE,
                                 "likelihood_ratio": round(factorgraph.ON_GENE / factorgraph.OFF_GENE, 3),
                                 "incidental_pp": factorgraph.INCIDENTAL_PP},
                   "acmg_bands": [[t, c.name] for t, c in BANDS],
                   "calibration": {"method": "isotonic", "folds": 5, "stratified": True,
                                   "shuffle": True, "seed": 0, "fit_on": "training folds only",
                                   "reported_on": "held-out folds only"},
                   "bootstrap": {"resamples": 1000, "interval": "percentile 95%"},
                   "timesplit_after": "2021-05-01"}, fh, indent=2)

    # ---- deposit: benchmarks -----------------------------------------------------------
    print("deposit: benchmarks")
    b = os.path.join(dep, "benchmarks")
    for src, sub in (
            ("bench/phase_r_variant_metrics.json", "variant_arm"),
            ("bench/track1b_erepo_metrics.json", "variant_arm"),
            ("bench/track1_metrics.json", "variant_arm"),
            ("bench/track1_variant_headtohead.csv", "variant_arm"),
            ("bench/data/genebe_erepo.jsonl", "variant_arm"),
            ("bench/data/genebe_h4set.jsonl", "variant_arm"),
            ("bench/data/erepo_bleeding.tsv", "variant_arm"),
            ("bench/track3_metrics.json", "trustworthiness"),
            ("bench/track3_risk_coverage.csv", "trustworthiness"),
            ("eval/gene_only_baseline.json", "diagnosis_arm"),
            ("eval/phenotype_tool_comparison.json", "diagnosis_arm"),
            ("eval/lirical_arm.json", "diagnosis_arm"),
            ("eval/cases/curated_cases.yaml", "diagnosis_arm"),
            ("eval/hpo_feature_crosswalk.yaml", "diagnosis_arm"),
            ("eval/data/lirical_cases.tsv", "diagnosis_arm/lirical_run"),
            ("eval/data/lirical_ranked.tsv", "diagnosis_arm/lirical_run"),
            ("eval/data/phenopackets", "diagnosis_arm/phenopackets"),
            ("bench/phase_r_gene_term_sensitivity.json", "gene_term"),
            ("eval/data/bleeding_subset.jsonl", "coupling_poc"),
            ("eval/champ_chbmp_ps1_revel.tsv", "champ_chbmp")):
        if copy(os.path.join(ROOT, src), os.path.join(b, sub)):
            print(f"  + {src}")

    if not skip_vm:
        print("deposit: VM artifacts")
        fetch_vm_artifacts(os.path.join(b, "diagnosis_arm", "lirical_run"))

    # ---- deposit: harness scripts (so every number regenerates) ------------------------
    print("deposit: harnesses")
    code = os.path.join(dep, "code", "harnesses")
    for src in ("bench/phase_r_variant.py", "bench/phase_r_gene_term_sensitivity.py",
                "bench/track1_variant_headtohead.py", "bench/track1b_erepo_headtohead.py",
                "bench/track3_trustworthiness.py", "eval/gene_only_baseline.py",
                "eval/phenotype_tool_comparison.py", "eval/lirical_arm.py",
                "eval/erepo_genomewide.py", "eval/curated_case_benchmark.py",
                "eval/coupling_poc.py", "rules/vcep/partition.py", "jointdx/factorgraph.py",
                "figures/make_figures.py", "figures/make_tables.py"):
        copy(os.path.join(ROOT, src), code)

    # ---- deposit: docs ------------------------------------------------------------------
    print("deposit: docs")
    d = os.path.join(dep, "docs")
    for fn in sorted(os.listdir(os.path.join(ROOT, "docs"))):
        if fn.endswith(".md"):
            copy(os.path.join(ROOT, "docs", fn), d)
    copy(os.path.join(ROOT, "docs", "DISCERN_OSF_PreRegistration_v1.md"),
         os.path.join(dep, "preregistration"))

    # ---- deposit: data manifests ---------------------------------------------------------
    print("deposit: data manifests")
    dm = os.path.join(dep, "data_manifests")
    os.makedirs(dm, exist_ok=True)
    for key, spec in DATA_MANIFESTS.items():
        entry = dict(spec)
        entry["files"] = []
        for rel in spec["files"]:
            p = os.path.join(ROOT, rel)
            if os.path.exists(p):
                entry["files"].append({"path": rel, "bytes": os.path.getsize(p), "md5": md5(p)})
            else:
                entry["files"].append({"path": rel, "note": "not present in this checkout"})
        with open(os.path.join(dm, f"{key}_manifest.json"), "w", encoding="utf-8") as fh:
            json.dump(entry, fh, indent=2)
    write_readmes(pkg, dep, sha, short)
    return pkg, dep, sha, short, dirty


DEPOSIT_README = """# DISCERN Paper 1 - supporting data, configurations, and results

**Version {version}** | Code commit `{sha}` | Repository <https://github.com/ahmedanees-m/discern>

This is Record B of a two-record deposit. Record A is the code, archived automatically from a
GitHub release. This record holds everything needed to check or regenerate the numbers in the
manuscript: the derived benchmark outputs, the frozen configuration, the variant and case
identifiers, the figures and tables, the harness scripts, and a manifest of every third-party
source with its version so it can be re-pulled.

## What is deliberately not here

**No third-party databases.** ClinVar, gnomAD, the ClinGen Evidence Repository, REVEL,
AlphaMissense, the HPO release files and the CDC CHAMP/CHBMP lists are referenced by manifest, not
redistributed. `data_manifests/` gives each one a URL, a version, and md5 checksums of the exact
files used, so the pinned inputs can be recovered. This keeps the deposit small and license-clean
while remaining reproducible.

**No patient-level data.** None was used in this work. A local consanguinity-enriched cohort exists
for Paper 2 under controlled access and appears nowhere here (Gate G7).

**No copyrighted case-report text.** The 42-case diagnosis benchmark ships as PMIDs, extracted HPO
terms, gene, and expected diagnosis only.

## File map

```
{tree}
```

## Where each manuscript claim is evidenced

| Claim | File |
|---|---|
| Variant AUROC and CIs, DeLong test | `benchmarks/variant_arm/phase_r_variant_metrics.json` |
| Calibration, raw and calibrated, all tools | `benchmarks/variant_arm/phase_r_variant_metrics.json` |
| The intrinsic-evidence ceiling and its attribution | same file, `ceiling_attribution` |
| Per-criterion kappa, GeneBe circularity | `benchmarks/variant_arm/track1b_erepo_metrics.json` |
| Safety interlock, diagnosis calibration | `benchmarks/trustworthiness/track3_metrics.json` |
| Diagnosis vs its baselines, strata, McNemar | `benchmarks/diagnosis_arm/gene_only_baseline.json` |
| HPO representability (13 of 48 features) | `benchmarks/diagnosis_arm/phenotype_tool_comparison.json` |
| LIRICAL arms and the fix-invariance assertion | `benchmarks/diagnosis_arm/lirical_arm.json` |
| LIRICAL raw per-case output | `benchmarks/diagnosis_arm/lirical_run/` |
| Gene-term sensitivity sweep | `benchmarks/gene_term/phase_r_gene_term_sensitivity.json` |
| The 33.2% inflation figure | regenerate with `code/harnesses/erepo_genomewide.py` |
| The evidence partition itself | `configs/partition_map.json` |
| Every discrimination likelihood ratio with its PMID | `tables/tableS3_cluster_likelihood_ratios.csv` |

## Reproducing

1. Clone the repository at the commit above, or install the code record.
2. Re-pull the third-party sources at the versions in `data_manifests/` and verify the md5s.
3. Run the harnesses in `code/harnesses/`. Each writes the JSON it is named for.
4. Regenerate the display items: `python -m figures.make_figures` and `python -m figures.make_tables`.
5. `pytest` locks every headline number; `tests/test_phase_r.py` and `tests/test_docs_claims.py`
   fail if a value drifts or a retired claim is reasserted.

## Read this first

`docs/DISCERN_PhaseR_Results_v1.md` is the authoritative record of what the pre-submission analysis
found, **including the results that cost the paper claims**: the calibration-uniqueness claim is
retired, the risk-coverage claim is withdrawn, the diagnosis accuracy number is in-sample and
subordinate to a 93% phenotype-blind baseline, and a factual error about the evidence partition was
corrected. Where an older document in `docs/` disagrees with it, the Phase R record wins; those
documents carry a banner saying so.

## Licences

Derived results and documentation: CC BY 4.0 (`LICENSE-data.txt`). Scripts: MIT
(`LICENSE-code.txt`). Third-party sources retain their own terms and are not redistributed.
"""

CITATION_CFF = """cff-version: 1.2.0
title: "DISCERN Paper 1 - supporting data, configurations, and results"
message: "If you use this deposit, please cite it and the associated code record."
type: dataset
version: "{version}"
authors:
  - family-names: "Mahaboob Ali"
    given-names: "Anees Ahmed"
    alias: "ahmedanees-m"
repository-code: "https://github.com/ahmedanees-m/discern"
keywords:
  - variant interpretation
  - ACMG/AMP
  - variant of uncertain significance
  - inherited bleeding disorders
  - inherited platelet disorders
  - differential diagnosis
  - calibration
  - ClinGen
  - reproducibility
license: CC-BY-4.0
abstract: >-
  Derived benchmark outputs, frozen configurations, variant and case identifiers, figures, tables,
  and harness scripts supporting DISCERN Paper 1. Third-party databases are referenced by versioned
  manifest rather than redistributed. Contains no patient-level data and no reproduced article text.
"""

LICENSE_DATA = """Creative Commons Attribution 4.0 International (CC BY 4.0)

Applies to the derived results, configurations, figures, tables and documentation in this deposit.

You are free to share and adapt this material for any purpose, provided you give appropriate
credit, link to the licence, and indicate if changes were made.

Full text: https://creativecommons.org/licenses/by/4.0/legalcode

NOT COVERED BY THIS LICENCE
Third-party databases referenced in data_manifests/ (ClinGen eRepo, ClinVar, gnomAD, REVEL,
AlphaMissense, the Human Phenotype Ontology, CDC CHAMP/CHBMP, GA4GH Phenopacket Store) are not
redistributed here and retain their own terms. Re-pull them from the URLs and versions given.
"""

CODE_README = """# Code record

The software is archived as a separate Zenodo record, minted from a tagged GitHub release.

- Repository: https://github.com/ahmedanees-m/discern
- Release tag: `{version}-paper1`
- Commit: `{sha}`
- Code DOI: [insert once the release is cut and Zenodo has snapshotted it]

The repository is not duplicated here. `harnesses/` carries copies of the scripts that produce the
benchmark outputs in this deposit, so the derivation is legible without a clone, but the release is
the citable artefact.
"""


def _tree(root, prefix="", depth=0, max_depth=2):
    out = []
    entries = sorted(os.listdir(root))
    dirs = [e for e in entries if os.path.isdir(os.path.join(root, e))]
    files = [e for e in entries if not os.path.isdir(os.path.join(root, e))]
    for d in dirs:
        n = sum(len(f) for _, _, f in os.walk(os.path.join(root, d)))
        out.append(f"{prefix}{d}/  ({n} files)")
        if depth < max_depth:
            out += _tree(os.path.join(root, d), prefix + "    ", depth + 1, max_depth)
    for f in files:
        out.append(f"{prefix}{f}")
    return out


def write_readmes(pkg, dep, sha, short):
    with open(os.path.join(dep, "LICENSE-data.txt"), "w", encoding="utf-8") as fh:
        fh.write(LICENSE_DATA)
    src_lic = os.path.join(ROOT, "LICENSE")
    if os.path.exists(src_lic):
        shutil.copy2(src_lic, os.path.join(dep, "LICENSE-code.txt"))
    with open(os.path.join(dep, "CITATION.cff"), "w", encoding="utf-8") as fh:
        fh.write(CITATION_CFF.format(version=VERSION))
    os.makedirs(os.path.join(dep, "code"), exist_ok=True)
    with open(os.path.join(dep, "code", "README.md"), "w", encoding="utf-8") as fh:
        fh.write(CODE_README.format(version=VERSION, sha=sha))
    with open(os.path.join(dep, "README.md"), "w", encoding="utf-8") as fh:
        fh.write(DEPOSIT_README.format(version=VERSION, sha=sha,
                                       tree="\n".join(_tree(dep))))
    with open(os.path.join(pkg, "README.md"), "w", encoding="utf-8") as fh:
        fh.write(PACKAGE_README.format(version=VERSION, sha=sha, short=short,
                                       tree="\n".join(_tree(pkg, max_depth=1))))


PACKAGE_README = """# DISCERN Paper 1 - manuscript files

**Version {version}** | Code commit `{sha}` | Repository <https://github.com/ahmedanees-m/discern>

Everything needed to submit DISCERN Paper 1, assembled from this repository, the planning documents,
and the VM that ran the external comparator. Regenerate the whole tree with:

```
python -m deploy.build_submission_package
```

## Contents

```
{tree}
```

| Item | What it is |
|---|---|
| `DISCERN_Paper1_Manuscript_v2.md` | the manuscript. All 21 references verified; the only placeholders left are the by-line and the three DOIs |
| `DISCERN_Paper1_Submission_Package.md` | the submission plan: corrections, consistency checks, figure and table specifications, journal-required statements, submission sequence |
| `DISCERN_PreSubmission_Analysis_Plan_PhaseR.md` | the Phase R plan with its dispositions filled in |
| `DISCERN_Zenodo_Deposit_Manifest.md` | the original deposit design |
| `figures/` | Figures 1-6 and S1-S4, vector PDF plus 300 dpi PNG, generated from committed JSON |
| `tables/` | Tables 1-3 and S1-S8 as CSV, generated from the same sources |
| `zenodo_deposit/` | the supporting-data record, ready to upload |

## Version

This package is **{version}**, not v1.0.0. The work is complete and internally consistent, but the
by-line, the three DOIs, and the journal-required statements are outstanding, and no release has
been tagged. v1.0.0 is the right label once the paper is submitted and the identifiers exist.

## Before submitting

Ordered, because several steps depend on the one before:

1. Fill the by-line: authors, order, corresponding author, ORCIDs, exact unit name.
2. Write the journal-required statements listed in the submission package section 6 (ethics, consent,
   competing interests, funding, CRediT contributions, acknowledgements, abbreviations).
3. Time-stamp the OSF pre-registration -> **OSF DOI**. This also clears Gate G12 and must precede any
   cohort analysis.
4. Enable the repository in Zenodo, then cut GitHub release `{version}-paper1` -> **code DOI**.
   Order matters: Zenodo only snapshots releases created after the integration is switched on.
5. Upload `zenodo_deposit/` as a Zenodo Dataset record -> **data DOI**.
6. Insert all three DOIs into the manuscript's Data and Code Availability section.
7. Re-run the reproducibility checklist from a clean clone at the tag.
8. Post the preprint, then submit.

## A note on scope

The manuscript is deliberate about what it does and does not claim. Two claims were retired during
pre-submission analysis and one was withdrawn; the diagnosis accuracy figure is in-sample and is
reported against a phenotype-blind baseline it does not significantly beat. That framing is the
paper's strongest asset and should survive editing. `docs/DISCERN_PhaseR_Results_v1.md` inside the
deposit records every disposition, including the unfavourable ones.
"""


def write_manifest(dep):
    rows = []
    for base, _dirs, files in os.walk(dep):
        for f in sorted(files):
            if f == "MANIFEST.md":
                continue
            p = os.path.join(base, f)
            rows.append((os.path.relpath(p, dep).replace("\\", "/"), os.path.getsize(p), md5(p)))
    rows.sort()
    total = sum(r[1] for r in rows)
    with open(os.path.join(dep, "MANIFEST.md"), "w", encoding="utf-8") as fh:
        fh.write("# MANIFEST\n\nEvery file in this deposit with its size and md5 checksum.\n\n")
        fh.write(f"{len(rows)} files, {total / 1e6:.1f} MB total.\n\n")
        fh.write("| file | bytes | md5 |\n|---|---|---|\n")
        for rel, size, h in rows:
            fh.write(f"| `{rel}` | {size} | `{h}` |\n")
    return len(rows), total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(ROOT), "manuscript_files"))
    ap.add_argument("--skip-vm", action="store_true")
    args = ap.parse_args()
    pkg, dep, sha, short, dirty = build(args.out, args.skip_vm)
    n, total = write_manifest(dep)
    print(f"\npackage : {pkg}")
    print(f"deposit : {dep}")
    print(f"commit  : {sha}{' (DIRTY - commit before minting a DOI)' if dirty else ''}")
    print(f"manifest: {n} files, {total / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
