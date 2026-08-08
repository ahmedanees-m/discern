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
PKG_HOME = os.path.join(os.path.dirname(ROOT), "manuscript_files")

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


# What the VM holds that nothing else does: the raw per-case LIRICAL output, and the script that
# produced it. The 381 MB LIRICAL distribution and its data bundle stay on the VM - they are the
# tool, not our result, and data_manifests/ records the version needed to re-pull them.
VM_ARTIFACTS = [
    ("lirical_raw_outputs.tar.gz", "lirical_raw_outputs.tar.gz"),
    ("cmds.sh", "lirical_run_commands.sh"),
]


def fetch_vm_artifacts(dest):
    """Pull the LIRICAL raw outputs and run script. Paramiko only; no rclone, no host installs."""
    try:
        import paramiko
    except ImportError:
        print("  ! paramiko unavailable - skipping VM collection")
        return []
    host = os.environ.get("VM_HOST", "10.30.158.35")
    user = os.environ.get("VM_USER", "anees_22phd0670")
    pw = os.environ.get("VM_PASSWORD")
    if not pw:
        print("  ! VM_PASSWORD not set - skipping VM collection")
        return []
    os.makedirs(dest, exist_ok=True)
    got = []
    # Retry: SFTP over this link occasionally short-reads, and one optional artefact failing must
    # not abort the whole package build. Size is verified against the remote stat, not assumed.
    for name, local_name in VM_ARTIFACTS:
        remote = f"/home/{user}/lirical_run/{name}"
        local = os.path.join(dest, local_name)
        for attempt in (1, 2, 3):
            t = None
            try:
                # Constructing the Transport already opens the socket, so it must be inside the
                # guard: an unreachable VM is a warning, never a reason to abandon the package.
                t = paramiko.Transport((host, 22))
                t.connect(username=user, password=pw)
                sftp = paramiko.SFTPClient.from_transport(t)
                expect = sftp.stat(remote).st_size
                # Chunked read rather than sftp.get(): the prefetch path short-reads reproducibly
                # on this link, and a truncated archive is worse than none.
                with sftp.open(remote, "rb") as rf, open(local, "wb") as wf:
                    rf.prefetch(expect)
                    while True:
                        block = rf.read(32768)
                        if not block:
                            break
                        wf.write(block)
                size = os.path.getsize(local)
                if size == expect:
                    print(f"  + VM: {local_name} ({size / 1e6:.1f} MB)")
                    got.append(local)
                    break
                print(f"  ! VM transfer short ({size} != {expect}), attempt {attempt}")
            except Exception as exc:                                     # noqa: BLE001
                print(f"  ! VM fetch of {name} failed on attempt {attempt}: {exc}")
            finally:
                if t is not None:
                    t.close()
        else:
            print(f"  ! VM artefact {name} not collected; the deposit is complete without it (the "
                  f"derived LIRICAL rankings are already in "
                  f"benchmarks/diagnosis_arm/lirical_run/lirical_ranked.tsv)")
            if os.path.exists(local):
                os.remove(local)
    return got


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
    # Two sources. The planning documents live on the author's Drive; the submission documents were
    # written in the package itself and have no other home, so a build to a different --out must
    # collect them from the canonical package rather than silently produce an empty manuscript
    # directory.
    print("manuscript layer")
    for fn in ("DISCERN_Paper1_Manuscript_v2.md", "DISCERN_Paper1_Submission_Package.md",
               "DISCERN_PreSubmission_Analysis_Plan_PhaseR.md", "DISCERN_Zenodo_Deposit_Manifest.md"):
        if copy(os.path.join(DRIVE, fn), pkg):
            print(f"  + {fn}")
    for fn in MANUSCRIPT_DOCS + ARCHIVE_DOCS:
        if os.path.exists(os.path.join(pkg, fn)):
            continue
        for cand in (os.path.join(PKG_HOME, fn),
                     os.path.join(PKG_HOME, "01_manuscript", fn),
                     os.path.join(PKG_HOME, "archive", fn)):
            if copy(cand, pkg):
                print(f"  + {fn}")
                break
        else:
            print(f"  ! {fn} not found in {PKG_HOME}")

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
    # Figures, tables and the compiled supplement are display items, published with the article.
    # None is deposited: every one regenerates from the JSON and YAML in this deposit using the
    # four generator scripts in code/harnesses/. The rule is to deposit what cannot be regenerated
    # from the deposit, so any that survive an earlier build are removed here.
    for sub in ("figures", "tables", "supplemental"):
        stale = os.path.join(dep, sub)
        if os.path.isdir(stale):
            try:
                shutil.rmtree(stale)
                print(f"  - removed {sub}/ from the deposit (display items, not deposited)")
            except OSError as exc:
                print(f"  ! {sub}/ is still in the deposit and could not be removed "
                      f"({exc.strerror}); remove it before uploading")

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
            ("bench/calibration_folds.json", "variant_arm"),
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
                "core/stats.py", "figures/make_excel_tables.py", "figures/make_document_s1.py",
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
    # The compiled supplement and the Excel display items travel with the record, and also sit at
    # the package root: these are the files the journal is uploaded, not just archived copies.
    dest = os.path.join(pkg, "supplemental")
    copy(os.path.join(ROOT, "figures", "Document_S1_Supplemental_Information.pdf"), dest)
    copy(os.path.join(ROOT, "figures", "excel"), dest, name="excel_tables")
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
    # Lay the package out before writing the READMEs, since the package README embeds a file map
    # of the tree as it will actually be delivered.
    lay_out_submission(pkg)
    write_readmes(pkg, dep, sha, short)
    return pkg, dep, sha, short, dirty


DEPOSIT_README = """# DISCERN Paper 1 - supporting data, configurations, and results

**Version {version}** | Code commit `{sha}` | Repository <https://github.com/ahmedanees-m/discern>

This is Record B of a two-record deposit. Record A is the code, archived automatically from a
GitHub release. This record holds everything needed to check or regenerate the numbers in the
manuscript: the derived benchmark outputs, the frozen configuration, the variant and case
identifiers, the harness scripts, and a manifest of every third-party source with its version
so it can be re-pulled. It deliberately carries no figures or tables; see below for why.

## What is deliberately not here

**No third-party databases.** ClinVar, gnomAD, the ClinGen Evidence Repository, REVEL,
AlphaMissense, the HPO release files and the CDC CHAMP/CHBMP lists are referenced by manifest, not
redistributed. `data_manifests/` gives each one a URL, a version, and md5 checksums of the exact
files used, so the pinned inputs can be recovered. This keeps the deposit small and license-clean
while remaining reproducible.

**No patient-level data.** None was used in this work. A local consanguinity-enriched cohort exists
for Paper 2 under controlled access and appears nowhere here.

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
| The out-of-fold assignment behind the calibration | `benchmarks/variant_arm/calibration_folds.json` |
| The intrinsic-evidence ceiling and its attribution | same file, `ceiling_attribution` |
| Per-criterion kappa, GeneBe circularity | `benchmarks/variant_arm/track1b_erepo_metrics.json` |
| Safety interlock, diagnosis calibration | `benchmarks/trustworthiness/track3_metrics.json` |
| Diagnosis vs its baselines, strata, McNemar | `benchmarks/diagnosis_arm/gene_only_baseline.json` |
| HPO representability (13 of 48 features) | `benchmarks/diagnosis_arm/phenotype_tool_comparison.json` |
| LIRICAL arms and the fix-invariance assertion | `benchmarks/diagnosis_arm/lirical_arm.json` |
| LIRICAL raw per-case output | `benchmarks/diagnosis_arm/lirical_run/` |
| The exact LIRICAL commands, one per case | `benchmarks/diagnosis_arm/lirical_run/lirical_run_commands.sh` |
| Gene-term sensitivity sweep | `benchmarks/gene_term/phase_r_gene_term_sensitivity.json` |
| The 33.2% inflation figure | regenerate with `code/harnesses/erepo_genomewide.py` |
| The evidence partition itself | `configs/partition_map.json` |
| Every discrimination likelihood ratio with its PMID | `configs/cluster_definitions/` |

## Why there are no figures or tables here

The rule applied to this deposit is to include what cannot be regenerated from it, and to exclude
what can. Every figure and every table in the article is a display item rendered from the JSON and
YAML files in this record, by `make_figures.py` and `make_tables.py`, both of which are included in
`code/harnesses/`. Depositing the rendered output as well would add copies that any reader can
reproduce in one command, and would leave two versions to fall out of step. The figures and tables
are published with the article.

The benchmark outputs are a different case and are all present. They cannot be regenerated without
re-pulling roughly a hundred gigabytes of third-party databases, so they are the surface on which
the article's numbers can actually be checked.

## Reproducing

1. Clone the repository at the commit above, or install the code record.
2. Re-pull the third-party sources at the versions in `data_manifests/` and verify the md5s.
3. Run the harnesses in `code/harnesses/`. Each writes the JSON it is named for.
4. Regenerate the display items from a clone of the code record: `python -m figures.make_figures`,
   `python -m figures.make_tables`, `python -m figures.make_excel_tables` and
   `python -m figures.make_document_s1`. Copies of those four scripts are in `code/harnesses/` here
   so they can be read without a clone, but run them from the code record, not from this directory,
   since they import the package.
5. Re-run the test suite **from the code record**. `tests/` belongs to Record A and is deliberately
   not duplicated here. `tests/test_phase_r.py` and `tests/test_docs_claims.py` fail if a reported
   value drifts or a retired claim is reasserted.

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
identifiers:
  - type: doi
    value: "10.17605/OSF.IO/NUZT9"
    description: "Pre-registered analysis plan for the disease-variant coupling"
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
    try:
        entries = sorted(os.listdir(root))
    except OSError as exc:
        # A directory that cannot be listed is worth recording in the tree, but it must not take
        # the build down with it: the file map is descriptive, not load-bearing.
        return [f"{prefix}[unreadable: {exc.strerror}]"]
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
    src_lic = os.path.join(ROOT, "LICENSE")
    if os.path.exists(src_lic):
        shutil.copy2(src_lic, os.path.join(dep, "LICENSE-code.txt"))
    os.makedirs(os.path.join(dep, "code"), exist_ok=True)

    # Render every file before opening any of them. Opening for write truncates, so rendering
    # inside the `with` would leave an empty file behind if the substitution raised - which is
    # exactly how a zero-byte README once shipped in the deposit.
    written = {
        os.path.join(dep, "LICENSE-data.txt"): LICENSE_DATA,
        os.path.join(dep, "CITATION.cff"): CITATION_CFF.format(version=VERSION),
        os.path.join(dep, "code", "README.md"): CODE_README.format(version=VERSION, sha=sha),
        os.path.join(dep, "README.md"): DEPOSIT_README.format(
            version=VERSION, sha=sha, tree="\n".join(_tree(dep))),
        os.path.join(pkg, "README.md"): PACKAGE_README.format(
            version=VERSION, sha=sha, short=short, tree="\n".join(_tree(pkg, max_depth=1))),
    }
    for path, text in written.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)


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

The numbered directories are the upload order.

| Item | What it is |
|---|---|
| `01_manuscript/` | `DISCERN_Paper1_BMC_HumanGenomics_SUBMISSION.md`, the manuscript to submit, in Human Genomics format (structured abstract, Background/Methods/Results/Discussion/Conclusions, the mandatory Declarations block, BMC Vancouver references with web links numbered), and its cover letter. The HGG Advances variant of the same science is in `archive/`: unstructured abstract, Material and methods, Web resources, every display item cited. Alongside it, the cover letter. The remaining placeholders in both are the by-line, the funding sanction number, and the three DOIs |
| `02_figures/` | `Figure1.pdf` to `Figure6.pdf`, vector, with 300 dpi PNG alongside each. Upload individually |
| `03_tables/` | `Table1.csv` to `Table3.csv`. These are embedded in the manuscript as formatted tables; the CSVs are the machine-readable source |
| `04_additional_files/` | exactly what is uploaded as supplemental material, named as the portal expects: `Additional_file_1.pdf` (the compiled supplement) through `Additional_file_5.xlsx`. Its README maps each file to the S-numbers used in the text. `components/` holds the individual figures and tables compiled into Additional file 1, for editing rather than upload |
| `05_graphical_abstract/` | `graphical_abstract.png`, 920 x 300 px landscape, required by BMC Human Genomics at submission |
| `zenodo_deposit/` | the supporting-data record, ready to upload to Zenodo. `MANIFEST.md` inside it lists every file with its md5 and a one-line description. It deliberately contains no figures or tables: those are display items published with the article, and each regenerates from the data in the deposit |
| `archive/` | superseded drafts and the planning documents, kept for provenance and not part of the submission: the v3 draft this was assembled from, the preceding generic-journal v2, the display-item and deposit audits, the pre-submission analysis plan, and the original deposit design |

## Version

This package is **{version}**, not v1.0.0. The work is complete and internally consistent, but the
by-line, the three DOIs, and the journal-required statements are outstanding, and no release has
been tagged. v1.0.0 is the right label once the paper is submitted and the identifiers exist.

## Before submitting

Target journal is **Human Genomics** (BMC, Springer Nature). Ordered, because several steps depend on the one
before, and every item here needs the author - none can be produced from the repository:

1. Fill the by-line: authors, order, corresponding author, ORCIDs, exact unit name, and the
   institutional PIN and email the submission form asks for.
2. Supply the VIT funding sanction number, then complete the remaining declarations: ethics,
   consent, competing interests, CRediT contributions, acknowledgements.
3. OSF pre-registration time-stamped: **10.17605/OSF.IO/NUZT9** (https://osf.io/nuzt9/). Gate G12 is cleared; no cohort
   analysis may precede this registration.
4. Enable the repository in Zenodo, then cut GitHub release `{version}-paper1` -> **code DOI**.
   Order matters: Zenodo only snapshots releases created after the integration is switched on.
5. Upload `zenodo_deposit/` as a Zenodo Dataset record -> **data DOI**.
6. Insert all three DOIs into the manuscript's Data and Code Availability section.
7. Confirm the reference list is in BMC/Springer Vancouver style (no italics, no issue
   numbers, abbreviated end pages, web links numbered with access dates), and run each
   reference through Retraction Watch.
8. Re-run the reproducibility checklist from a clean clone at the tag.
9. Nominate suggested reviewers, and request any APC waiver or discount **in the submission
   form itself**: BMC cannot consider requests made during review or after acceptance.
10. Post the preprint, then submit `DISCERN_Paper1_BMC_HumanGenomics_SUBMISSION.md` as the
    manuscript, `02_figures/` as the display items, `04_additional_files/` as the
    Additional files, and `05_graphical_abstract/` as the required graphical abstract.

## A note on scope

The manuscript is deliberate about what it does and does not claim. Two claims were retired during
pre-submission analysis and one was withdrawn; the diagnosis accuracy figure is in-sample and is
reported against a phenotype-blind baseline it does not significantly beat. That framing is the
paper's strongest asset and should survive editing. `docs/DISCERN_PhaseR_Results_v1.md` inside the
deposit records every disposition, including the unfavourable ones.
"""


# The submission tree. Display items are renamed to what a Cell Press portal expects, because
# "Figure3.pdf" is what the upload form asks for and "fig3_per_criterion_kappa.pdf" is not.
# Numbered directories put the upload in order: manuscript, then figures, then tables, then
# supplemental, then the data record.
MAIN_FIGURES = ["fig1_architecture_and_partition", "fig2_discrimination_and_calibration",
                "fig3_per_criterion_kappa", "fig4_intrinsic_ceiling",
                "fig5_clinvar_circularity", "fig6_diagnosis_baselines"]
SUPP_FIGURES = ["figS1_gene_term_sensitivity", "figS2_safety_matrix",
                "figS3_diagnosis_calibration", "figS4_champ_chbmp_recall",
                "figS5_worked_example"]
MAIN_TABLES = ["table1_variant_head_to_head", "table2_decision_quality_matched_coverage",
               "table3_diagnosis_vs_baselines"]
SUPP_TABLES = [f"tableS{i}_" for i in range(1, 10)]
MANUSCRIPT_DOCS = ["DISCERN_Paper1_BMC_HumanGenomics_SUBMISSION.md",
                   "DISCERN_Paper1_BMC_HumanGenomics_COVER_LETTER.md"]
# The HGG Advances variant of the same science, kept for provenance rather than submission.
HGG_DOCS = ["DISCERN_Paper1_HGG_Advances_SUBMISSION.md",
            "DISCERN_Paper1_HGG_Advances_COVER_LETTER.md"]

# BMC publishes supplementary material as separately uploaded, individually cited Additional
# files, so the compiled supplement becomes Additional file 1 and each large table follows it in
# citation order. The mapping to the S-numbering used in the text is written into the directory.
ADDITIONAL_FILES = [
    ("Document_S1_Supplemental_Information.pdf", "Additional_file_1.pdf",
     "Supplementary methods, Figures S1-S5, and Tables S1, S2, S4, S6, S8"),
    ("Table_S3.xlsx", "Additional_file_2.xlsx", "Table S3, cluster likelihood ratios"),
    ("Table_S5.xlsx", "Additional_file_3.xlsx", "Table S5, the curated published-case benchmark"),
    ("Table_S7.xlsx", "Additional_file_4.xlsx", "Table S7, third-party data sources"),
    ("Table_S9.xlsx", "Additional_file_5.xlsx", "Table S9, software and database versions"),
]
ARCHIVE_DOCS = HGG_DOCS + ["DISCERN_Paper1_HGG_Advances_v3.md", "DISCERN_Paper1_Manuscript_v2.md",
                "DISCERN_Display_Items_and_Deposit_Audit.md",
                "DISCERN_Final_PreSubmission_Report.md",
                "DISCERN_Paper1_Submission_Package.md",
                "DISCERN_PreSubmission_Analysis_Plan_PhaseR.md",
                "DISCERN_Zenodo_Deposit_Manifest.md"]


def lay_out_submission(pkg):
    """Arrange the package as the journal receives it, renaming display items to portal names."""
    d_ms = os.path.join(pkg, "01_manuscript")
    d_fig = os.path.join(pkg, "02_figures")
    d_tab = os.path.join(pkg, "03_tables")
    d_sup = os.path.join(pkg, "04_additional_files")
    d_ga = os.path.join(pkg, "05_graphical_abstract")
    d_sup_f = os.path.join(d_sup, "components", "figures")
    d_sup_t = os.path.join(d_sup, "components", "tables")
    d_arc = os.path.join(pkg, "archive")
    for d in (d_ms, d_fig, d_tab, d_sup_f, d_sup_t, d_ga, d_arc):
        os.makedirs(d, exist_ok=True)

    src_f, src_t = os.path.join(pkg, "figures"), os.path.join(pkg, "tables")
    n = 0
    for i, stem in enumerate(MAIN_FIGURES, 1):
        for ext in ("pdf", "png"):
            if copy(os.path.join(src_f, f"{stem}.{ext}"), d_fig, name=f"Figure{i}.{ext}"):
                n += 1
    for i, stem in enumerate(SUPP_FIGURES, 1):
        for ext in ("pdf", "png"):
            if copy(os.path.join(src_f, f"{stem}.{ext}"), d_sup_f, name=f"FigureS{i}.{ext}"):
                n += 1
    for i, stem in enumerate(MAIN_TABLES, 1):
        if copy(os.path.join(src_t, f"{stem}.csv"), d_tab, name=f"Table{i}.csv"):
            n += 1
    for fn in sorted(os.listdir(src_t)) if os.path.isdir(src_t) else []:
        if fn.startswith("tableS"):
            i = fn.split("_")[0].replace("tableS", "")
            if copy(os.path.join(src_t, fn), d_sup_t, name=f"TableS{i}.csv"):
                n += 1

    # The additional files are the supplemental upload itself, named as the portal expects and in
    # citation order; the components that were compiled into Additional file 1 sit beneath them.
    sup_src = os.path.join(pkg, "supplemental")
    xls = os.path.join(sup_src, "excel_tables")
    rows = []
    for src_name, dest_name, what in ADDITIONAL_FILES:
        src = (os.path.join(sup_src, src_name) if src_name.endswith(".pdf")
               else os.path.join(xls, src_name))
        if copy(src, d_sup, name=dest_name):
            rows.append(f"| `{dest_name}` | {what} |")
            n += 1
        else:
            print(f"  ! additional file source missing: {src_name}")
    with open(os.path.join(d_sup, "README.md"), "w", encoding="utf-8") as fh:
        fh.write("# Additional files\n\nUploaded individually, in citation order. The manuscript "
                 "refers to the items inside them by their S-numbers, which are kept because "
                 "readers and reviewers use them; this table maps the two.\n\n"
                 "| File | Contents |\n|---|---|\n" + "\n".join(rows) + "\n\n"
                 "Additional file 1 is the same compiled supplement that a Cell Press submission "
                 "would call Document S1.\n")

    # The graphical abstract is a separate, required upload rather than a numbered figure.
    ga = os.path.join(pkg, "figures", "graphical_abstract.png")
    if copy(ga, d_ga):
        n += 1
    else:
        print("  ! graphical_abstract.png not found; BMC requires it at submission")

    for fn in MANUSCRIPT_DOCS:
        copy(os.path.join(pkg, fn), d_ms)

    # BMC accepts DOC/DOCX/RTF/TeX for the main text and will not take markdown or PDF, so the
    # submission manuscript is also emitted as Word, with the double spacing and the line and page
    # numbering the guidelines ask for. Optional dependency: a missing python-docx must not fail
    # the whole package.
    md = os.path.join(d_ms, MANUSCRIPT_DOCS[0])
    if os.path.exists(md):
        try:
            from deploy.make_manuscript_docx import convert
            out = convert(md, md.replace(".md", ".docx"))
            print(f"  + {os.path.basename(out)}")
        except ImportError:
            print("  ! python-docx unavailable; the .docx BMC requires was not generated")
    for fn in ARCHIVE_DOCS:
        copy(os.path.join(pkg, fn), d_arc)

    # The flat copies were staging for this layout; remove them so there is one home per file.
    for path in [os.path.join(pkg, x) for x in ("figures", "tables", "supplemental")]:
        shutil.rmtree(path, ignore_errors=True)
    for fn in MANUSCRIPT_DOCS + ARCHIVE_DOCS:
        flat = os.path.join(pkg, fn)
        if os.path.exists(flat):
            os.remove(flat)
    print(f"submission layout: {n} display items placed")
    return n


# Descriptions for the manifest. Exact relative paths win; otherwise the longest matching path
# prefix applies, so a directory rule covers everything under it without listing each file.
FILE_NOTES = {
    "README.md": "what this deposit is, how it is organised, and how to reproduce it",
    "MANIFEST.md": "this file",
    "CITATION.cff": "machine-readable citation metadata",
    "LICENSE-code.txt": "licence covering the code record (MIT)",
    "LICENSE-data.txt": "licence covering the data record (CC BY 4.0)",
    "benchmarks/diagnosis_arm/lirical_run/lirical_run_commands.sh":
        "the exact LIRICAL 2.4.1 invocation for all 23 cases, one line per case, as run on the VM",
    "benchmarks/diagnosis_arm/lirical_run/lirical_raw_outputs.tar.gz":
        "raw per-case LIRICAL output, unmodified, as produced by the commands above",
}
DIR_NOTES = {
    "benchmarks/variant_arm": "variant-arm benchmark output: metrics, per-variant scores, and the "
                              "published calibration fold assignment",
    "benchmarks/diagnosis_arm": "diagnosis-arm benchmark output, including the LIRICAL comparison "
                                "and its per-case rankings",
    "benchmarks": "benchmark harness output, as written by the harnesses themselves",
    "code": "the analysis code as committed, with the harnesses that produce every reported number",
    "configs/vcep_specs": "VCEP specifications the engine applies, one file per panel",
    "configs/cluster_definitions": "confusable-disease cluster definitions and their likelihood "
                                   "ratios, each entry sourced to a PMID",
    "configs": "engine configuration as run",
    "data_manifests": "third-party sources with URL, version and role; none is redistributed here",
    "docs": "the Phase R record, the claims map, and the supplemental methods and legends",
    "preregistration": "the pre-registered analysis plan and its gates",
}


def describe(rel):
    if rel in FILE_NOTES:
        return FILE_NOTES[rel]
    best = ""
    for prefix in DIR_NOTES:
        if rel.startswith(prefix + "/") and len(prefix) > len(best):
            best = prefix
    return DIR_NOTES.get(best, "")


def write_manifest(dep):
    rows = []
    for base, _dirs, files in os.walk(dep):
        for f in sorted(files):
            if f == "MANIFEST.md":
                continue
            p = os.path.join(base, f)
            rel = os.path.relpath(p, dep).replace("\\", "/")
            rows.append((rel, os.path.getsize(p), md5(p), describe(rel)))
    rows.sort()
    total = sum(r[1] for r in rows)
    undescribed = [r[0] for r in rows if not r[3]]
    with open(os.path.join(dep, "MANIFEST.md"), "w", encoding="utf-8") as fh:
        fh.write("# MANIFEST\n\nEvery file in this deposit with its size, md5 checksum, and what it "
                 "is.\n\nVerify an archive is intact by recomputing the checksums:\n\n"
                 "```\nfind . -type f ! -name MANIFEST.md -exec md5sum {} + | sort -k2\n```\n\n")
        fh.write(f"{len(rows)} files, {total / 1e6:.1f} MB total.\n\n")
        fh.write("| file | bytes | md5 | description |\n|---|---|---|---|\n")
        for rel, size, h, note in rows:
            fh.write(f"| `{rel}` | {size} | `{h}` | {note} |\n")
    if undescribed:
        print(f"  ! {len(undescribed)} file(s) have no manifest description, e.g. {undescribed[:3]}")
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
