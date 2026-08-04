"""Run lint and tests against exactly what continuous integration sees.

A working tree is not the repository. `app/`, the root `Dockerfile`, `tests/test_app.py` and
`manuscript/` are gitignored on purpose, so a local `ruff check .` and a local `pytest` both cover
files CI never receives - and, less obviously, ruff's own import classification changes when a
package directory is present, which once made a file lint clean locally and fail in CI.

This materialises the tracked-file set into a temporary directory and runs the CI commands there.
Anything it reports is what CI will report.

Run:  python -m deploy.check_like_ci
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main() -> int:
    tracked = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True)
    if tracked.returncode != 0:
        print("not a git repository", file=sys.stderr)
        return 2
    files = [f for f in tracked.stdout.splitlines() if f.strip()]

    tmp = tempfile.mkdtemp(prefix="discern-ci-")
    try:
        for rel in files:
            src = os.path.join(ROOT, rel)
            if not os.path.isfile(src):
                continue
            dst = os.path.join(tmp, rel.replace("/", os.sep))
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(src, dst)
        print(f"materialised {len(files)} tracked files into a clean tree")

        rc = 0
        for label, cmd in (("lint (ruff)", [sys.executable, "-m", "ruff", "check", "."]),
                           ("tests", [sys.executable, "-m", "pytest", "-q"])):
            r = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True)
            tail = (r.stdout or r.stderr).strip().splitlines()
            print(f"\n== {label}: {'pass' if r.returncode == 0 else 'FAIL'} ==")
            print("\n".join(tail[-15:]))
            rc = rc or r.returncode
        print("\nthis is what CI will report" if rc == 0 else "\nCI would fail on the above")
        return rc
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
