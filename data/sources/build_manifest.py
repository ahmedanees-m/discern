"""Build / verify data/manifest.json (plan Step 0.4).

Walks ``data/raw`` and ``data/processed``, records each file's size and md5, and
writes the result into ``data/manifest.json`` so the reference store is
provenance-logged and reproducible. Per-source download scripts live alongside
this file and register their URL/version/date into the same manifest.

Run: ``python -m data.sources.build_manifest``
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
MANIFEST = DATA / "manifest.json"


def _md5(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def scan() -> dict:
    files: dict[str, dict] = {}
    for sub in ("raw", "processed"):
        base = DATA / sub
        if not base.exists():
            continue
        for p in sorted(base.rglob("*")):
            if p.is_file():
                rel = p.relative_to(DATA).as_posix()
                files[rel] = {"size": p.stat().st_size, "md5": _md5(p)}
    return files


def main() -> None:
    manifest = {"manifest_version": "0.1.0", "sources": {}}
    if MANIFEST.exists():
        manifest = json.loads(MANIFEST.read_text())
    manifest["files"] = scan()
    MANIFEST.write_text(json.dumps(manifest, indent=2))
    print(f"manifest: {len(manifest['files'])} files indexed -> {MANIFEST}")


if __name__ == "__main__":
    main()
