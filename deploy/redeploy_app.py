"""Rebuild and redeploy the DISCERN web app on the VM.

The app is deliberately not in the public repository (`app/` and the root `Dockerfile` are
gitignored), so the deployment ships the working tree over SFTP rather than cloning from GitHub.
Everything server-side runs in Docker; nothing is installed on the VM host.

The sequence is written to be safe to re-run and safe to abandon half way:

  1. verify the local tree first - a broken build is never uploaded
  2. tag the running image as a rollback point before anything is replaced
  3. upload a content-addressed bundle and build the new image beside the old one
  4. smoke-test the new image on a scratch port, still not touching the live container
  5. swap only once the new image has answered correctly
  6. verify the live endpoint, and roll back automatically if it does not come up

Run:  python -m deploy.redeploy_app [--dry-run] [--keep N]

Environment: VM_HOST, VM_USER, VM_PASSWORD (never read from code or git).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tarfile
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

APP_NAME = "discern-app"
IMAGE = "discern-app"
LIVE_PORT = 8500
CONTAINER_PORT = 7860
SCRATCH_PORT = 8599                    # used only for the pre-swap smoke test
REMOTE_ROOT = "discern-app-src"

# Mirrors .dockerignore: what the image build would discard anyway, kept off the wire.
EXCLUDE_DIRS = {".git", ".github", "__pycache__", ".pytest_cache", ".ruff_cache", "data",
                "docs", "figures", "sim", ".venv", "venv", "_local", "manuscript", "node_modules"}
EXCLUDE_SUFFIX = (".pyc", ".pyo", ".vcf", ".gz", ".parquet", ".coverage")
KEEP_ANYWAY = ("app/static/kb.json",)   # small and required, despite a broad suffix rule


def sh(cmd, cwd=ROOT):
    return subprocess.run(cmd, shell=isinstance(cmd, str), cwd=cwd, capture_output=True,
                          text=True)


def _client():
    import paramiko
    host = os.environ.get("VM_HOST", "10.30.158.35")
    user = os.environ.get("VM_USER", "anees_22phd0670")
    pw = os.environ.get("VM_PASSWORD")
    if not pw:
        sys.exit("VM_PASSWORD is not set; refusing to guess credentials.")
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(host, username=user, password=pw, timeout=30, banner_timeout=30, auth_timeout=30)
    return c


def remote(c, cmd, timeout=900, check=True):
    _in, out, err = c.exec_command(cmd, timeout=timeout)
    _in.close()
    o = out.read().decode("utf-8", "replace")
    e = err.read().decode("utf-8", "replace")
    rc = out.channel.recv_exit_status()
    if check and rc != 0:
        raise RuntimeError(f"remote command failed ({rc}): {cmd}\n{o}\n{e}")
    return o, e, rc


# ---------------------------------------------------------------------------------------------
def verify_local():
    """Never upload a tree whose own tests do not pass."""
    print("1. verifying the local tree")
    lint = sh([sys.executable, "-m", "ruff", "check", "."])
    if lint.returncode != 0:
        sys.exit("lint failed locally; not deploying\n" + lint.stdout)
    tests = sh([sys.executable, "-m", "pytest", "-q", "tests/test_app.py"])
    if tests.returncode != 0:
        sys.exit("app tests failed locally; not deploying\n" + tests.stdout[-2000:])
    print("   lint clean, app tests pass")

    sys.path.insert(0, ROOT)
    from fastapi.testclient import TestClient

    from app.server import app
    cli = TestClient(app)
    checks = [
        ("F8 must lead haemophilia A", {"gene": "F8", "features": {"prolonged_aptt": True,
                                                                   "delayed_bleeding": True}},
         lambda r: r["leading_id"] == "hemophilia_a"),
        ("F9 must lead haemophilia B", {"gene": "F9", "features": {"prolonged_aptt": True,
                                                                   "delayed_bleeding": True}},
         lambda r: r["leading_id"] == "hemophilia_b"),
        ("GP1BA + DDAVP must hard-stop for type 2B",
         {"gene": "GP1BA", "planned_tx": "ddavp",
          "features": {"ripa_low_dose_enhanced": True, "ripa_mixing_platelet_origin": True}},
         lambda r: r["leading_id"] == "ptvwd"
         and any("HARD STOP" in f["message"] and "2B" in f["message"] for f in r["safety"])),
    ]
    for label, payload, ok in checks:
        r = cli.post("/api/analyze", json=payload).json()
        if not ok(r):
            sys.exit(f"local behaviour check failed: {label}\n{json.dumps(r)[:400]}")
        print(f"   ok  {label}")
    return checks


def build_bundle():
    print("2. packaging the working tree")
    fd, path = tempfile.mkstemp(suffix=".tar.gz", prefix="discern-app-")
    os.close(fd)
    n = 0
    with tarfile.open(path, "w:gz") as tar:
        for base, dirs, files in os.walk(ROOT):
            dirs[:] = [d for d in dirs if d not in EXCLUDE_DIRS]
            for f in files:
                p = os.path.join(base, f)
                rel = os.path.relpath(p, ROOT).replace(os.sep, "/")
                if rel not in KEEP_ANYWAY:
                    if f.endswith(EXCLUDE_SUFFIX) or rel.startswith("bench/data/"):
                        continue
                tar.add(p, arcname=rel)
                n += 1
    print(f"   {n} files, {os.path.getsize(path) / 1e6:.1f} MB")
    return path


def upload(c, local_path, remote_path):
    """Chunked put: sftp.put short-reads on this link, and a truncated bundle builds a broken image."""
    sftp = c.open_sftp()
    try:
        size = os.path.getsize(local_path)
        with open(local_path, "rb") as rf, sftp.open(remote_path, "wb") as wf:
            wf.set_pipelined(True)
            while True:
                block = rf.read(32768)
                if not block:
                    break
                wf.write(block)
        got = sftp.stat(remote_path).st_size
        if got != size:
            raise RuntimeError(f"upload truncated: {got} != {size}")
    finally:
        sftp.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true", help="verify and package, do not touch the VM")
    ap.add_argument("--keep", type=int, default=3, help="rollback images to retain")
    args = ap.parse_args()

    verify_local()
    bundle = build_bundle()
    sha = sh(["git", "rev-parse", "--short", "HEAD"]).stdout.strip() or "nogit"
    dirty = bool(sh(["git", "status", "--porcelain"]).stdout.strip())
    tag = f"{sha}{'-dirty' if dirty else ''}-{time.strftime('%Y%m%d%H%M%S')}"
    print(f"   image tag: {IMAGE}:{tag}")
    if args.dry_run:
        print(f"\ndry run: bundle left at {bundle}")
        return

    c = _client()
    try:
        print("3. uploading and building beside the running container")
        remote(c, f"rm -rf ~/{REMOTE_ROOT} && mkdir -p ~/{REMOTE_ROOT}")
        upload(c, bundle, f"{REMOTE_ROOT}.tar.gz")
        remote(c, f"tar xzf ~/{REMOTE_ROOT}.tar.gz -C ~/{REMOTE_ROOT} && rm ~/{REMOTE_ROOT}.tar.gz")

        prev = remote(c, f"docker inspect {APP_NAME} --format '{{{{.Image}}}}' 2>/dev/null || true",
                      check=False)[0].strip()
        if prev:
            remote(c, f"docker tag {prev} {IMAGE}:rollback-{time.strftime('%Y%m%d%H%M%S')}",
                   check=False)
            print("   tagged the running image as a rollback point")

        out, err, _ = remote(c, f"cd ~/{REMOTE_ROOT} && docker build -q -t {IMAGE}:{tag} .",
                             timeout=1800)
        print(f"   built {IMAGE}:{tag} ({out.strip()[:20]})")

        print("4. smoke-testing the new image on a scratch port")
        remote(c, f"docker rm -f {APP_NAME}-probe 2>/dev/null || true", check=False)
        remote(c, f"docker run -d --name {APP_NAME}-probe -p {SCRATCH_PORT}:{CONTAINER_PORT} "
                  f"{IMAGE}:{tag}")
        ok = False
        try:
            for _ in range(20):
                time.sleep(3)
                h, _e, rc = remote(c, f"curl -s -m 5 http://127.0.0.1:{SCRATCH_PORT}/healthz",
                                   check=False)
                if rc == 0 and "ok" in h:
                    ok = True
                    break
            if not ok:
                raise RuntimeError("new image never became healthy on the scratch port")
            probe, _e, _rc = remote(
                c, f"curl -s -m 20 -X POST http://127.0.0.1:{SCRATCH_PORT}/api/analyze "
                   "-H 'Content-Type: application/json' "
                   """-d '{"gene":"F8","features":{"prolonged_aptt":true,"delayed_bleeding":true}}'""")
            lead = json.loads(probe).get("leading_id")
            if lead != "hemophilia_a":
                raise RuntimeError(f"new image still mis-routes F8 (leading={lead}); not swapping")
            print("   new image healthy and routes F8 correctly")
        finally:
            remote(c, f"docker rm -f {APP_NAME}-probe 2>/dev/null || true", check=False)

        print("5. swapping the live container")
        remote(c, f"docker rm -f {APP_NAME} 2>/dev/null || true", check=False)
        remote(c, f"docker tag {IMAGE}:{tag} {IMAGE}:latest")
        remote(c, f"docker run -d --name {APP_NAME} --restart unless-stopped "
                  f"-p {LIVE_PORT}:{CONTAINER_PORT} {IMAGE}:{tag}")

        print("6. verifying the live endpoint")
        live = False
        for _ in range(20):
            time.sleep(3)
            h, _e, rc = remote(c, f"curl -s -m 5 http://127.0.0.1:{LIVE_PORT}/healthz", check=False)
            if rc == 0 and "ok" in h:
                live = True
                break
        if not live:
            print("   ! live endpoint did not come up - rolling back")
            rb = remote(c, f"docker images {IMAGE} --format '{{{{.Tag}}}}' | grep '^rollback-' | head -1",
                        check=False)[0].strip()
            if rb:
                remote(c, f"docker rm -f {APP_NAME} 2>/dev/null || true", check=False)
                remote(c, f"docker run -d --name {APP_NAME} --restart unless-stopped "
                          f"-p {LIVE_PORT}:{CONTAINER_PORT} {IMAGE}:{rb}")
                sys.exit(f"rolled back to {IMAGE}:{rb}")
            sys.exit("deployment failed and no rollback image was found")
        print("   live and healthy")

        keep, _e, _rc = remote(
            c, f"docker images {IMAGE} --format '{{{{.Tag}}}}' | grep '^rollback-' | tail -n +{args.keep + 1}",
            check=False)
        for old in [t for t in keep.split() if t]:
            remote(c, f"docker rmi {IMAGE}:{old} 2>/dev/null || true", check=False)
        remote(c, "docker image prune -f >/dev/null 2>&1 || true", check=False)
        print(f"\ndeployed {IMAGE}:{tag} on port {LIVE_PORT}")
    finally:
        c.close()
        os.unlink(bundle)


if __name__ == "__main__":
    main()
