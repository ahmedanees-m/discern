"""VM orchestration over SSH/SFTP.

Division of labor: the laptop authors code and drives the VM; the VM runs all
heavy work inside Docker. This helper opens an SSH session and provides SFTP
sync (we use SFTP for data transfer, not rclone).

Secrets policy: host / user / auth are read from the environment, never from
code or git. Supported auth, in order of preference:
  1. SSH agent / default key   (set VM_HOST, VM_USER)
  2. explicit private key file  (VM_KEY=/path/to/id_ed25519)
  3. password                   (VM_PASSWORD=...)  -- convenient but key auth is preferred

Example env (PowerShell):
    $env:VM_HOST="your.vm.hostname"; $env:VM_USER="your_vm_user"; $env:VM_PASSWORD="..."
"""
from __future__ import annotations

import os
from contextlib import contextmanager

import paramiko


def vm() -> paramiko.SSHClient:
    """Open an SSH connection to the VM using env-configured credentials."""
    host = os.environ["VM_HOST"]
    user = os.environ["VM_USER"]
    port = int(os.environ.get("VM_PORT", "22"))
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.load_system_host_keys()
    kwargs: dict = {"hostname": host, "port": port, "username": user, "timeout": 30}
    if os.environ.get("VM_KEY"):
        kwargs["key_filename"] = os.environ["VM_KEY"]
    elif os.environ.get("VM_PASSWORD"):
        kwargs["password"] = os.environ["VM_PASSWORD"]
    # else: rely on the SSH agent / default key
    c.connect(**kwargs)
    return c


@contextmanager
def session():
    c = vm()
    try:
        yield c
    finally:
        c.close()


def run(cmd: str) -> tuple[str, str]:
    """Run a single command on the VM; return (stdout, stderr)."""
    with session() as c:
        _, out, err = c.exec_command(cmd)
        return out.read().decode("utf-8", "replace"), err.read().decode("utf-8", "replace")


def put(local: str, remote: str) -> None:
    """Upload one file laptop -> VM via SFTP."""
    with session() as c:
        sftp = c.open_sftp()
        try:
            sftp.put(local, remote)
        finally:
            sftp.close()


def get(remote: str, local: str) -> None:
    """Download one file VM -> laptop via SFTP."""
    with session() as c:
        sftp = c.open_sftp()
        try:
            sftp.get(remote, local)
        finally:
            sftp.close()


# Example (do not auto-run):
#   run("cd ~/discern && docker compose -f deploy/compose.vm.yml up -d --build")
