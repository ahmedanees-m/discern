# Security policy

## Reporting a vulnerability
Please report suspected security issues privately to the maintainer
(`ahmedanees-m` on GitHub) rather than opening a public issue. You will get an
acknowledgement, and fixes are prioritized over new features.

## Never commit secrets
The GitHub token and the NVIDIA/Nemotron API key live **outside** this repository
and are matched by `.gitignore` patterns (`*token*.txt`, `*API_key*.txt`,
`**/nvidia_api_key*`, `.env`, ...). If you ever add a credential file, confirm
`git status` does not list it before committing.

## Credentials at runtime
All credentials are read from environment variables; none are hard-coded:

| Variable | Use |
|---|---|
| `NVIDIA_API_KEY` | cloud Nemotron (NVIDIA NIM) |
| `VM_HOST`, `VM_USER`, `VM_PORT` | VM SSH target |
| `VM_KEY` *or* `VM_PASSWORD` | VM auth (key preferred) |
| `GITHUB_TOKEN` | release automation only |

Copy `deploy/.env.example` to `deploy/.env` (git-ignored) and fill in real values
locally. Never commit real values.

## VM access
Key-based SSH is recommended: generate a keypair, add the public key to the VM
`~/.ssh/authorized_keys`, and connect via the SSH agent. Password auth is
supported by `deploy/remote.py` (`VM_PASSWORD`) for convenience but should be
rotated and replaced with a key. The compute host sits behind a private network;
reach hosted services through an SSH tunnel or the authenticated proxy, never by
exposing them publicly.

## Patient data
This is decision support, not a medical device. **No real patient data** is stored
in the repository or in any public or hosted artifact. Public demos run
on synthetic and public data only. No protected health information is written to
logs or the audit ledger in hosted contexts.
