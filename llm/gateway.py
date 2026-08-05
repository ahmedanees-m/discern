"""LLM gateway.

Deployment decision for this project: the reasoning/soft-task LLM is the
**cloud-hosted NVIDIA Nemotron** endpoint (NVIDIA NIM, OpenAI-compatible) rather
than a self-hosted Ollama model. This removes the 16 GB on-VM GPU constraint for
LLM inference and frees the A4000 for other work.

Guardrails (unchanged from the plan):
  * The LLM is used ONLY for soft tasks (notes->HPO, literature proposals,
    phrasing a deterministic plan). It NEVER produces the classification or points.
  * Every call is logged to the audit ledger with model + version.
  * RAG answers must carry citations; the rule engine assigns points.

Secrets: the API key is read from the ``NVIDIA_API_KEY`` environment variable.
It is never hardcoded, committed, or baked into an image.
"""
from __future__ import annotations

import os

# OpenAI-compatible client; import is deferred so the module loads without the
# dependency present (laptop dev / CI on tiny fixtures).
try:  # pragma: no cover - thin import guard
    from openai import OpenAI
except Exception:  # pragma: no cover
    OpenAI = None  # type: ignore

NVIDIA_BASE_URL = os.environ.get("NVIDIA_BASE_URL", "https://integrate.api.nvidia.com/v1")

# Default cloud model: Nemotron Super 49B (Llama-3.3-based), served free from the
# build.nvidia.com catalog. Overridable via NEMOTRON_MODEL. (The older
# llama-3.1-nemotron-70b now 404s, hence pinning to the 3.3 Super 49B build.)
NEMOTRON_MODEL = os.environ.get("NEMOTRON_MODEL", "nvidia/llama-3.3-nemotron-super-49b-v1")

# Task -> model routing. Both tasks default to cloud Nemotron; an operator may
# point "fast" at a local OpenAI-compatible endpoint via OLLAMA_URL if desired.
ROUTES = {
    "fast": os.environ.get("DISCERN_FAST_MODEL", NEMOTRON_MODEL),
    "reason": os.environ.get("DISCERN_REASON_MODEL", NEMOTRON_MODEL),
}


def _client(task: str) -> OpenAI:
    if OpenAI is None:  # pragma: no cover
        raise RuntimeError("openai package not installed in this environment")
    # "fast" may optionally use a local Ollama endpoint; "reason" always cloud.
    if task == "fast" and os.environ.get("OLLAMA_URL"):
        return OpenAI(base_url=os.environ["OLLAMA_URL"], api_key="ollama")
    key = os.environ.get("NVIDIA_API_KEY")
    if not key:
        raise RuntimeError("NVIDIA_API_KEY not set in the environment")
    return OpenAI(base_url=NVIDIA_BASE_URL, api_key=key)


def complete(task: str, prompt: str, *, audit=None, **kw) -> str:
    """Run a soft-task completion. ``task`` is one of ROUTES.

    Pass an ``AuditLedger`` as ``audit`` to log model + version for the call.
    """
    if task not in ROUTES:
        raise ValueError(f"unknown task route {task!r}; expected one of {list(ROUTES)}")
    model = ROUTES[task]
    client = _client(task)
    resp = client.chat.completions.create(
        model=model,
        messages=[{"role": "user", "content": prompt}],
        **kw,
    )
    out = resp.choices[0].message.content or ""
    if audit is not None:
        audit.append("llm_call", {"task": task, "model": model, "chars_in": len(prompt),
                                  "chars_out": len(out)})
    return out
