"""Shared Rich live-status display for the evaluation pipeline.

Any module can call :func:`update_status` to push a one-line description of
what it's doing right now, and :func:`update_best_score` to broadcast the
best overall fitness seen so far. The CLI wires both up via
:func:`set_status_callback` / :func:`set_best_score_callback`; when no display
is attached the calls are no-ops.

Usage (in problem.py, ribonn.py, etc.)::

    from chainofcustody.progress import update_status, update_best_score
    update_status("RiboNN  loading 50 models into GPU…")
    update_best_score(0.72)
"""

from __future__ import annotations

from typing import Callable

_status_callback: Callable[[str], None] | None = None
_best_score_callback: Callable[[float], None] | None = None


def set_status_callback(fn: Callable[[str], None] | None) -> None:
    """Register a callback that receives status strings.  Pass ``None`` to clear."""
    global _status_callback
    _status_callback = fn


def set_best_score_callback(fn: Callable[[float], None] | None) -> None:
    """Register a callback that receives the best overall fitness after each generation."""
    global _best_score_callback
    _best_score_callback = fn


def update_status(message: str) -> None:
    """Push a status message to whatever display is currently registered."""
    if _status_callback is not None:
        _status_callback(message)


def update_best_score(score: float) -> None:
    """Push the best overall fitness score to whatever display is registered."""
    if _best_score_callback is not None:
        _best_score_callback(score)
