"""Shared Rich live-status display for the evaluation pipeline.

Any module can call :func:`update_status` to push a one-line description of
what it's doing right now. The CLI wires the display up via
:func:`set_live_status`; when no display is attached the calls are no-ops.

Usage (in problem.py, ribonn.py, etc.)::

    from chainofcustody.progress import update_status
    update_status("RiboNN  loading 50 models into GPU…")
"""

from __future__ import annotations

from typing import Callable

_status_callback: Callable[[str], None] | None = None


def set_status_callback(fn: Callable[[str], None] | None) -> None:
    """Register a callback that receives status strings.  Pass ``None`` to clear."""
    global _status_callback
    _status_callback = fn


def update_status(message: str) -> None:
    """Push a status message to whatever display is currently registered."""
    if _status_callback is not None:
        _status_callback(message)
