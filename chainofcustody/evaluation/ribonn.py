"""Metric 4: Translation efficiency prediction via RiboNN (Sanofi).

RiboNN is an external deep-CNN tool that predicts per-tissue translation
efficiency from mRNA sequences.  It is optional — when unavailable the metric
degrades gracefully to GREY / 0.5.

Configure the path to a local RiboNN clone via the ``RIBONN_DIR`` environment
variable.  The directory must contain ``src/main.py`` and will auto-download
pretrained weights on first run.

Reference: Karollus et al., Nature Biotechnology (2024).
Repo: https://github.com/Sanofi-Public/RiboNN
"""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from pathlib import Path

from chainofcustody.sequence import mRNASequence

logger = logging.getLogger(__name__)

_RIBONN_DIR: Path | None = None


def _get_ribonn_dir() -> Path | None:
    """Resolve the RiboNN installation directory (cached)."""
    global _RIBONN_DIR
    if _RIBONN_DIR is not None:
        return _RIBONN_DIR if _RIBONN_DIR.exists() else None

    env = os.environ.get("RIBONN_DIR")
    if env:
        p = Path(env)
        if p.is_dir() and (p / "src" / "main.py").exists():
            _RIBONN_DIR = p
            return p
    return None


def _unavailable(reason: str = "RiboNN not available") -> dict:
    return {
        "available": False,
        "mean_te": None,
        "per_tissue": None,
        "status": "GREY",
        "message": reason,
    }


def _te_status(mean_te: float) -> str:
    if mean_te >= 1.5:
        return "GREEN"
    if mean_te >= 1.0:
        return "AMBER"
    return "RED"


def score_ribonn(parsed: mRNASequence) -> dict:
    """Predict translation efficiency using RiboNN.

    Returns a dict with keys: ``available``, ``mean_te``, ``per_tissue``,
    ``status``, ``message``.
    """
    ribonn_dir = _get_ribonn_dir()
    if ribonn_dir is None:
        return _unavailable()

    try:
        return _run_ribonn(parsed, ribonn_dir)
    except Exception as exc:
        logger.warning("RiboNN prediction failed: %s", exc)
        return _unavailable(f"RiboNN error: {exc}")


def _run_ribonn(parsed: mRNASequence, ribonn_dir: Path) -> dict:
    """Write input, invoke RiboNN subprocess, parse output."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Write input in RiboNN's tab-separated format (Option A: split regions)
        input_file = tmp / "prediction_input.txt"
        input_file.write_text(
            "tx_id\tutr5_sequence\tcds_sequence\tutr3_sequence\n"
            f"query\t{parsed.utr5}\t{parsed.cds}\t{parsed.utr3}\n"
        )

        # Symlink data dir so RiboNN finds the input where it expects
        data_link = ribonn_dir / "data" / "prediction_input.txt"
        data_link_existed = data_link.exists()
        if not data_link_existed:
            data_link.parent.mkdir(parents=True, exist_ok=True)

        # Copy input into RiboNN's expected location
        import shutil
        backup = None
        if data_link.exists():
            backup = tmp / "backup_input.txt"
            shutil.copy2(data_link, backup)
        shutil.copy2(input_file, data_link)

        try:
            result = subprocess.run(
                ["python3", "-m", "src.main", "--predict", "human"],
                cwd=str(ribonn_dir),
                capture_output=True,
                text=True,
                timeout=300,
            )
        finally:
            # Restore original input if we overwrote it
            if backup and backup.exists():
                shutil.copy2(backup, data_link)
            elif not data_link_existed and data_link.exists():
                data_link.unlink()

        if result.returncode != 0:
            return _unavailable(f"RiboNN exited with code {result.returncode}: {result.stderr[:200]}")

        return _parse_output(ribonn_dir)


def _parse_output(ribonn_dir: Path) -> dict:
    """Read and parse RiboNN's prediction_output.txt."""
    output_file = ribonn_dir / "results" / "human" / "prediction_output.txt"
    if not output_file.exists():
        return _unavailable("RiboNN output file not found")

    import csv
    with output_file.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)

    if not rows:
        return _unavailable("RiboNN produced empty output")

    row = rows[0]

    # Extract mean TE
    mean_te = float(row.get("mean_predicted_TE", 0))

    # Extract per-tissue predictions (columns starting with "predicted_")
    per_tissue = {}
    for key, val in row.items():
        if key.startswith("predicted_") and key != "mean_predicted_TE":
            tissue = key.removeprefix("predicted_")
            try:
                per_tissue[tissue] = float(val)
            except (ValueError, TypeError):
                continue

    return {
        "available": True,
        "mean_te": round(mean_te, 4),
        "per_tissue": per_tissue if per_tissue else None,
        "status": _te_status(mean_te),
        "message": f"RiboNN predicted mean TE = {mean_te:.4f}",
    }
