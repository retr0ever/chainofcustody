"""Metric 4: Translation efficiency prediction via RiboNN (Sanofi).

RiboNN is a deep-CNN tool that predicts per-tissue translation efficiency from
mRNA sequences. It is bundled as a git submodule at ``vendor/RiboNN`` and its
Python dependencies are installed in the project venv.

Pretrained weights are downloaded automatically on first run.

Reference: Karollus et al., Nature Biotechnology (2024).
Repo: https://github.com/Sanofi-Public/RiboNN
"""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from chainofcustody.sequence import mRNASequence

# Path to the RiboNN submodule, relative to this file:
#   chainofcustody/evaluation/ribonn.py  → parents[2] = repo root
_RIBONN_DIR = Path(__file__).parents[2] / "vendor" / "RiboNN"

# RiboNN hard-codes this input path relative to its own root (see vendor/RiboNN/src/main.py)
_INPUT_FILE = "data/prediction_input1.txt"
_OUTPUT_FILE = "results/human/prediction_output.txt"


def _te_status(mean_te: float) -> str:
    if mean_te >= 1.5:
        return "GREEN"
    if mean_te >= 1.0:
        return "AMBER"
    return "RED"


def score_ribonn(parsed: mRNASequence) -> dict:
    """Predict translation efficiency using RiboNN.

    Returns a dict with keys: ``mean_te``, ``per_tissue``, ``status``, ``message``.

    Raises:
        RuntimeError: if the RiboNN subprocess fails or produces no output.
    """
    return _run_ribonn(parsed, _RIBONN_DIR)


def _run_ribonn(parsed: mRNASequence, ribonn_dir: Path) -> dict:
    """Write input, invoke RiboNN subprocess, parse output."""
    target = ribonn_dir / _INPUT_FILE
    target.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        backup = Path(tmpdir) / "backup.txt"

        # Back up any pre-existing input file
        if target.exists():
            shutil.copy2(target, backup)

        target.write_text(
            "tx_id\tutr5_sequence\tcds_sequence\tutr3_sequence\n"
            f"query\t{parsed.utr5}\t{parsed.cds}\t{parsed.utr3}\n"
        )

        try:
            result = subprocess.run(
                [sys.executable, "-m", "src.main", "--predict", "human"],
                cwd=str(ribonn_dir),
                capture_output=True,
                text=True,
                timeout=300,
            )
        finally:
            if backup.exists():
                shutil.copy2(backup, target)
            elif target.exists():
                target.unlink()

    if result.returncode != 0:
        raise RuntimeError(
            f"RiboNN exited with code {result.returncode}: {result.stderr[:300]}"
        )

    return _parse_output(ribonn_dir)


def _parse_output(ribonn_dir: Path) -> dict:
    """Read and parse RiboNN's prediction_output.txt."""
    output_file = ribonn_dir / _OUTPUT_FILE

    if not output_file.exists():
        raise RuntimeError(f"RiboNN output file not found: {output_file}")

    with output_file.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    if not rows:
        raise RuntimeError("RiboNN produced empty output")

    row = rows[0]
    mean_te = float(row["mean_predicted_TE"])

    per_tissue = {}
    for key, val in row.items():
        if key.startswith("predicted_") and key != "mean_predicted_TE":
            try:
                per_tissue[key.removeprefix("predicted_")] = float(val)
            except (ValueError, TypeError):
                continue

    return {
        "mean_te": round(mean_te, 4),
        "per_tissue": per_tissue or None,
        "status": _te_status(mean_te),
        "message": f"RiboNN predicted mean TE = {mean_te:.4f}",
    }
