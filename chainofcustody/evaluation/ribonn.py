"""Metric 4: Translation efficiency prediction via RiboNN (Sanofi).

RiboNN is a deep-CNN tool that predicts per-tissue translation efficiency from
mRNA sequences. It is bundled as a git submodule at ``vendor/RiboNN`` and its
Python dependencies are installed in the project venv.

Pretrained weights are downloaded automatically on first run.

Reference: Karollus et al., Nature Biotechnology (2024).
Repo: https://github.com/Sanofi-Public/RiboNN
"""

from __future__ import annotations

import contextlib
import sys
import tempfile
from pathlib import Path

import pandas as pd

from chainofcustody.sequence import mRNASequence

# Path to the RiboNN submodule, relative to this file:
#   chainofcustody/evaluation/ribonn.py  →  parents[2] = repo root
_RIBONN_DIR = Path(__file__).parents[2] / "vendor" / "RiboNN"

_SPECIES = "human"
_TOP_K = 5


def _ensure_importable() -> None:
    """Add vendor/RiboNN to sys.path so src.* modules can be imported."""
    ribonn_str = str(_RIBONN_DIR)
    if ribonn_str not in sys.path:
        sys.path.insert(0, ribonn_str)


def _te_status(mean_te: float) -> str:
    if mean_te >= 1.5:
        return "GREEN"
    if mean_te >= 1.0:
        return "AMBER"
    return "RED"


def score_ribonn(parsed: mRNASequence) -> dict:
    """Predict translation efficiency using RiboNN.

    Returns a dict with keys: ``mean_te``, ``per_tissue``, ``status``,
    ``message``.

    Raises:
        RuntimeError: if RiboNN prediction fails.
    """
    return _score_with_dir(parsed, _RIBONN_DIR)


def _score_with_dir(parsed: mRNASequence, ribonn_dir: Path) -> dict:
    """Run RiboNN from *ribonn_dir*. Separated for testability."""
    _ensure_importable()
    from src.predict import predict_using_nested_cross_validation_models  # noqa: PLC0415

    run_df = pd.read_csv(ribonn_dir / "models" / _SPECIES / "runs.csv")

    # RiboNNDataModule still needs a file path — write to system tmp, not into
    # the submodule, so the submodule working tree stays clean.
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=True) as fh:
        fh.write("tx_id\tutr5_sequence\tcds_sequence\tutr3_sequence\n")
        fh.write(f"query\t{parsed.utr5}\t{parsed.cds}\t{parsed.utr3}\n")
        fh.flush()
        input_path = fh.name

        # contextlib.chdir is required: model weight paths in predict.py are
        # relative strings (e.g. "models/human/{run_id}/state_dict.pth").
        with contextlib.chdir(ribonn_dir):
            predictions = predict_using_nested_cross_validation_models(
                input_path,
                _SPECIES,
                run_df,
                top_k_models_to_use=_TOP_K,
                batch_size=32,
                num_workers=0,  # 0 avoids fork/spawn issues on macOS
            )

    return _aggregate(predictions)


def _aggregate(predictions: pd.DataFrame) -> dict:
    """Average across folds and return the scored result dict."""
    predicted_cols = [c for c in predictions.columns if c.startswith("predicted_")]

    row = (
        predictions
        .groupby(
            ["tx_id", "utr5_sequence", "cds_sequence", "utr3_sequence"],
            as_index=False,
        )[predicted_cols]
        .mean()
        .iloc[0]
    )

    mean_te = float(row[predicted_cols].mean())
    per_tissue = {
        c.removeprefix("predicted_TE_"): round(float(row[c]), 4)
        for c in predicted_cols
    }

    return {
        "mean_te": round(mean_te, 4),
        "per_tissue": per_tissue or None,
        "status": _te_status(mean_te),
        "message": f"RiboNN predicted mean TE = {mean_te:.4f}",
    }
