import pandas as pd
import pytest

from chainofcustody.sequence import mRNASequence
from chainofcustody.evaluation.ribonn import (
    _RIBONN_DIR,
    _null_result,
    _te_status,
    score_ribonn,
    score_ribonn_batch,
)

_PARSED = mRNASequence(utr5="AAAAAA", cds="AUGCCCAAGUAA", utr3="CCCGGG")


# ── Submodule path ───────────────────────────────────────────────────────────

def test_submodule_dir_points_to_vendor():
    assert _RIBONN_DIR.name == "RiboNN"
    assert (_RIBONN_DIR / "src" / "main.py").exists()


# ── Status thresholds ────────────────────────────────────────────────────────

def test_te_status_green():
    assert _te_status(2.0, 1.0) == "GREEN"   # target >= 1.5 and diff >= 0.5
    assert _te_status(1.5, 0.5) == "GREEN"


def test_te_status_amber():
    assert _te_status(1.0, 0.8) == "AMBER"   # target >= 1.0 and diff >= 0.0
    assert _te_status(1.2, 1.1) == "AMBER"


def test_te_status_red():
    assert _te_status(0.5, 0.8) == "RED"     # target < 1.0
    assert _te_status(0.0, 0.0) == "RED"


# ── _null_result ─────────────────────────────────────────────────────────────

def test_null_result_shape():
    r = _null_result()
    assert r["mean_te"] == 0.0
    assert r["target_te"] == 0.0
    assert r["mean_off_target_te"] == 0.0
    assert r["status"] == "RED"
    assert r["per_tissue"] is None
    assert r["target_cell_type"] == "megakaryocytes"


# ── score_ribonn (single-sequence, mocked predictor) ─────────────────────────

def _make_fake_result(target_te: float = 2.0, mean_off: float = 0.9) -> dict:
    return {
        "mean_te": round((target_te + mean_off) / 2, 4),
        "target_cell_type": "megakaryocytes",
        "target_te": target_te,
        "mean_off_target_te": mean_off,
        "per_tissue": {"HeLa": 0.9, "megakaryocytes": target_te},
        "status": _te_status(target_te, mean_off),
        "message": f"RiboNN: megakaryocytes TE = {target_te:.4f}, mean off-target = {mean_off:.4f}",
    }


def test_score_ribonn_delegates_to_predictor(mocker):
    """score_ribonn should call get_predictor().predict_batch and return its result."""
    fake_predictor = mocker.MagicMock()
    fake_predictor.predict_batch.return_value = [_make_fake_result(2.0, 0.9)]
    mocker.patch("chainofcustody.evaluation.ribonn.get_predictor", return_value=fake_predictor)

    result = score_ribonn(_PARSED)

    fake_predictor.predict_batch.assert_called_once_with([_PARSED], target_cell_type="megakaryocytes")
    assert result["target_te"] == pytest.approx(2.0, abs=1e-4)
    assert result["status"] == "GREEN"


def test_score_ribonn_batch_delegates_to_predictor(mocker):
    """score_ribonn_batch should call get_predictor().predict_batch."""
    seqs = [_PARSED, _PARSED]
    fake_results = [_make_fake_result(2.0, 0.9), _make_fake_result(0.5, 0.8)]
    fake_predictor = mocker.MagicMock()
    fake_predictor.predict_batch.return_value = fake_results
    mocker.patch("chainofcustody.evaluation.ribonn.get_predictor", return_value=fake_predictor)

    results = score_ribonn_batch(seqs)

    fake_predictor.predict_batch.assert_called_once_with(seqs, target_cell_type="megakaryocytes")
    assert len(results) == 2
    assert results[0]["status"] == "GREEN"
    assert results[1]["status"] == "RED"


def test_score_ribonn_result_keys(mocker):
    """Result dict must contain the expected keys."""
    fake_predictor = mocker.MagicMock()
    fake_predictor.predict_batch.return_value = [_make_fake_result()]
    mocker.patch("chainofcustody.evaluation.ribonn.get_predictor", return_value=fake_predictor)

    result = score_ribonn(_PARSED)
    assert {"mean_te", "target_cell_type", "target_te", "mean_off_target_te", "per_tissue", "status", "message"}.issubset(result.keys())
    assert all(not k.startswith("predicted_") for k in (result["per_tissue"] or {}))
    assert all(not k.startswith("TE_") for k in (result["per_tissue"] or {}))
