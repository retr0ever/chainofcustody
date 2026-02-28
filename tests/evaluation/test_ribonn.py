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
    assert _te_status(1.5) == "GREEN"
    assert _te_status(2.5) == "GREEN"


def test_te_status_amber():
    assert _te_status(1.0) == "AMBER"
    assert _te_status(1.49) == "AMBER"


def test_te_status_red():
    assert _te_status(0.5) == "RED"
    assert _te_status(0.0) == "RED"


# ── _null_result ─────────────────────────────────────────────────────────────

def test_null_result_shape():
    r = _null_result()
    assert r["mean_te"] == 0.0
    assert r["status"] == "RED"
    assert r["per_tissue"] is None


# ── score_ribonn (single-sequence, mocked predictor) ─────────────────────────

def _make_fake_result(mean_te: float = 1.5) -> dict:
    return {
        "mean_te": mean_te,
        "per_tissue": {"HeLa": 1.8, "HepG2": 1.2},
        "status": _te_status(mean_te),
        "message": f"RiboNN predicted mean TE = {mean_te:.4f}",
    }


def test_score_ribonn_delegates_to_predictor(mocker):
    """score_ribonn should call get_predictor().predict_batch and return its result."""
    fake_predictor = mocker.MagicMock()
    fake_predictor.predict_batch.return_value = [_make_fake_result(1.5)]
    mocker.patch("chainofcustody.evaluation.ribonn.get_predictor", return_value=fake_predictor)

    result = score_ribonn(_PARSED)

    fake_predictor.predict_batch.assert_called_once_with([_PARSED])
    assert result["mean_te"] == pytest.approx(1.5, abs=1e-4)
    assert result["status"] == "GREEN"


def test_score_ribonn_batch_delegates_to_predictor(mocker):
    """score_ribonn_batch should call get_predictor().predict_batch."""
    seqs = [_PARSED, _PARSED]
    fake_results = [_make_fake_result(1.5), _make_fake_result(0.8)]
    fake_predictor = mocker.MagicMock()
    fake_predictor.predict_batch.return_value = fake_results
    mocker.patch("chainofcustody.evaluation.ribonn.get_predictor", return_value=fake_predictor)

    results = score_ribonn_batch(seqs)

    fake_predictor.predict_batch.assert_called_once_with(seqs)
    assert len(results) == 2
    assert results[0]["status"] == "GREEN"
    assert results[1]["status"] == "RED"


def test_score_ribonn_result_keys(mocker):
    """Result dict must contain the expected keys."""
    fake_predictor = mocker.MagicMock()
    fake_predictor.predict_batch.return_value = [_make_fake_result()]
    mocker.patch("chainofcustody.evaluation.ribonn.get_predictor", return_value=fake_predictor)

    result = score_ribonn(_PARSED)
    assert set(result.keys()) == {"mean_te", "per_tissue", "status", "message"}
    assert all(not k.startswith("predicted_") for k in (result["per_tissue"] or {}))
    assert all(not k.startswith("TE_") for k in (result["per_tissue"] or {}))
