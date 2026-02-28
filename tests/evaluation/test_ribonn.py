import pandas as pd
import pytest

from chainofcustody.sequence import mRNASequence
from chainofcustody.evaluation.ribonn import (
    _RIBONN_DIR,
    _aggregate,
    _score_with_dir,
    _te_status,
    score_ribonn,
)

_PARSED = mRNASequence(utr5="AAAAAA", cds="AUGCCCAAGUAA", utr3="CCCGGG")


def _mock_predictions(hela: float = 1.8, hepg2: float = 1.2, folds: int = 2) -> pd.DataFrame:
    """Build a fake multi-fold predictions DataFrame like RiboNN returns."""
    rows = []
    for fold in range(folds):
        rows.append({
            "tx_id": "query",
            "utr5_sequence": "AAAAAA",
            "cds_sequence": "AUGCCCAAGUAA",
            "utr3_sequence": "CCCGGG",
            "predicted_TE_HeLa": hela,
            "predicted_TE_HepG2": hepg2,
            "fold": fold,
        })
    return pd.DataFrame(rows)


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


# ── _aggregate (pure function) ───────────────────────────────────────────────

def test_aggregate_mean_te():
    df = _mock_predictions(hela=1.8, hepg2=1.2, folds=3)
    result = _aggregate(df)
    assert result["mean_te"] == pytest.approx(1.5, abs=1e-4)
    assert result["status"] == "GREEN"
    assert "HeLa" in result["per_tissue"]
    assert result["per_tissue"]["HeLa"] == pytest.approx(1.8)


def test_aggregate_no_available_key():
    result = _aggregate(_mock_predictions())
    assert "available" not in result


def test_aggregate_tissue_names_stripped():
    result = _aggregate(_mock_predictions())
    assert all(not k.startswith("predicted_") for k in result["per_tissue"])
    assert all(not k.startswith("TE_") for k in result["per_tissue"])


def test_aggregate_amber():
    df = _mock_predictions(hela=1.1, hepg2=1.3)
    assert _aggregate(df)["status"] == "AMBER"


def test_aggregate_red():
    df = _mock_predictions(hela=0.5, hepg2=0.8)
    assert _aggregate(df)["status"] == "RED"


# ── _score_with_dir (mocked prediction fn) ──────────────────────────────────

def _fake_runs_csv(tmp_path: "Path") -> pd.DataFrame:
    """Create a minimal runs.csv and return as DataFrame."""
    runs = pd.DataFrame({"run_id": ["abc123"], "params.test_fold": [0]})
    (tmp_path / "models" / "human").mkdir(parents=True)
    runs.to_csv(tmp_path / "models" / "human" / "runs.csv", index=False)
    return runs


def test_score_with_dir_success(tmp_path, mocker):
    _fake_runs_csv(tmp_path)

    from chainofcustody.evaluation.ribonn import _ensure_importable
    _ensure_importable()
    import src.predict as src_predict

    mocker.patch.object(
        src_predict,
        "predict_using_nested_cross_validation_models",
        side_effect=lambda *a, **kw: _mock_predictions(),
    )

    result = _score_with_dir(_PARSED, tmp_path)

    assert result["mean_te"] == pytest.approx(1.5, abs=1e-4)
    assert result["status"] == "GREEN"
    assert "available" not in result


def test_score_with_dir_propagates_exception(tmp_path, mocker):
    _fake_runs_csv(tmp_path)

    from chainofcustody.evaluation.ribonn import _ensure_importable
    _ensure_importable()
    import src.predict as src_predict

    mocker.patch.object(
        src_predict,
        "predict_using_nested_cross_validation_models",
        side_effect=RuntimeError("model weights missing"),
    )

    with pytest.raises(RuntimeError, match="model weights missing"):
        _score_with_dir(_PARSED, tmp_path)
