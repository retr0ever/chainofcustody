import pandas as pd
import numpy as np
import pytest

from chainofcustody.sequence import mRNASequence
from chainofcustody.evaluation.ribonn import (
    _RIBONN_DIR,
    _MAX_UTR5_LEN,
    _MAX_CDS_UTR3_LEN,
    _null_result,
    _te_status,
    _encode_sequences_vectorized,
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


# ── Encoding ─────────────────────────────────────────────────────────────────

def test_encode_full_sequence():
    """The full transcript (5'UTR + CDS + UTR3) must be encoded."""
    utr5 = "AUGCCC"   # 6 nt
    cds  = "AUGCCCAAGUAA"  # 12 nt
    utr3 = "CCCCCCCCCC"   # 10 nt
    seq = mRNASequence(utr5=utr5, cds=cds, utr3=utr3)
    tensor, valid = _encode_sequences_vectorized([seq])

    assert valid[0] is True
    arr = tensor.numpy()  # (1, 5, _PADDED_LEN)

    # CDS starts at _MAX_UTR5_LEN; first few positions should be non-zero
    cds_start = _MAX_UTR5_LEN
    assert arr[0, :4, cds_start : cds_start + len(cds)].sum() > 0, (
        "CDS region should be encoded"
    )


def test_encode_codon_start_mask_is_set():
    """Channel 4 must mark the first nucleotide of every CDS codon."""
    utr5 = "AAAAAA"       # 6 nt
    cds  = "AUGCCCAAGUAA" # 12 nt → 4 codons → 4 codon-start positions
    seq = mRNASequence(utr5=utr5, cds=cds, utr3="")
    tensor, valid = _encode_sequences_vectorized([seq])

    assert valid[0] is True
    arr = tensor.numpy()
    cds_start = _MAX_UTR5_LEN
    cds_len = len(cds)
    expected_positions = list(range(cds_start, cds_start + cds_len - 3 + 1, 3))
    assert arr[0, 4, expected_positions].sum() == len(expected_positions)
    # No codon-start marks outside the CDS
    mask_sum = arr[0, 4, :].sum()
    assert mask_sum == len(expected_positions)


def test_encode_utr5_is_right_aligned():
    """The 5'UTR must be right-aligned so it ends at position _MAX_UTR5_LEN - 1."""
    utr5 = "AAUU"   # 4 nt
    seq = mRNASequence(utr5=utr5, cds="AUGAAGUAA", utr3="")
    tensor, valid = _encode_sequences_vectorized([seq])

    assert valid[0]
    arr = tensor.numpy()
    utr5_len = len(utr5)

    # The 4 nt should occupy columns [_MAX_UTR5_LEN - 4, _MAX_UTR5_LEN - 1]
    occupied = arr[0, :4, _MAX_UTR5_LEN - utr5_len : _MAX_UTR5_LEN]
    assert occupied.sum() == utr5_len, "Each of the 4 nt positions should have exactly one hot channel"

    # Positions before the 5'UTR should be zero
    assert np.all(arr[0, :4, : _MAX_UTR5_LEN - utr5_len] == 0.0)


def test_encode_utr5_too_long_is_invalid():
    """Sequences whose 5'UTR exceeds _MAX_UTR5_LEN must be marked invalid."""
    long_utr5 = "A" * (_MAX_UTR5_LEN + 1)
    seq = mRNASequence(utr5=long_utr5, cds="AUGAAGUAA", utr3="")
    tensor, valid = _encode_sequences_vectorized([seq])

    assert valid[0] is False
    assert np.all(tensor.numpy()[0] == 0.0)


def test_encode_cds_too_long_is_invalid():
    """Sequences whose CDS+UTR3 exceeds _MAX_CDS_UTR3_LEN must be marked invalid."""
    long_cds = "AUGCCC" * 2000 + "UAA"   # well over _MAX_CDS_UTR3_LEN
    seq = mRNASequence(utr5="AAAAAA", cds=long_cds, utr3="")
    tensor, valid = _encode_sequences_vectorized([seq])

    assert valid[0] is False
    assert np.all(tensor.numpy()[0] == 0.0)
