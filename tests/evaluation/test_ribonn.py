import subprocess

import pytest

from chainofcustody.sequence import mRNASequence
from chainofcustody.evaluation.ribonn import score_ribonn, _te_status, _parse_output


_PARSED = mRNASequence(utr5="AAAAAA", cds="AUGCCCAAGUAA", utr3="CCCGGG")


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


# ── Graceful fallback ────────────────────────────────────────────────────────

def test_unavailable_when_no_env_var(monkeypatch):
    monkeypatch.delenv("RIBONN_DIR", raising=False)
    # Reset cached dir
    import chainofcustody.evaluation.ribonn as mod
    mod._RIBONN_DIR = None

    result = score_ribonn(_PARSED)

    assert result["available"] is False
    assert result["mean_te"] is None
    assert result["status"] == "GREY"


def test_unavailable_when_dir_missing(monkeypatch, tmp_path):
    monkeypatch.setenv("RIBONN_DIR", str(tmp_path / "nonexistent"))
    import chainofcustody.evaluation.ribonn as mod
    mod._RIBONN_DIR = None

    result = score_ribonn(_PARSED)

    assert result["available"] is False
    assert result["status"] == "GREY"


def test_unavailable_when_no_main_py(monkeypatch, tmp_path):
    # Dir exists but no src/main.py
    monkeypatch.setenv("RIBONN_DIR", str(tmp_path))
    import chainofcustody.evaluation.ribonn as mod
    mod._RIBONN_DIR = None

    result = score_ribonn(_PARSED)

    assert result["available"] is False
    assert result["status"] == "GREY"


# ── Subprocess mocking ───────────────────────────────────────────────────────

def test_successful_prediction(monkeypatch, tmp_path, mocker):
    # Set up a fake RiboNN dir
    src_dir = tmp_path / "src"
    src_dir.mkdir()
    (src_dir / "main.py").write_text("")
    (tmp_path / "data").mkdir()

    monkeypatch.setenv("RIBONN_DIR", str(tmp_path))
    import chainofcustody.evaluation.ribonn as mod
    mod._RIBONN_DIR = None

    # Mock subprocess.run to succeed
    mocker.patch("chainofcustody.evaluation.ribonn.subprocess.run", return_value=subprocess.CompletedProcess(
        args=[], returncode=0, stdout="", stderr="",
    ))

    # Create fake output
    results_dir = tmp_path / "results" / "human"
    results_dir.mkdir(parents=True)
    (results_dir / "prediction_output.txt").write_text(
        "tx_id\tutr5_sequence\tcds_sequence\tutr3_sequence\tpredicted_HeLa\tpredicted_HepG2\tmean_predicted_TE\n"
        "query\tAAAAAA\tAUGCCCAAGUAA\tCCCGGG\t1.8\t1.2\t1.5\n"
    )

    result = score_ribonn(_PARSED)

    assert result["available"] is True
    assert result["mean_te"] == 1.5
    assert result["status"] == "GREEN"
    assert result["per_tissue"]["HeLa"] == 1.8
    assert result["per_tissue"]["HepG2"] == 1.2


def test_subprocess_failure(monkeypatch, tmp_path, mocker):
    src_dir = tmp_path / "src"
    src_dir.mkdir()
    (src_dir / "main.py").write_text("")
    (tmp_path / "data").mkdir()

    monkeypatch.setenv("RIBONN_DIR", str(tmp_path))
    import chainofcustody.evaluation.ribonn as mod
    mod._RIBONN_DIR = None

    mocker.patch("chainofcustody.evaluation.ribonn.subprocess.run", return_value=subprocess.CompletedProcess(
        args=[], returncode=1, stdout="", stderr="Model not found",
    ))

    result = score_ribonn(_PARSED)

    assert result["available"] is False
    assert result["status"] == "GREY"


def test_subprocess_timeout(monkeypatch, tmp_path, mocker):
    src_dir = tmp_path / "src"
    src_dir.mkdir()
    (src_dir / "main.py").write_text("")
    (tmp_path / "data").mkdir()

    monkeypatch.setenv("RIBONN_DIR", str(tmp_path))
    import chainofcustody.evaluation.ribonn as mod
    mod._RIBONN_DIR = None

    mocker.patch(
        "chainofcustody.evaluation.ribonn.subprocess.run",
        side_effect=subprocess.TimeoutExpired(cmd="python3", timeout=300),
    )

    result = score_ribonn(_PARSED)

    assert result["available"] is False
    assert result["status"] == "GREY"
