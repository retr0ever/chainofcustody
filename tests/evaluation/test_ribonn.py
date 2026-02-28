import subprocess

import pytest

from chainofcustody.sequence import mRNASequence
from chainofcustody.evaluation.ribonn import score_ribonn, _te_status, _parse_output, _RIBONN_DIR


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


# ── Subprocess mocking ───────────────────────────────────────────────────────

def _fake_ribonn_dir(tmp_path):
    """Create a minimal fake RiboNN directory layout."""
    (tmp_path / "src").mkdir()
    (tmp_path / "src" / "main.py").write_text("")
    (tmp_path / "data").mkdir()
    return tmp_path


def test_successful_prediction(tmp_path, mocker):
    ribonn_dir = _fake_ribonn_dir(tmp_path)

    mocker.patch("chainofcustody.evaluation.ribonn.subprocess.run", return_value=subprocess.CompletedProcess(
        args=[], returncode=0, stdout="", stderr="",
    ))

    results_dir = ribonn_dir / "results" / "human"
    results_dir.mkdir(parents=True)
    (results_dir / "prediction_output.txt").write_text(
        "tx_id\tutr5_sequence\tcds_sequence\tutr3_sequence\tpredicted_HeLa\tpredicted_HepG2\tmean_predicted_TE\n"
        "query\tAAAAAA\tAUGCCCAAGUAA\tCCCGGG\t1.8\t1.2\t1.5\n"
    )

    result = score_ribonn.__wrapped__(_PARSED) if hasattr(score_ribonn, "__wrapped__") else \
        __import__("chainofcustody.evaluation.ribonn", fromlist=["_run_ribonn"])._run_ribonn(_PARSED, ribonn_dir)

    assert result["mean_te"] == 1.5
    assert result["status"] == "GREEN"
    assert result["per_tissue"]["HeLa"] == 1.8
    assert result["per_tissue"]["HepG2"] == 1.2
    assert "available" not in result


def test_subprocess_failure_raises(tmp_path, mocker):
    ribonn_dir = _fake_ribonn_dir(tmp_path)

    mocker.patch("chainofcustody.evaluation.ribonn.subprocess.run", return_value=subprocess.CompletedProcess(
        args=[], returncode=1, stdout="", stderr="Model not found",
    ))

    import chainofcustody.evaluation.ribonn as mod
    with pytest.raises(RuntimeError, match="RiboNN exited with code 1"):
        mod._run_ribonn(_PARSED, ribonn_dir)


def test_subprocess_timeout_propagates(tmp_path, mocker):
    ribonn_dir = _fake_ribonn_dir(tmp_path)

    mocker.patch(
        "chainofcustody.evaluation.ribonn.subprocess.run",
        side_effect=subprocess.TimeoutExpired(cmd="python3", timeout=300),
    )

    import chainofcustody.evaluation.ribonn as mod
    with pytest.raises(subprocess.TimeoutExpired):
        mod._run_ribonn(_PARSED, ribonn_dir)


def test_missing_output_raises(tmp_path, mocker):
    ribonn_dir = _fake_ribonn_dir(tmp_path)

    mocker.patch("chainofcustody.evaluation.ribonn.subprocess.run", return_value=subprocess.CompletedProcess(
        args=[], returncode=0, stdout="", stderr="",
    ))

    import chainofcustody.evaluation.ribonn as mod
    with pytest.raises(RuntimeError, match="output file not found"):
        mod._run_ribonn(_PARSED, ribonn_dir)


# ── Input file backup/restore ────────────────────────────────────────────────

def test_existing_input_is_restored_after_run(tmp_path, mocker):
    ribonn_dir = _fake_ribonn_dir(tmp_path)
    input_path = ribonn_dir / "data" / "prediction_input1.txt"
    original_content = "original data\n"
    input_path.write_text(original_content)

    mocker.patch("chainofcustody.evaluation.ribonn.subprocess.run", return_value=subprocess.CompletedProcess(
        args=[], returncode=1, stdout="", stderr="fail",
    ))

    import chainofcustody.evaluation.ribonn as mod
    with pytest.raises(RuntimeError):
        mod._run_ribonn(_PARSED, ribonn_dir)

    assert input_path.read_text() == original_content
