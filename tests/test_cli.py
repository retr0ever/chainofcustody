import pytest
from click.testing import CliRunner

from chainofcustody.cli import main


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def mock_optimize_run(mocker):
    import numpy as np
    from chainofcustody.optimization.problem import N_OBJECTIVES
    mock = mocker.patch("chainofcustody.cli.run")
    mock.return_value = (
        np.array([[0, 1, 2, 3]] * 3),
        np.array([[0.3] * N_OBJECTIVES] * 3),
    )
    return mock


@pytest.fixture
def mock_scoring(mocker):
    """Mock the evaluation pipeline so optimize tests don't need ViennaRNA etc."""
    mock_report = {
        "sequence_info": {"total_length": 4, "utr5_length": 0, "cds_length": 3, "utr3_length": 0, "num_codons": 1},
        "codon_scores": {"cai": 0.8, "gc_content": {"overall": 50, "cds": 50}, "liver_selectivity": 0.0},
        "mirna_scores": {"detargeting": {"miR-122-5p": {"utr3_sites": 0, "total_sites": 0}}},
        "structure_scores": {"utr5_accessibility": {"mfe": None, "status": "GREY"}, "global_mfe": {"mfe": -1.0, "mfe_per_nt": -0.25}},
        "manufacturing_scores": {"total_violations": 0, "overall_pass": True, "gc_windows": {"pass": True, "violations": []}, "homopolymers": {"pass": True, "violations": []}, "restriction_sites": {"pass": True, "violations": []}},
        "stability_scores": {"gc3": 0.5, "mfe_per_nt": -0.3, "au_rich_elements": 0, "stability_score": 0.7, "status": "GREEN"},
        "summary": {"codon_quality": "GREEN", "gc_content": "GREEN", "mir122_detargeting": "RED", "utr5_accessibility": "GREY", "manufacturability": "GREEN", "stability": "GREEN"},
    }
    mock_fitness = {
        "scores": {m: {"value": 0.8, "weight": 0.1, "weighted": 0.08, "status": "GREEN"} for m in ["codon_quality", "gc_content", "mir122_detargeting", "utr5_accessibility", "manufacturability", "stability"]},
        "overall": 0.8,
        "suggestions": [],
    }
    mocker.patch("chainofcustody.cli.score_sequence", return_value=mock_report)
    mocker.patch("chainofcustody.cli.compute_fitness", return_value=mock_fitness)


def test_help(runner):
    result = runner.invoke(main, ["--help"])

    assert result.exit_code == 0
    assert "NSGA3" in result.output


def test_summary_output(runner, mock_optimize_run, mock_scoring):
    result = runner.invoke(main, ["--seq-len", "4", "--pop-size", "10", "--n-gen", "2"])

    assert result.exit_code == 0
    assert "Candidate Ranking" in result.output or "Pareto front" in result.output
    mock_optimize_run.assert_called_once_with(seq_len=4, pop_size=10, n_gen=2, mutation_rate=0.01, seed=None, n_workers=None)


def test_json_output(runner, mock_optimize_run, mock_scoring):
    result = runner.invoke(main, ["--seq-len", "4", "--output", "json"])

    assert result.exit_code == 0
    assert "fitness" in result.output


def test_custom_mutation_rate(runner, mock_optimize_run, mock_scoring):
    runner.invoke(main, ["--mutation-rate", "0.05"])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["mutation_rate"] == 0.05


def test_seed_is_passed(runner, mock_optimize_run, mock_scoring):
    runner.invoke(main, ["--seed", "42"])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["seed"] == 42


def test_workers_is_passed(runner, mock_optimize_run, mock_scoring):
    runner.invoke(main, ["--workers", "4"])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["n_workers"] == 4


def test_csv_output(runner, mock_optimize_run, mock_scoring, tmp_path):
    csv_file = tmp_path / "results.csv"
    result = runner.invoke(main, ["--seq-len", "4", "--csv", str(csv_file)])

    assert result.exit_code == 0
    assert csv_file.exists()

    import csv as csv_mod
    rows = list(csv_mod.DictReader(csv_file.open()))
    assert len(rows) == 3
    assert set(rows[0].keys()) >= {"rank", "label", "sequence", "overall"}
