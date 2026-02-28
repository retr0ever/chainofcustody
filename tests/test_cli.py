import pytest
from click.testing import CliRunner

from chainofcustody.cli import main

_GENE = "TP53"
_CDS = "AUGCCCAAGUAA"  # minimal valid CDS: AUG + codon + stop
_UTR3 = "GAGTAGUCCC"
_UTR5_MIN = 4
_UTR5_MAX = 10
_OFF_TARGET_CELL_TYPE = "Dendritic_cell"


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def mock_get_cds(mocker):
    """Mock Ensembl CDS lookup so tests don't need network access."""
    mock = mocker.patch("chainofcustody.cli.get_canonical_cds")
    mock.return_value = _CDS.replace("U", "T")  # returns DNA like the real function
    return mock


@pytest.fixture
def mock_generate_utr3(mocker):
    """Mock three_prime 3'UTR generation so tests don't need the expression DB."""
    mock = mocker.patch("chainofcustody.cli.generate_utr3")
    mock.return_value = _UTR3
    return mock


@pytest.fixture
def mock_optimize_run(mocker):
    import numpy as np
    from chainofcustody.optimization.problem import METRIC_NAMES, N_OBJECTIVES
    mock = mocker.patch("chainofcustody.cli.run")
    mock_history = [
        {"generation": g, "sequence": "ACGU" + _CDS + _UTR3, **{m: 0.8 for m in METRIC_NAMES}, "overall": 0.8}
        for g in range(1, 4)
    ]
    # Column 0 = length (4), columns 1..10 = nucleotides
    X_row = np.array([_UTR5_MIN] + [0, 1, 2, 3] + [0] * (_UTR5_MAX - _UTR5_MIN))
    mock.return_value = (
        np.array([X_row] * 3),
        np.array([[0.3] * N_OBJECTIVES] * 3),
        mock_history,
    )
    return mock


@pytest.fixture
def mock_scoring(mocker):
    """Mock the evaluation pipeline so optimize tests don't need ViennaRNA etc."""
    mock_report = {
        "sequence_info": {"total_length": 4, "utr5_length": 0, "cds_length": 3, "utr3_length": 0, "num_codons": 1},
        "structure_scores": {"utr5_accessibility": {"mfe": None, "status": "GREY"}, "global_mfe": {"mfe": -1.0, "mfe_per_nt": -0.25}},
        "manufacturing_scores": {"total_violations": 0, "overall_pass": True, "gc_windows": {"pass": True, "violations": []}, "homopolymers": {"pass": True, "violations": []}, "restriction_sites": {"pass": True, "violations": []}},
        "stability_scores": {"gc3": 0.5, "mfe_per_nt": -0.3, "au_rich_elements": 0, "stability_score": 0.7, "status": "GREEN"},
        "summary": {"utr5_accessibility": "GREY", "manufacturability": "GREEN", "stability": "GREEN"},
    }
    mock_fitness = {
        "scores": {m: {"value": 0.8, "weight": 0.1, "weighted": 0.08, "status": "GREEN"} for m in ["utr5_accessibility", "manufacturability", "stability"]},
        "overall": 0.8,
        "suggestions": [],
    }
    mocker.patch("chainofcustody.cli.score_parsed", return_value=mock_report)
    mocker.patch("chainofcustody.cli.compute_fitness", return_value=mock_fitness)


def test_help(runner):
    result = runner.invoke(main, ["--help"])

    assert result.exit_code == 0
    assert "NSGA3" in result.output


def test_summary_output(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring):
    result = runner.invoke(main, [
        "--gene", _GENE,
        "--utr5-min", str(_UTR5_MIN), "--utr5-max", str(_UTR5_MAX),
        "--pop-size", "10", "--n-gen", "2",
    ])

    assert result.exit_code == 0, result.output
    assert "Candidate Ranking" in result.output or "Pareto front" in result.output
    _, kwargs = mock_optimize_run.call_args
    assert {k: kwargs[k] for k in ("utr5_min", "utr5_max", "pop_size", "n_gen", "mutation_rate", "seed", "n_workers")} == {
        "utr5_min": _UTR5_MIN, "utr5_max": _UTR5_MAX,
        "pop_size": 10, "n_gen": 2, "mutation_rate": 0.01, "seed": None, "n_workers": None,
    }


def test_json_output(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring):
    result = runner.invoke(main, ["--gene", _GENE, "--output", "json"])

    assert result.exit_code == 0
    assert "fitness" in result.output


def test_custom_mutation_rate(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring):
    runner.invoke(main, ["--gene", _GENE, "--mutation-rate", "0.05"])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["mutation_rate"] == 0.05


def test_seed_is_passed(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring):
    runner.invoke(main, ["--gene", _GENE, "--seed", "42"])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["seed"] == 42


def test_workers_is_passed(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring):
    runner.invoke(main, ["--gene", _GENE, "--workers", "4"])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["n_workers"] == 4


def test_csv_output(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring, tmp_path):
    csv_file = tmp_path / "history.csv"
    result = runner.invoke(main, ["--gene", _GENE, "--csv", str(csv_file)])

    assert result.exit_code == 0
    assert csv_file.exists()
    assert "History written to" in result.output

    import csv as csv_mod
    rows = list(csv_mod.DictReader(csv_file.open()))
    assert len(rows) == 3
    assert set(rows[0].keys()) >= {"generation", "sequence", "overall"}


def test_gene_not_found(runner, mocker, mock_generate_utr3):
    from chainofcustody.cds import GeneNotFoundError
    mocker.patch("chainofcustody.cli.get_canonical_cds", side_effect=GeneNotFoundError("BADGENE not found"))
    result = runner.invoke(main, ["--gene", "BADGENE"])

    assert result.exit_code == 1
    assert "Error" in result.output


def test_utr5_min_gt_max_rejected(runner, mock_get_cds):
    result = runner.invoke(main, ["--gene", _GENE, "--utr5-min", "50", "--utr5-max", "10"])

    assert result.exit_code == 1
    assert "Error" in result.output


def test_off_target_cell_type_is_passed(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring):
    runner.invoke(main, ["--gene", _GENE, "--off-target-cell-type", "Hepatocyte_derived"])

    mock_generate_utr3.assert_called_once_with("Hepatocyte_derived")


def test_off_target_cell_type_default(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring):
    runner.invoke(main, ["--gene", _GENE])

    mock_generate_utr3.assert_called_once_with(_OFF_TARGET_CELL_TYPE)


def test_off_target_cell_type_not_found(runner, mock_get_cds, mocker):
    mocker.patch("chainofcustody.cli.generate_utr3", side_effect=ValueError("Unknown cell type 'BADTYPE'"))
    result = runner.invoke(main, ["--gene", _GENE, "--off-target-cell-type", "BADTYPE"])

    assert result.exit_code == 1
    assert "Error" in result.output
