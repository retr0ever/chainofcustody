import pytest
from click.testing import CliRunner

from chainofcustody.cli import main

_GENE = "TP53"
_CDS = "AUGCCCAAGUAA"  # minimal valid CDS: AUG + codon + stop
_UTR3 = "GAGTAGUCCC"
_UTR5_MIN = 4
_UTR5_MAX = 10
_TARGET = "Fibroblast"
_RIBONN_TARGET = "fibroblast"


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
    """Mock the evaluation pipeline so optimize tests don't need DRfold2/RiboNN."""
    mock_report = {
        "sequence_info": {"total_length": 4, "full_length": 127, "cap5_length": 3, "utr5_length": 0, "cds_length": 3, "utr3_length": 0, "num_codons": 1},
        "structure_scores": {"utr5_accessibility": {"mfe": None, "mfe_per_nt": None, "status": "GREY"}, "global_mfe": {"mfe": -1.0, "mfe_per_nt": -0.25}},
        "manufacturing_scores": {"total_violations": 0, "utr5_violations": 0, "overall_pass": True, "gc_windows": {"pass": True, "violations": []}, "homopolymers": {"pass": True, "violations": []}, "restriction_sites": {"pass": True, "violations": []}, "uorfs": {"pass": True, "count": 0, "positions": [], "violations": []}},
        "stability_scores": {"gc3": 0.5, "mfe_per_nt": -0.3, "stability_score": 0.7, "status": "GREEN"},
        "ribonn_scores": {
            "mean_te": 1.8,
            "target_cell_type": _RIBONN_TARGET,
            "target_te": 2.2,
            "mean_off_target_te": 0.9,
            "per_tissue": None,
            "status": "GREEN",
            "message": f"RiboNN: {_RIBONN_TARGET} TE = 2.2000, mean off-target = 0.9000",
        },
        "summary": {"utr5_accessibility": "GREY", "manufacturability": "GREEN", "stability": "GREEN", "specificity": "GREY"},
    }
    mock_fitness = {
        "scores": {m: {"value": 0.8, "weight": 0.1, "weighted": 0.08, "status": "GREEN"} for m in ["utr5_accessibility", "manufacturability", "stability", "specificity"]},
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
    assert {k: kwargs[k] for k in ("utr5_min", "utr5_max", "pop_size", "n_gen", "mutation_rate", "seed", "n_workers", "initial_length", "max_length_delta", "seed_from_data", "gradient_seed_steps")} == {
        "utr5_min": _UTR5_MIN, "utr5_max": _UTR5_MAX,
        "pop_size": 10, "n_gen": 2, "mutation_rate": 0.05, "seed": None, "n_workers": None,
        "initial_length": 200, "max_length_delta": 50,
        "seed_from_data": True, "gradient_seed_steps": 0,
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


def test_target_default(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring):
    """Default --target uses Fibroblast; 3'UTR generated for seed-map name, RiboNN gets mapped name."""
    result = runner.invoke(main, ["--gene", _GENE])

    assert result.exit_code == 0, result.output
    mock_generate_utr3.assert_called_once_with(_TARGET)
    _, kwargs = mock_optimize_run.call_args
    assert kwargs["target_cell_type"] == _RIBONN_TARGET


def test_target_custom(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring):
    """Explicit --target is forwarded correctly to both 3'UTR generation and optimizer."""
    result = runner.invoke(main, ["--gene", _GENE, "--target", "Fibroblast"])

    assert result.exit_code == 0, result.output
    mock_generate_utr3.assert_called_once_with("Fibroblast")
    _, kwargs = mock_optimize_run.call_args
    assert kwargs["target_cell_type"] == "fibroblast"


def test_target_invalid(runner, mock_get_cds, mock_generate_utr3):
    """An unrecognised --target exits with an error before touching any IO."""
    result = runner.invoke(main, ["--gene", _GENE, "--target", "NotARealCellType"])

    assert result.exit_code == 1
    assert "Error" in result.output


def test_target_error_propagates(runner, mock_get_cds, mocker):
    """A ValueError from generate_utr3 is shown as an error."""
    mocker.patch("chainofcustody.cli.generate_utr3", side_effect=ValueError("No miRNAs found"))
    result = runner.invoke(main, ["--gene", _GENE, "--target", _TARGET])

    assert result.exit_code == 1
    assert "Error" in result.output


def test_no_off_target_option(runner):
    """The old --off-target-cell-type flag must no longer exist."""
    result = runner.invoke(main, ["--help"])

    assert "off-target-cell-type" not in result.output
    assert "off_target_cell_type" not in result.output


def test_seed_from_data_default_is_true(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring):
    """--seed-from-data is enabled by default."""
    runner.invoke(main, ["--gene", _GENE])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["seed_from_data"] is True


def test_no_seed_from_data_flag(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring):
    """--no-seed-from-data disables MOESM3 warm-start seeding."""
    runner.invoke(main, ["--gene", _GENE, "--no-seed-from-data"])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["seed_from_data"] is False


def test_gradient_seed_steps_default_is_zero(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring):
    """--gradient-seed-steps defaults to 0 (disabled)."""
    runner.invoke(main, ["--gene", _GENE])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["gradient_seed_steps"] == 0


def test_gradient_seed_steps_passed(runner, mock_get_cds, mock_generate_utr3, mock_optimize_run, mock_scoring):
    """--gradient-seed-steps value is forwarded to run()."""
    runner.invoke(main, ["--gene", _GENE, "--gradient-seed-steps", "100"])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["gradient_seed_steps"] == 100
