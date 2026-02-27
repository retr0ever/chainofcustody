import pytest
from click.testing import CliRunner

from chainofcustody.cli import main


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture(autouse=True)
def mock_get_canonical_cds(mocker):
    # Patch at the `initial` package level — where the lazy import resolves to
    return mocker.patch("chainofcustody.initial.get_canonical_cds")


def test_fetch_prints_cds(runner, mock_get_canonical_cds):
    mock_get_canonical_cds.return_value = "ATGCCCAAA"

    result = runner.invoke(main, ["fetch", "BRCA1"])

    assert result.exit_code == 0
    assert "ATGCCCAAA" in result.output
    mock_get_canonical_cds.assert_called_once_with("BRCA1")


def test_fetch_gene_not_found_exits_with_error(runner, mock_get_canonical_cds):
    from chainofcustody.initial import GeneNotFoundError
    mock_get_canonical_cds.side_effect = GeneNotFoundError("Gene 'XYZ' not found")

    result = runner.invoke(main, ["fetch", "XYZ"])

    assert result.exit_code == 1
    assert "XYZ" in result.output


def test_fetch_missing_argument_exits_with_error(runner):
    result = runner.invoke(main, ["fetch"])

    assert result.exit_code != 0


def test_no_subcommand_shows_help(runner):
    result = runner.invoke(main, ["--help"])

    assert result.exit_code == 0
    assert "fetch" in result.output
    assert "optimize" in result.output


# ── optimize ──────────────────────────────────────────────────────────────────

@pytest.fixture
def mock_optimize_run(mocker):
    import numpy as np
    mock = mocker.patch("chainofcustody.optimization.run")
    mock.return_value = (
        np.array([[0, 1, 2, 3]] * 3),
        np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]]),
    )
    return mock


def test_optimize_summary_output(runner, mock_optimize_run):
    result = runner.invoke(main, ["optimize", "--seq-len", "4", "--pop-size", "10", "--n-gen", "2"])

    assert result.exit_code == 0
    assert "Pareto front" in result.output
    mock_optimize_run.assert_called_once_with(seq_len=4, pop_size=10, n_gen=2, mutation_rate=0.01, seed=None)


def test_optimize_json_output(runner, mock_optimize_run):
    result = runner.invoke(main, ["optimize", "--seq-len", "4", "--output", "json"])

    assert result.exit_code == 0
    assert "sequence" in result.output
    assert "objectives" in result.output


def test_optimize_custom_mutation_rate(runner, mock_optimize_run):
    runner.invoke(main, ["optimize", "--mutation-rate", "0.05"])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["mutation_rate"] == 0.05


def test_optimize_seed_is_passed(runner, mock_optimize_run):
    runner.invoke(main, ["optimize", "--seed", "42"])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["seed"] == 42
