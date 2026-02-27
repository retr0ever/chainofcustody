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
    from chainofcustody.optimization.problem import N_OBJECTIVES
    mock = mocker.patch("chainofcustody.optimization.run")
    mock.return_value = (
        np.array([[0, 1, 2, 3]] * 3),
        np.array([[0.3] * N_OBJECTIVES] * 3),
    )
    return mock


@pytest.fixture
def mock_scoring(mocker):
    """Mock the evaluation pipeline so optimize tests don't need ViennaRNA etc."""
    mock_score = mocker.patch("chainofcustody.cli.score_sequence", create=True)
    mock_fitness = mocker.patch("chainofcustody.cli.compute_fitness", create=True)

    # Patch at the module level where the optimize command imports them
    mock_score = mocker.patch("chainofcustody.evaluation.report.score_sequence")
    mock_score.return_value = {
        "sequence_info": {"total_length": 4, "utr5_length": 0, "cds_length": 3, "utr3_length": 0, "num_codons": 1},
        "codon_scores": {"cai": 0.8, "gc_content": {"overall": 50, "cds": 50}, "liver_selectivity": 0.0},
        "mirna_scores": {"detargeting": {"miR-122-5p": {"utr3_sites": 0, "total_sites": 0}}},
        "structure_scores": {"utr5_accessibility": {"mfe": None, "status": "GREY"}, "global_mfe": {"mfe": -1.0, "mfe_per_nt": -0.25}},
        "manufacturing_scores": {"total_violations": 0, "overall_pass": True, "gc_windows": {"pass": True, "violations": []}, "homopolymers": {"pass": True, "violations": []}, "restriction_sites": {"pass": True, "violations": []}},
        "stability_scores": {"gc3": 0.5, "mfe_per_nt": -0.3, "au_rich_elements": 0, "stability_score": 0.7, "status": "GREEN"},
        "summary": {"codon_quality": "GREEN", "gc_content": "GREEN", "mir122_detargeting": "RED", "utr5_accessibility": "GREY", "manufacturability": "GREEN", "stability": "GREEN"},
    }
    return mock_score


def test_optimize_summary_output(runner, mock_optimize_run, mock_scoring):
    result = runner.invoke(main, ["optimize", "--seq-len", "4", "--pop-size", "10", "--n-gen", "2"])

    assert result.exit_code == 0
    assert "Candidate Ranking" in result.output or "Pareto front" in result.output
    mock_optimize_run.assert_called_once_with(seq_len=4, pop_size=10, n_gen=2, mutation_rate=0.01, seed=None)


def test_optimize_json_output(runner, mock_optimize_run, mock_scoring):
    result = runner.invoke(main, ["optimize", "--seq-len", "4", "--output", "json"])

    assert result.exit_code == 0
    assert "fitness" in result.output


def test_optimize_custom_mutation_rate(runner, mock_optimize_run, mock_scoring):
    runner.invoke(main, ["optimize", "--mutation-rate", "0.05"])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["mutation_rate"] == 0.05


def test_optimize_seed_is_passed(runner, mock_optimize_run, mock_scoring):
    runner.invoke(main, ["optimize", "--seed", "42"])

    _, kwargs = mock_optimize_run.call_args
    assert kwargs["seed"] == 42
