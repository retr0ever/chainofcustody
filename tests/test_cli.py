import pytest
from click.testing import CliRunner

from chainofcustody.initial import GeneNotFoundError
from chainofcustody.cli import main


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture(autouse=True)
def mock_get_canonical_cds(mocker):
    return mocker.patch("chainofcustody.cli.get_canonical_cds")


def test_cli_prints_cds(runner, mock_get_canonical_cds):
    mock_get_canonical_cds.return_value = "ATGCCCAAA"

    result = runner.invoke(main, ["BRCA1"])

    assert result.exit_code == 0
    assert "ATGCCCAAA" in result.output
    mock_get_canonical_cds.assert_called_once_with("BRCA1")


def test_cli_gene_not_found_exits_with_error(runner, mock_get_canonical_cds):
    mock_get_canonical_cds.side_effect = GeneNotFoundError("Gene 'XYZ' not found")

    result = runner.invoke(main, ["XYZ"])

    assert result.exit_code == 1
    assert "XYZ" in result.output


def test_cli_missing_argument_exits_with_error(runner):
    result = runner.invoke(main, [])

    assert result.exit_code != 0
