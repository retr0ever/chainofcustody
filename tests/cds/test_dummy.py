import pytest
from unittest.mock import MagicMock, patch

from chainofcustody.cds import GeneNotFoundError, get_canonical_cds


@pytest.fixture(autouse=True)
def mock_mygene(mocker):
    mg_mock = mocker.patch("chainofcustody.cds.lookup._mg")
    mg_mock.query.return_value = {
        "hits": [{"ensembl": {"gene": "ENSG00000012048"}}]
    }
    return mg_mock


@pytest.fixture
def mock_requests(mocker):
    return mocker.patch("chainofcustody.cds.lookup.requests.get")


def _make_response(json_data: dict, status_code: int = 200) -> MagicMock:
    response = MagicMock()
    response.status_code = status_code
    response.json.return_value = json_data
    response.raise_for_status = MagicMock()
    return response


def test_get_canonical_cds_returns_sequence(mock_requests):
    mock_requests.side_effect = [
        _make_response({"canonical_transcript": "ENST00000357654"}),
        _make_response({"seq": "ATGCCCAAA"}),
    ]

    result = get_canonical_cds("BRCA1")

    assert result == "ATGCCCAAA"
    assert mock_requests.call_count == 2


def test_get_canonical_cds_gene_not_found(mock_mygene):
    mock_mygene.query.return_value = {"hits": []}

    with pytest.raises(GeneNotFoundError):
        get_canonical_cds("UNKNOWN_GENE_XYZ")


def test_get_canonical_cds_no_canonical_transcript(mock_requests):
    mock_requests.return_value = _make_response({"id": "ENSG00000012048"})

    with pytest.raises(GeneNotFoundError, match="No canonical transcript"):
        get_canonical_cds("BRCA1")


def test_get_canonical_cds_multiple_ensembl_entries(mock_mygene, mock_requests):
    """mygene can return a list when a symbol maps to multiple Ensembl entries."""
    mock_mygene.query.return_value = {
        "hits": [{"ensembl": [{"gene": "ENSG00000012048"}, {"gene": "ENSG00000099999"}]}]
    }
    mock_requests.side_effect = [
        _make_response({"canonical_transcript": "ENST00000357654"}),
        _make_response({"seq": "ATGAAATTT"}),
    ]

    result = get_canonical_cds("BRCA1")

    assert result == "ATGAAATTT"
