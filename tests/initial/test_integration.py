import pytest

from chainofcustody.initial import get_canonical_cds


@pytest.mark.integration
def test_pou5f1_canonical_cds():
    cds = get_canonical_cds("POU5F1")

    assert len(cds) == 1083
    assert cds[:12] == "ATGGCGGGACAC"
    assert cds[-3:] == "TGA"  # stop codon
