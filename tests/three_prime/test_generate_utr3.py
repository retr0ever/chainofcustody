from chainofcustody.three_prime.generate_UTR3 import generate_mrna_sponge_utr


def test_generate_mrna_sponge_utr_basic() -> None:
    """generate_mrna_sponge_utr builds a valid 3'UTR with correct sites."""
    mir122_3p = "AACGCCAUUAUCACACUAAAUA"
    mir21_5p = "UAGCUUAUCAGACUGAUGUUGA"

    result = generate_mrna_sponge_utr([mir122_3p, mir21_5p], num_sites=4)

    assert len(result["single_sites"]) == 2
    assert "full_utr" in result

    site1 = result["single_sites"][0]
    site2 = result["single_sites"][1]
    # Sites should alternate in the cassette
    assert site1 in result["full_utr"]
    assert site2 in result["full_utr"]


def test_generate_mrna_sponge_utr_single_mirna_string() -> None:
    """Passing a single miRNA string (not a list) should work."""
    mir122_3p = "AACGCCAUUAUCACACUAAAUA"

    result = generate_mrna_sponge_utr(mir122_3p, num_sites=3)

    assert len(result["single_sites"]) == 1
    assert result["single_sites"][0] in result["full_utr"]
