"""Tests for three_prime.optimize_utr3 — DnaChisel spacer-only optimisation."""

import pytest

from chainofcustody.three_prime.generate_UTR3 import generate_mrna_sponge_utr
from chainofcustody.three_prime.optimize_utr3 import optimize_utr3_spacers


# Build a reusable test mRNA with a 3'UTR sponge cassette
MOCK_5UTR = "GGG"
MOCK_CDS = "AUGAAAGCUGCUUAA"  # Met-Lys-Ala-Ala-Stop
MIR122 = "AACGCCAUUAUCACACUAAAUA"
NUM_SITES = 4

_sponge = generate_mrna_sponge_utr(MIR122, num_sites=NUM_SITES)
SPONGE_SITES = _sponge["single_sites"]
UTR3 = _sponge["full_utr"]
FULL_SEQ = MOCK_5UTR + MOCK_CDS + UTR3
PREFIX = MOCK_5UTR + MOCK_CDS


def test_returns_requested_variant_count() -> None:
    results = optimize_utr3_spacers(
        FULL_SEQ, sponge_sites=SPONGE_SITES, num_sites=NUM_SITES, n_variants=3,
    )
    assert len(results) == 3
    for seq in results:
        assert isinstance(seq, str)
        assert "T" not in seq, "Output must be RNA (no T)"


def test_preserves_frozen_prefix() -> None:
    results = optimize_utr3_spacers(
        FULL_SEQ, sponge_sites=SPONGE_SITES, num_sites=NUM_SITES, n_variants=3,
    )
    for seq in results:
        assert seq[: len(PREFIX)] == PREFIX, "5'UTR + CDS was modified"


def test_preserves_sponge_sites() -> None:
    """Every miRNA sponge site must appear intact in every variant."""
    results = optimize_utr3_spacers(
        FULL_SEQ, sponge_sites=SPONGE_SITES, num_sites=NUM_SITES, n_variants=3,
    )
    site = SPONGE_SITES[0].upper()
    for seq in results:
        assert seq.upper().count(site) == NUM_SITES, (
            f"Expected {NUM_SITES} occurrences of sponge site, "
            f"got {seq.upper().count(site)}"
        )


def test_raises_on_no_utr3() -> None:
    no_utr3 = "AUGAAAUAA"
    with pytest.raises(ValueError, match="no 3'UTR"):
        optimize_utr3_spacers(
            no_utr3, sponge_sites=SPONGE_SITES, num_sites=NUM_SITES, n_variants=1,
        )
