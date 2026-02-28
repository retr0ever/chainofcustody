"""Tests for chainofcustody.three_prime.generate_utr3."""

import pytest

from chainofcustody.three_prime.generate_utr3 import generate_mrna_sponge_utr

# ── helpers ──────────────────────────────────────────────────────────────────

_MIR122 = "UGGAGUGUGACAAUGGUGUUUG"  # 22-nt example
_MIR21  = "UAGCUUAUCAGACUGAUGUUGA"  # 22-nt example
_RNA_BASES = set("AUGC")


def _is_rna(seq: str) -> bool:
    return set(seq.upper()) <= _RNA_BASES


# ── generate_mrna_sponge_utr ─────────────────────────────────────────────────

class TestGenerateMrnaSpongeUtr:

    def test_returns_dict_with_expected_keys(self):
        result = generate_mrna_sponge_utr([_MIR122])
        assert set(result.keys()) == {"single_sites", "full_utr"}

    def test_single_site_count_matches_input(self):
        result = generate_mrna_sponge_utr([_MIR122, _MIR21])
        assert len(result["single_sites"]) == 2

    def test_full_utr_is_non_empty_string(self):
        result = generate_mrna_sponge_utr([_MIR122])
        assert isinstance(result["full_utr"], str)
        assert len(result["full_utr"]) > 0

    def test_full_utr_starts_with_stop_codon(self):
        result = generate_mrna_sponge_utr([_MIR122])
        assert result["full_utr"].startswith("UAA")

    def test_full_utr_is_rna(self):
        result = generate_mrna_sponge_utr([_MIR122])
        # lowercase spacers + uppercase sites — normalise before checking
        assert _is_rna(result["full_utr"].upper())

    def test_single_sites_are_rna(self):
        result = generate_mrna_sponge_utr([_MIR122, _MIR21])
        for site in result["single_sites"]:
            assert _is_rna(site), f"Non-RNA characters in site: {site!r}"

    def test_num_sites_controls_cassette_count(self):
        """More sites → longer full UTR."""
        short = generate_mrna_sponge_utr([_MIR122], num_sites=4)
        long_ = generate_mrna_sponge_utr([_MIR122], num_sites=16)
        assert len(long_["full_utr"]) > len(short["full_utr"])

    def test_single_string_input_accepted(self):
        """Passing a bare string instead of a list should work."""
        result_str  = generate_mrna_sponge_utr(_MIR122)
        result_list = generate_mrna_sponge_utr([_MIR122])
        assert result_str["full_utr"] == result_list["full_utr"]

    def test_single_site_length_equals_mirna_length(self):
        """Site = rc(miRNA) with bulge mutated → same length as input miRNA."""
        result = generate_mrna_sponge_utr([_MIR122])
        assert len(result["single_sites"][0]) == len(_MIR122)

    def test_multiple_mirnas_alternate_in_cassette(self):
        """With 2 miRNAs and 4 sites the cassette should contain each site ×2."""
        result = generate_mrna_sponge_utr([_MIR122, _MIR21], num_sites=4)
        site0 = result["single_sites"][0]
        site1 = result["single_sites"][1]
        assert result["full_utr"].count(site0) == 2
        assert result["full_utr"].count(site1) == 2

    def test_dna_input_is_normalised_to_rna(self):
        """T in input should be silently converted to U."""
        dna_seq = _MIR122.replace("U", "T")
        result_dna = generate_mrna_sponge_utr(dna_seq)
        result_rna = generate_mrna_sponge_utr(_MIR122)
        assert result_dna["full_utr"] == result_rna["full_utr"]

    def test_num_sites_one(self):
        result = generate_mrna_sponge_utr([_MIR122], num_sites=1)
        assert result["full_utr"].count(result["single_sites"][0]) == 1

    def test_seed_match_embedded_in_site(self):
        """The last 8 nt of the reverse complement should appear in the site."""
        complement = {"A": "U", "U": "A", "G": "C", "C": "G"}
        rc = "".join(complement[b] for b in reversed(_MIR122))
        seed_match = rc[-8:]
        result = generate_mrna_sponge_utr([_MIR122])
        assert seed_match in result["single_sites"][0]

    def test_empty_sequence_raises_value_error(self):
        with pytest.raises(ValueError, match="too short"):
            generate_mrna_sponge_utr([""])

    def test_too_short_sequence_raises_value_error(self):
        """Sequences under 12 nt cannot form seed, bulge, and 3'-match regions."""
        with pytest.raises(ValueError, match="too short"):
            generate_mrna_sponge_utr(["AUGCAUGC"])  # 8 nt
