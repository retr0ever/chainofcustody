"""3'UTR mRNA sponge sequence generator."""

from __future__ import annotations


def generate_mrna_sponge_utr(
    mirna_sequences: list[str] | str,
    num_sites: int = 16,
) -> dict[str, list[str] | str]:
    """Generate a 3'UTR sequence with alternating, bulged miRNA sponge sites.

    Parameters
    ----------
    mirna_sequences : list[str] or str
        Mature miRNA sequences (5' to 3') in RNA alphabet (A/U/G/C).
        A single string is accepted for backward compatibility.
    num_sites : int
        Total number of sponge site repeats. Default is 16.

    Returns
    -------
    dict
        ``single_sites`` — list of per-miRNA bulged site sequences.
        ``full_utr``     — assembled 3'UTR string (RNA, mixed case) including
                           stop codon, sponge cassette, and poly-A signal.

    Raises
    ------
    ValueError
        If any miRNA sequence is shorter than 12 nucleotides (too short to
        decompose into seed, bulge, and 3′-match regions).
    """
    if isinstance(mirna_sequences, str):
        mirna_sequences = [mirna_sequences]

    def _reverse_complement(seq: str) -> str:
        complement = {"A": "U", "U": "A", "G": "C", "C": "G"}
        return "".join(complement[base] for base in reversed(seq))

    def _create_mismatch(seq: str) -> str:
        mismatch_map = {"A": "C", "U": "G", "G": "U", "C": "A"}
        return "".join(mismatch_map[base] for base in seq)

    sponge_sites: list[str] = []
    for mirna_seq in mirna_sequences:
        mirna = mirna_seq.upper().replace("T", "U")
        if len(mirna) < 12:
            raise ValueError(
                f"miRNA sequence too short ({len(mirna)} nt): '{mirna}'. "
                "At least 12 nucleotides are required to form seed, bulge, "
                "and 3′-match regions."
            )
        rc_mirna = _reverse_complement(mirna)

        seed_match = rc_mirna[-8:]
        bulge_rc = rc_mirna[-12:-8]
        three_prime_match = rc_mirna[:-12]

        bulge_mismatch = _create_mismatch(bulge_rc)

        site = three_prime_match + bulge_mismatch + seed_match
        sponge_sites.append(site)

    spacers = [
        "aauu", "ucga", "caag", "auac", "gaau",
        "cuua", "uuca", "agcu", "uacg", "gaua",
        "cuac", "acuc", "uguu", "caua", "ucuu", "agau",
    ]

    cassette = ""
    for i in range(num_sites):
        current_site = sponge_sites[i % len(sponge_sites)]
        cassette += current_site
        if i < num_sites - 1:
            cassette += spacers[i % len(spacers)]

    stop_codon = "UAA"
    lead_in = "gcauac"
    lead_out = "gauc"
    # Synthetic poly-A signal with upstream regulatory elements; sequence
    # derived from empirical mRNA stabilisation constructs.
    poly_a_signal = (
        "CUCAGGUGCAGGCUGCCUAUCAGAAGGUGGUGGCUGGUGUGGCCAAUGCCCUGGCUCACAAAUACCACUGAGAUC"
        "UUUUUCCCUCUGCCAAAAAUUAUGGGGACAUCAUGAAGCCCCUUGAGCAUCUGACUUCUGGCUAAUAAAGGAAAU"
        "UUAUUUUCAUUGCAAUAGUGUGUUGGAAUUUUUUGUGUCUCUCACUCGGAAGGACAUAUGGGAGGGCAAAUCAUU"
        "UAAAACAUCAGAAUGAGUAUUUGGUUUAGAGUUUGGCA"
    )

    final_utr = f"{stop_codon}{lead_in}{cassette}{lead_out}{poly_a_signal}"

    return {
        "single_sites": sponge_sites,
        "full_utr": final_utr,
    }
