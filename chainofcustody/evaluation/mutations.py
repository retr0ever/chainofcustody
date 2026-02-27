"""Mutation operators for evolutionary optimisation. All preserve the protein sequence."""

import random

from .parser import ParsedSequence, parse_sequence, STOP_CODONS
from .codons import CODON_TO_AA, AA_TO_CODONS, HUMAN_CODON_TABLE, _relative_adaptiveness
from .mirna import _get_seed_target, MIRNA_LIBRARY


def synonymous_swap(seq: str, parsed: ParsedSequence, rate: float = 0.05) -> str:
    """
    Replace random codons with synonymous alternatives.
    Preserves the protein sequence. Biased toward higher-CAI codons.

    Args:
        seq: Full DNA sequence.
        parsed: Parsed sequence with region boundaries.
        rate: Fraction of codons to mutate (0.0 to 1.0).

    Returns: New sequence with codon swaps in the CDS.
    """
    codons = list(parsed.codons)
    w = _relative_adaptiveness(HUMAN_CODON_TABLE)
    n_to_swap = max(1, int(len(codons) * rate))
    indices = random.sample(range(len(codons)), min(n_to_swap, len(codons)))

    for idx in indices:
        codon = codons[idx]
        aa = CODON_TO_AA.get(codon)
        if not aa:
            continue
        synonyms = [c for c in AA_TO_CODONS.get(aa, []) if c != codon and c not in STOP_CODONS]
        if not synonyms:
            continue
        # Weight toward higher adaptiveness
        weights = [max(w.get(c, 0.01), 0.01) for c in synonyms]
        codons[idx] = random.choices(synonyms, weights=weights, k=1)[0]

    new_cds = "".join(codons)
    return parsed.utr5 + new_cds + parsed.utr3


def crossover(seq_a: str, seq_b: str) -> tuple[str, str]:
    """
    Region-level crossover between two parent sequences.
    Produces two children by swapping CDS halves (codon-aligned).

    Both parents must have the same protein sequence (same CDS length).
    Falls back to full-CDS swap if lengths differ.
    """
    parsed_a = parse_sequence(seq_a)
    parsed_b = parse_sequence(seq_b)

    if len(parsed_a.cds) == len(parsed_b.cds):
        # Codon-aligned crossover at a random point in the CDS
        n_codons = len(parsed_a.codons)
        xpoint = random.randint(1, n_codons - 1)
        nt_point = xpoint * 3

        cds_child1 = parsed_a.cds[:nt_point] + parsed_b.cds[nt_point:]
        cds_child2 = parsed_b.cds[:nt_point] + parsed_a.cds[nt_point:]

        child1 = parsed_a.utr5 + cds_child1 + parsed_b.utr3
        child2 = parsed_b.utr5 + cds_child2 + parsed_a.utr3
    else:
        # Different CDS lengths — swap whole regions
        child1 = parsed_a.utr5 + parsed_a.cds + parsed_b.utr3
        child2 = parsed_b.utr5 + parsed_b.cds + parsed_a.utr3

    return child1, child2


def insert_mir122_sites(seq: str, parsed: ParsedSequence, n_sites: int = 3, spacing: int = 12) -> str:
    """
    Insert miR-122-5p seed match sites into the 3'UTR.
    If the 3'UTR already has sites, only adds enough to reach n_sites.

    Args:
        seq: Full DNA sequence.
        parsed: Parsed sequence.
        n_sites: Target number of miR-122 sites in 3'UTR.
        spacing: Nucleotides between sites.

    Returns: Sequence with miR-122 sites inserted at end of 3'UTR.
    """
    # miR-122 seed target (reverse complement of nt 2-8)
    mir122_rna = MIRNA_LIBRARY["miR-122-5p"]
    seed_target = _get_seed_target(mir122_rna)

    # Count existing sites in 3'UTR
    existing = parsed.utr3.count(seed_target)
    to_add = max(0, n_sites - existing)

    if to_add == 0:
        return seq

    # Build insert: seed sites with spacer between them
    spacer = "A" * spacing
    insert = spacer.join([seed_target] * to_add)

    # Insert before the last 10nt of 3'UTR (or at end if short)
    utr3 = parsed.utr3
    if len(utr3) > 10:
        new_utr3 = utr3[:-10] + spacer + insert + spacer + utr3[-10:]
    else:
        new_utr3 = utr3 + spacer + insert

    return parsed.utr5 + parsed.cds + new_utr3
