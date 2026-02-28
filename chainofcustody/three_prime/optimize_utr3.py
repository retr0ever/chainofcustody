"""Optimize only the spacer sequences in a 3'UTR sponge cassette using DnaChisel.

Everything is frozen (5'UTR, CDS, miRNA sponge sites, lead-in/out, poly-A signal)
except the short spacer segments between sponge sites.  Multiple variants are
produced by randomising the spacers before each DnaChisel pass.

All internal handling is RNA (A/U/G/C).  Conversion to DNA happens only at the
DnaChisel boundary and is converted back immediately.
"""

from __future__ import annotations

import random

from dnachisel import (
    AvoidChanges,
    AvoidHairpins,
    AvoidPattern,
    DnaOptimizationProblem,
    EnforceGCContent,
    Location,
)

from chainofcustody.evaluation.parser import clean_sequence, parse_sequence

# ---------------------------------------------------------------------------
# RNA ↔ DNA boundary helpers
# ---------------------------------------------------------------------------

_RNA_TO_DNA = str.maketrans("U", "T")
_DNA_TO_RNA = str.maketrans("T", "U")
_NUCLEOTIDES = "ACGT"  # DNA alphabet — used for randomisation before DnaChisel


def _rna_to_dna(seq: str) -> str:
    return seq.translate(_RNA_TO_DNA)


def _dna_to_rna(seq: str) -> str:
    return seq.translate(_DNA_TO_RNA)


# ---------------------------------------------------------------------------
# Spacer position detection
# ---------------------------------------------------------------------------

def _find_spacer_positions(
    full_sequence: str,
    utr3_start: int,
    sponge_sites: list[str],
    num_sites: int,
) -> list[tuple[int, int]]:
    """Return (start, end) positions of each spacer in *full_sequence*.

    Walks the 3'UTR cassette layout:
        stop_codon + lead_in + [site spacer site spacer … site] + lead_out + …

    The sponge sites are matched in order; everything between two consecutive
    site occurrences that is not a site is a spacer.
    """
    # The cassette starts after stop_codon(3) + lead_in(6) inside the 3'UTR
    cassette_offset = utr3_start + 3 + 6  # "UAA" + "gcauac"

    spacer_positions: list[tuple[int, int]] = []
    cursor = cassette_offset

    for i in range(num_sites):
        site = sponge_sites[i % len(sponge_sites)].upper()
        site_len = len(site)
        # Skip over the sponge site
        cursor += site_len

        if i < num_sites - 1:
            # Next site tells us where the spacer ends
            next_site = sponge_sites[(i + 1) % len(sponge_sites)].upper()
            next_site_pos = full_sequence.upper().find(next_site, cursor)
            if next_site_pos == -1:
                # Fallback: assume fixed 4-nt spacer (default from generate_UTR3)
                spacer_positions.append((cursor, cursor + 4))
                cursor += 4
            else:
                spacer_positions.append((cursor, next_site_pos))
                cursor = next_site_pos

    return spacer_positions


# ---------------------------------------------------------------------------
# Core optimiser
# ---------------------------------------------------------------------------

def optimize_utr3_spacers(
    sequence: str,
    sponge_sites: list[str],
    num_sites: int = 16,
    n_variants: int = 5,
    gc_min: float = 0.25,
    gc_max: float = 0.75,
    utr5_end: int | None = None,
    cds_end: int | None = None,
    seed: int | None = None,
) -> list[str]:
    """Optimize only the spacers between miRNA sponge sites with DnaChisel.

    Everything else — 5'UTR, CDS, sponge sites, lead-in/out, poly-A signal —
    is frozen via ``AvoidChanges``.

    Parameters
    ----------
    sequence : str
        Full mRNA sequence (DNA or RNA) including a 3'UTR built by
        ``generate_mrna_sponge_utr``.
    sponge_sites : list[str]
        The ``single_sites`` list returned by ``generate_mrna_sponge_utr``.
    num_sites : int
        Number of sponge sites used when building the 3'UTR.
    n_variants : int
        Number of optimised output sequences to produce.
    gc_min, gc_max : float
        GC bounds applied to each spacer region.
    utr5_end, cds_end : int | None
        Optional manual boundary overrides forwarded to the parser.
    seed : int | None
        Random seed for reproducibility.

    Returns
    -------
    list[str]
        ``n_variants`` full RNA sequences with optimised spacers.
    """
    rng = random.Random(seed)

    # ── 1. Parse to find the 3'UTR boundary ──────────────────────────────
    rna_seq = clean_sequence(sequence)
    parsed = parse_sequence(rna_seq, utr5_end=utr5_end, cds_end=cds_end)
    utr3_start = len(parsed.utr5) + len(parsed.cds)

    if not parsed.utr3:
        raise ValueError("Sequence has no 3'UTR region to optimise")

    # ── 2. Locate spacer positions ───────────────────────────────────────
    spacer_spans = _find_spacer_positions(
        rna_seq, utr3_start, sponge_sites, num_sites,
    )
    if not spacer_spans:
        raise ValueError("No spacers found in the 3'UTR cassette")

    total_spacer_nt = sum(end - start for start, end in spacer_spans)

    # ── 3. Build frozen indices (everything that is NOT a spacer) ────────
    spacer_indices: set[int] = set()
    for start, end in spacer_spans:
        spacer_indices.update(range(start, end))

    frozen_indices = sorted(set(range(len(rna_seq))) - spacer_indices)

    # ── 4. Build DnaChisel constraints & objectives ──────────────────────
    constraints = [
        AvoidChanges(indices=frozen_indices),
    ]
    # Per-spacer constraints
    for start, end in spacer_spans:
        loc = Location(start, end)
        constraints.extend([
            EnforceGCContent(mini=gc_min, maxi=gc_max, location=loc),
            AvoidPattern("AAAA", location=loc),
            AvoidPattern("TTTT", location=loc),
            AvoidPattern("GGGG", location=loc),
            AvoidPattern("CCCC", location=loc),
        ])

    objectives = [
        AvoidHairpins(stem_size=4, hairpin_window=60),
    ]

    # ── 5. Generate N variants ───────────────────────────────────────────
    full_dna = _rna_to_dna(rna_seq)
    variants: list[str] = []

    for _ in range(n_variants):
        # Start from original sequence, randomise only spacer positions
        seq_list = list(full_dna)
        for start, end in spacer_spans:
            for j in range(start, end):
                seq_list[j] = rng.choice(_NUCLEOTIDES)
        candidate_dna = "".join(seq_list)

        problem = DnaOptimizationProblem(
            sequence=candidate_dna,
            constraints=constraints,
            objectives=objectives,
        )
        problem.resolve_constraints()
        problem.optimize()

        result_rna = _dna_to_rna(problem.sequence)
        variants.append(result_rna)

    return variants
