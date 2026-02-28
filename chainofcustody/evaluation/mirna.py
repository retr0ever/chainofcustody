"""Metric 3: miRNA detargeting site scanning."""

import re
from dataclasses import dataclass

from chainofcustody.sequence import mRNASequence

# miRNA library: name -> mature sequence (RNA, 5'->3')
MIRNA_LIBRARY = {
    "miR-122-5p": "UGGAGUGUGACAAUGGUGUUUG",
    "miR-1-3p": "UGGAAUGUAAAGAAGUAUGUAU",
    "miR-208a-3p": "AUAAGACGAGCAAAAAGCUUGU",
    "miR-142-3p": "UGUAGUGUUUCCUACUUUAUGGA",
}

# Which miRNAs are "desired" (for detargeting) vs "accidental" (should NOT be present)
DESIRED_MIRNAS = {"miR-122-5p"}  # liver detargeting
WARNING_MIRNAS = {"miR-1-3p", "miR-208a-3p", "miR-142-3p"}  # accidental matches


@dataclass
class MirnaSiteHit:
    """A detected miRNA target site in the sequence."""
    mirna_name: str
    position: int  # 0-indexed position in full sequence
    region: str  # "5utr", "cds", or "3utr"
    match_type: str  # "8mer", "7mer-m8", "7mer-a1", "6mer", "full"
    seed_seq: str  # the matched seed sequence


def _reverse_complement(seq: str) -> str:
    """Reverse complement of an RNA sequence."""
    comp = {"A": "U", "U": "A", "G": "C", "C": "G"}
    return "".join(comp[nt] for nt in reversed(seq))


def _get_seed_target(mirna_rna: str) -> str:
    """
    Get the RNA target sequence for the miRNA seed region (nt 2-8).
    This is the reverse complement of the seed, in RNA.
    """
    seed_rna = mirna_rna[1:8]  # positions 2-8 (0-indexed: 1-7)
    return _reverse_complement(seed_rna)


def _get_full_target(mirna_rna: str) -> str:
    """Get the full reverse complement target sequence in RNA."""
    return _reverse_complement(mirna_rna)


def _classify_region(pos: int, parsed: mRNASequence) -> str:
    """Classify a position as 5'UTR, CDS, or 3'UTR."""
    if pos < parsed.cds_start:
        return "5utr"
    elif pos < parsed.cds_end:
        return "cds"
    else:
        return "3utr"


def scan_for_mirna(
    parsed: mRNASequence,
    mirna_name: str,
    mirna_seq: str,
) -> list[MirnaSiteHit]:
    """Scan sequence for target sites of a specific miRNA."""
    hits = []
    seq = str(parsed)

    # Full complementary match (most potent — triggers RISC cleavage)
    full_target = _get_full_target(mirna_seq)
    for m in re.finditer(re.escape(full_target), seq):
        hits.append(MirnaSiteHit(
            mirna_name=mirna_name,
            position=m.start(),
            region=_classify_region(m.start(), parsed),
            match_type="full",
            seed_seq=full_target,
        ))

    # Seed-based matches (7mer) — nt 2-8 of miRNA
    seed_target = _get_seed_target(mirna_seq)
    for m in re.finditer(re.escape(seed_target), seq):
        # Don't double-count if this is part of a full match
        is_in_full = any(
            h.match_type == "full" and h.position <= m.start() <= h.position + len(full_target)
            for h in hits
        )
        if not is_in_full:
            hits.append(MirnaSiteHit(
                mirna_name=mirna_name,
                position=m.start(),
                region=_classify_region(m.start(), parsed),
                match_type="7mer-seed",
                seed_seq=seed_target,
            ))

    return hits


def compute_site_spacing(hits: list[MirnaSiteHit]) -> list[int]:
    """Compute spacing (in nt) between consecutive miRNA sites."""
    if len(hits) < 2:
        return []
    positions = sorted(h.position for h in hits)
    return [positions[i+1] - positions[i] for i in range(len(positions) - 1)]


def score_mirna(parsed: mRNASequence) -> dict:
    """
    Run miRNA target site scanning.

    Returns dict with:
    - detargeting: miR-122 site details (count, positions, spacing, region breakdown)
    - warnings: accidental matches for target cell type miRNAs
    """
    result = {"detargeting": {}, "warnings": []}

    # Scan for desired detargeting miRNAs (miR-122)
    for mirna_name in DESIRED_MIRNAS:
        mirna_seq = MIRNA_LIBRARY[mirna_name]
        hits = scan_for_mirna(parsed, mirna_name, mirna_seq)

        utr3_hits = [h for h in hits if h.region == "3utr"]
        other_hits = [h for h in hits if h.region != "3utr"]
        spacing = compute_site_spacing(utr3_hits)

        result["detargeting"][mirna_name] = {
            "total_sites": len(hits),
            "utr3_sites": len(utr3_hits),
            "other_sites": len(other_hits),
            "positions": [h.position for h in hits],
            "utr3_positions": [h.position for h in utr3_hits],
            "match_types": [h.match_type for h in hits],
            "inter_site_spacing": spacing,
            "adequate_spacing": all(s >= 8 for s in spacing) if spacing else True,
        }

    # Scan for accidental matches with target cell type miRNAs
    for mirna_name in WARNING_MIRNAS:
        mirna_seq = MIRNA_LIBRARY[mirna_name]
        hits = scan_for_mirna(parsed, mirna_name, mirna_seq)
        if hits:
            result["warnings"].append({
                "mirna": mirna_name,
                "sites_found": len(hits),
                "positions": [h.position for h in hits],
                "regions": [h.region for h in hits],
                "message": f"Accidental {mirna_name} target site(s) detected — may cause unwanted silencing in target cell type",
            })

    return result
