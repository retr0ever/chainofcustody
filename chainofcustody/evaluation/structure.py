"""Metric 4: RNA secondary structure analysis via ViennaRNA."""

import RNA

from .parser import ParsedSequence


def fold_sequence(seq: str) -> tuple[str, float]:
    """Fold an RNA sequence. Returns (structure, MFE)."""
    structure, mfe = RNA.fold(seq)
    return structure, mfe


def check_utr5_accessibility(parsed: ParsedSequence) -> dict:
    """
    Check if the 5'UTR is accessible for ribosome loading.
    Strong secondary structure in the 5'UTR blocks the 43S pre-initiation complex.
    """
    utr5 = parsed.utr5
    if not utr5 or len(utr5) < 10:
        return {"mfe": None, "status": "no_utr5", "message": "No 5'UTR or too short to assess"}

    # Fold the 5'UTR + first 30nt of CDS (ribosome landing zone)
    landing_zone = utr5 + parsed.cds[:30]
    structure, mfe = fold_sequence(landing_zone)

    if mfe > -20:
        status = "GREEN"
        message = "5'UTR is accessible — low secondary structure"
    elif mfe > -30:
        status = "AMBER"
        message = "5'UTR has moderate structure — may partially impede ribosome loading"
    else:
        status = "RED"
        message = "5'UTR is highly structured — likely blocks ribosome scanning"

    return {
        "mfe": round(mfe, 2),
        "status": status,
        "message": message,
        "landing_zone_length": len(landing_zone),
    }


def check_mirna_site_accessibility(
    parsed: ParsedSequence,
    site_positions: list[int],
    site_length: int = 22,
    flank: int = 30,
) -> list[dict]:
    """
    Check if miRNA target sites are structurally accessible (not buried in hairpins).

    Args:
        site_positions: 0-indexed positions of miRNA sites in the full sequence.
        site_length: Length of the miRNA target site.
        flank: How many nt of context to include on each side for folding.
    """
    results = []
    seq = parsed.raw

    for pos in site_positions:
        # Extract local window around the site
        start = max(0, pos - flank)
        end = min(len(seq), pos + site_length + flank)
        window = seq[start:end]

        structure, mfe = fold_sequence(window)

        # Check if the seed region (first 8nt of site) is unpaired
        site_offset = pos - start
        seed_structure = structure[site_offset:site_offset + 8]
        paired_count = seed_structure.count("(") + seed_structure.count(")")
        unpaired_count = seed_structure.count(".")

        accessible = unpaired_count >= 5  # at least 5 of 8 seed positions unpaired

        results.append({
            "position": pos,
            "local_mfe": round(mfe, 2),
            "seed_structure": seed_structure,
            "seed_paired": paired_count,
            "seed_unpaired": unpaired_count,
            "accessible": accessible,
        })

    return results


def compute_global_mfe(parsed: ParsedSequence, max_length: int = 2000) -> dict:
    """
    Compute global MFE. For long sequences, fold in windows to avoid O(n^3) blowup.
    """
    seq = parsed.raw

    if len(seq) <= max_length:
        structure, mfe = fold_sequence(seq)
        return {
            "mfe": round(mfe, 2),
            "mfe_per_nt": round(mfe / len(seq), 4),
            "length": len(seq),
            "method": "full_fold",
        }

    # Window-based folding for long sequences
    window_size = 500
    step = 250
    mfe_values = []

    for i in range(0, len(seq) - window_size + 1, step):
        window = seq[i:i + window_size]
        _, window_mfe = fold_sequence(window)
        mfe_values.append(window_mfe)

    avg_mfe = sum(mfe_values) / len(mfe_values) if mfe_values else 0
    total_estimated_mfe = avg_mfe * (len(seq) / window_size)

    return {
        "mfe": round(total_estimated_mfe, 2),
        "mfe_per_nt": round(total_estimated_mfe / len(seq), 4),
        "length": len(seq),
        "method": "windowed_fold",
        "windows": len(mfe_values),
    }


def score_structure(parsed: ParsedSequence, mirna_site_positions: list[int] | None = None) -> dict:
    """Run all structure-related scoring."""
    result = {
        "utr5_accessibility": check_utr5_accessibility(parsed),
        "global_mfe": compute_global_mfe(parsed),
    }

    if mirna_site_positions:
        result["mirna_site_accessibility"] = check_mirna_site_accessibility(
            parsed, mirna_site_positions
        )

    return result
