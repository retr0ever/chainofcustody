"""Metric 1: RNA secondary structure analysis via ViennaRNA."""

from __future__ import annotations

import RNA

from chainofcustody.sequence import mRNASequence

# Cap for global-MFE folds during batch scoring.  Only the 5'UTR (≤100 nt)
# varies between optimizer candidates; the CDS and 3'UTR are fixed.  Folding
# the first 150 nt covers the entire variable region (5'UTR) plus the CDS
# start codon context while keeping each fold at O(150³) ≈ 12 ms instead
# of O(1800³) ≈ 4 s for a full mRNA.  For the final single-sequence report
# the full fold is used automatically (seq <= 2000 nt path in compute_global_mfe).
_GLOBAL_FOLD_CAP = 150

# When the 5'UTR is long (e.g. 1000 nt), folding the entire region becomes
# O(n³) ≈ several seconds per sequence.  Only the region immediately upstream
# of the AUG matters for ribosome scanning, so we fold the last
# _UTR5_FOLD_WINDOW nt (i.e. the AUG-proximal end) and score that.
_UTR5_FOLD_WINDOW = 200


def fold_sequence(seq: str) -> tuple[str, float]:
    """Fold an RNA sequence. Returns ``(dot_bracket, mfe_kcal_mol)``."""
    structure, mfe = RNA.fold(seq)
    return structure, float(mfe)


def fold_sequence_bounded(seq: str, cap: int = _GLOBAL_FOLD_CAP) -> tuple[str, float]:
    """Fold up to *cap* nt of *seq*, return ``(dot_bracket, mfe_kcal_mol)``.

    Used for global-MFE computations during batch optimisation where the
    variable region (5'UTR) is always within the first few hundred nt.
    Returns a length-scaled pseudo-MFE so callers can treat it like the
    full-sequence result.
    """
    if len(seq) <= cap:
        return fold_sequence(seq)
    structure, mfe = fold_sequence(seq[:cap])
    # Scale MFE linearly to the full sequence length for downstream normalisation
    scaled_mfe = mfe * len(seq) / cap
    return structure, scaled_mfe


def windowed_mfe_values(
    seq: str,
    window_size: int = 500,
    step: int = 250,
) -> list[float]:
    """Fold a sequence in overlapping windows and return each window's MFE."""
    return [
        fold_sequence(seq[i:i + window_size])[1]
        for i in range(0, len(seq) - window_size + 1, step)
    ]


def check_utr5_accessibility(parsed: mRNASequence) -> dict:
    """Check if the 5'UTR is accessible for ribosome loading.

    Folds the AUG-proximal end of the 5'UTR (last _UTR5_FOLD_WINDOW nt,
    including Kozak) and returns MFE/nt.  For short UTRs the full sequence is
    folded.  Limiting the window keeps ViennaRNA tractable even when the evolved
    5'UTR reaches 1000 nt, while still scoring the region that matters for
    ribosome scanning.  Excluding the CDS removes the dominant contribution of
    the fixed GC-rich start codon context.

    A less negative MFE/nt indicates a more open, accessible structure.
    """
    utr5 = parsed.utr5
    if not utr5 or len(utr5) < 10:
        return {
            "mfe": None,
            "mfe_per_nt": None,
            "status": "no_utr5",
            "message": "No 5'UTR or too short to assess",
        }

    fold_region = utr5[-_UTR5_FOLD_WINDOW:] if len(utr5) > _UTR5_FOLD_WINDOW else utr5
    structure, mfe = fold_sequence(fold_region)
    mfe_per_nt = mfe / len(fold_region)

    if mfe_per_nt >= -0.1:
        status = "GREEN"
        message = "5'UTR is accessible — weak secondary structure"
    elif mfe_per_nt >= -0.3:
        status = "AMBER"
        message = "5'UTR has moderate secondary structure"
    else:
        status = "RED"
        message = "5'UTR is highly structured — may impede ribosome scanning"

    return {
        "mfe": round(mfe, 2),
        "mfe_per_nt": round(mfe_per_nt, 4),
        "utr5_length": len(utr5),
        "fold_window": len(fold_region),
        "status": status,
        "message": message,
    }


def check_mirna_site_accessibility(
    parsed: mRNASequence,
    site_positions: list[int],
    site_length: int = 22,
    flank: int = 30,
) -> list[dict]:
    """Check if miRNA target sites are structurally accessible.

    Args:
        site_positions: 0-indexed positions of miRNA sites in the full sequence.
        site_length: Length of the miRNA target site.
        flank: How many nt of context to include on each side for folding.
    """
    results = []
    seq = str(parsed)

    for pos in site_positions:
        start = max(0, pos - flank)
        end = min(len(seq), pos + site_length + flank)
        window = seq[start:end]

        structure, mfe = fold_sequence(window)

        site_offset = pos - start
        seed_structure = structure[site_offset:site_offset + 8]
        paired_count = seed_structure.count("(") + seed_structure.count(")")
        unpaired_count = seed_structure.count(".")

        accessible = unpaired_count >= 5

        results.append({
            "position": pos,
            "local_mfe": round(mfe, 2),
            "seed_structure": seed_structure,
            "seed_paired": paired_count,
            "seed_unpaired": unpaired_count,
            "accessible": accessible,
        })

    return results


def compute_global_mfe(
    parsed: mRNASequence,
    max_length: int = 2000,
    _precomputed: tuple[str, float] | None = None,
) -> dict:
    """Compute the global MFE of the full mRNA sequence.

    If *_precomputed* is provided it is used directly, avoiding a second fold.
    For sequences longer than *max_length* nt (and no pre-computed value),
    folds in overlapping windows to avoid quadratic memory growth.
    """
    seq = str(parsed)

    if _precomputed is not None:
        structure, mfe = _precomputed
        return {
            "mfe": round(mfe, 2),
            "mfe_per_nt": round(mfe / len(seq), 4) if seq else 0.0,
            "length": len(seq),
            "method": "precomputed",
        }

    if len(seq) <= max_length:
        structure, mfe = fold_sequence(seq)
        return {
            "mfe": round(mfe, 2),
            "mfe_per_nt": round(mfe / len(seq), 4),
            "length": len(seq),
            "method": "full_fold",
        }

    mfe_values = windowed_mfe_values(seq)
    avg_mfe = sum(mfe_values) / len(mfe_values) if mfe_values else 0
    total_estimated_mfe = avg_mfe * (len(seq) / 500)

    return {
        "mfe": round(total_estimated_mfe, 2),
        "mfe_per_nt": round(total_estimated_mfe / len(seq), 4),
        "length": len(seq),
        "method": "windowed_fold",
        "windows": len(mfe_values),
    }


def score_structure(
    parsed: mRNASequence,
    mirna_site_positions: list[int] | None = None,
    _precomputed_global: tuple[str, float] | None = None,
) -> dict:
    """Run all structure-related scoring.

    *_precomputed_global* is an optional ``(dot_bracket, mfe)`` tuple for the
    full sequence — when supplied ``compute_global_mfe`` skips a second fold.
    """
    result = {
        "utr5_accessibility": check_utr5_accessibility(parsed),
        "global_mfe": compute_global_mfe(parsed, _precomputed=_precomputed_global),
    }

    if mirna_site_positions:
        result["mirna_site_accessibility"] = check_mirna_site_accessibility(
            parsed, mirna_site_positions
        )

    return result
