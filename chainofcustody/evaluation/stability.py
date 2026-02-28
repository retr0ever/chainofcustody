"""Metric 6: mRNA stability scoring."""

import re

from chainofcustody.sequence import mRNASequence
from chainofcustody.evaluation.structure import fold_sequence, windowed_mfe_values


# AU-rich element motif (AUUUA pentamer)
ARE_PATTERN = re.compile(r"AUUUA")


def compute_gc3(parsed: mRNASequence) -> float:
    """
    GC content at the 3rd codon position (wobble position).
    Higher GC3 correlates with greater mRNA stability in mammals.
    """
    codons = parsed.codons
    if not codons:
        return 0.0
    gc3_count = sum(1 for c in codons if len(c) == 3 and c[2] in "GC")
    return gc3_count / len(codons)


def compute_mfe_per_nt(parsed: mRNASequence, max_length: int = 2000) -> float:
    """
    Minimum free energy per nucleotide. More negative = more stable.
    For long sequences, uses windowed folding.
    """
    seq = str(parsed)

    if len(seq) <= max_length:
        _, mfe = fold_sequence(seq)
        return mfe / len(seq) if seq else 0.0

    mfe_values = windowed_mfe_values(seq)
    if not mfe_values:
        return 0.0
    return sum(mfe_values) / len(mfe_values) / 500


def count_au_rich_elements(parsed: mRNASequence) -> int:
    """
    Count AU-rich elements (AREs) in the 3'UTR.
    AREs (AUUUA pentamer) destabilise mRNA by recruiting exosome degradation.
    Fewer = more stable.
    """
    if not parsed.utr3:
        return 0
    return len(ARE_PATTERN.findall(parsed.utr3))


def score_stability(parsed: mRNASequence) -> dict:
    """
    Compute mRNA stability metrics.

    Returns dict with:
    - gc3: GC content at wobble position (0-1)
    - mfe_per_nt: thermodynamic stability (kcal/mol/nt, negative = stable)
    - au_rich_elements: count of ARE motifs in 3'UTR (fewer = better)
    - stability_score: combined normalised score (0-1, higher = more stable)
    - status: GREEN/AMBER/RED traffic light
    """
    gc3 = compute_gc3(parsed)
    mfe_per_nt = compute_mfe_per_nt(parsed)
    are_count = count_au_rich_elements(parsed)

    # Normalise each sub-metric to 0-1
    # GC3: optimal around 0.5-0.7, penalise extremes
    if 0.5 <= gc3 <= 0.7:
        gc3_norm = 1.0
    elif gc3 < 0.5:
        gc3_norm = max(0.0, gc3 / 0.5)
    else:
        gc3_norm = max(0.0, (1.0 - gc3) / 0.3)

    # MFE/nt: more negative = more stable. Typical range -0.2 to -0.5
    # -0.4 or lower = very stable (1.0), above -0.1 = unstable (0.0)
    if mfe_per_nt <= -0.4:
        mfe_norm = 1.0
    elif mfe_per_nt >= -0.1:
        mfe_norm = 0.0
    else:
        mfe_norm = (-mfe_per_nt - 0.1) / 0.3

    # AREs: 0 = perfect, 3+ = bad
    are_norm = max(0.0, 1.0 - are_count / 3)

    # Combined score (weighted)
    stability_score = 0.4 * gc3_norm + 0.4 * mfe_norm + 0.2 * are_norm

    # Traffic light
    if stability_score >= 0.7:
        status = "GREEN"
    elif stability_score >= 0.4:
        status = "AMBER"
    else:
        status = "RED"

    return {
        "gc3": round(gc3, 4),
        "mfe_per_nt": round(mfe_per_nt, 4),
        "au_rich_elements": are_count,
        "stability_score": round(stability_score, 4),
        "status": status,
    }
