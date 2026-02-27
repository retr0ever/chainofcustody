"""Metric 1: Codon quality + Metric 2: Cell-type selectivity scoring."""

import math
from collections import Counter

from python_codon_tables import get_codons_table

from .parser import ParsedSequence
from .data import get_liver_codon_weights, get_codon_weights

# Human codon usage table (NCBI Taxonomy ID 9606)
HUMAN_CODON_TABLE = get_codons_table("h_sapiens_9606")

# Genetic code: codon -> amino acid
CODON_TO_AA: dict[str, str] = {}
AA_TO_CODONS: dict[str, list[str]] = {}
for aa, codons_dict in HUMAN_CODON_TABLE.items():
    if aa == "*":
        continue
    for codon, freq in codons_dict.items():
        codon_dna = codon.upper().replace("U", "T")
        CODON_TO_AA[codon_dna] = aa
        AA_TO_CODONS.setdefault(aa, []).append(codon_dna)


def _relative_adaptiveness(codon_table: dict) -> dict[str, float]:
    """Compute relative adaptiveness (w_ij) for each codon: freq / max_freq_for_that_aa."""
    w = {}
    for aa, codons_dict in codon_table.items():
        if aa == "*":
            continue
        max_freq = max(codons_dict.values())
        if max_freq == 0:
            continue
        for codon, freq in codons_dict.items():
            codon_dna = codon.upper().replace("U", "T")
            w[codon_dna] = freq / max_freq if max_freq > 0 else 0
    return w


def compute_cai(parsed: ParsedSequence) -> float:
    """
    Codon Adaptation Index (Sharp & Li, 1987).
    Returns value between 0 and 1. Higher = better adapted to human codon usage.
    """
    w = _relative_adaptiveness(HUMAN_CODON_TABLE)
    codons = parsed.codons

    # Skip start codon (ATG) and stop codons
    scoring_codons = [c for c in codons if c in w and w[c] > 0]
    if not scoring_codons:
        return 0.0

    log_sum = sum(math.log(w[c]) for c in scoring_codons)
    return math.exp(log_sum / len(scoring_codons))


def compute_gc_content(parsed: ParsedSequence) -> dict[str, float]:
    """Compute GC content for full sequence and CDS separately."""
    def gc_pct(seq: str) -> float:
        if not seq:
            return 0.0
        gc = sum(1 for nt in seq if nt in "GC")
        return gc / len(seq) * 100

    return {
        "overall": round(gc_pct(parsed.raw), 1),
        "cds": round(gc_pct(parsed.cds), 1),
        "utr5": round(gc_pct(parsed.utr5), 1) if parsed.utr5 else None,
        "utr3": round(gc_pct(parsed.utr3), 1) if parsed.utr3 else None,
    }


def compute_rare_codon_clusters(parsed: ParsedSequence, window: int = 18, threshold: float = 0.1) -> list[dict]:
    """
    Find clusters of rare codons (%MinMax-style).
    Returns positions where a window of codons has unusually low average relative adaptiveness.
    """
    w = _relative_adaptiveness(HUMAN_CODON_TABLE)
    codons = parsed.codons
    clusters = []

    for i in range(len(codons) - window + 1):
        window_codons = codons[i:i + window]
        scores = [w.get(c, 0) for c in window_codons]
        avg = sum(scores) / len(scores) if scores else 0
        if avg < threshold:
            clusters.append({
                "codon_position": i,
                "nt_position": parsed.cds_start + i * 3,
                "avg_adaptiveness": round(avg, 4),
            })

    return clusters


def compute_liver_selectivity(parsed: ParsedSequence) -> float:
    """
    Score how well the sequence's codon choices correlate with liver cell TE.

    Returns average liver TE weight across all codons in the CDS.
    More negative = codons that are associated with lower liver translation (good for detargeting).
    More positive = codons that are associated with higher liver translation (bad for detargeting).
    """
    weights = get_liver_codon_weights()
    codons = parsed.codons
    scores = [weights.get(c, 0.0) for c in codons if c in weights]
    if not scores:
        return 0.0
    return sum(scores) / len(scores)


def compute_target_selectivity(parsed: ParsedSequence, target_cell_type: str) -> float:
    """
    Score how well the sequence's codon choices correlate with the target cell type's TE.
    Same logic as liver selectivity but for an arbitrary cell type column.
    """
    weights = get_codon_weights(target_cell_type)
    codons = parsed.codons
    scores = [weights.get(c, 0.0) for c in codons if c in weights]
    if not scores:
        return 0.0
    return sum(scores) / len(scores)


def score_codons(parsed: ParsedSequence, target_cell_type: str | None = None) -> dict:
    """Run all codon-related scoring. Returns a dict of results."""
    result = {
        "cai": round(compute_cai(parsed), 4),
        "gc_content": compute_gc_content(parsed),
        "rare_codon_clusters": compute_rare_codon_clusters(parsed),
        "liver_selectivity": round(compute_liver_selectivity(parsed), 4),
    }

    if target_cell_type:
        target_score = compute_target_selectivity(parsed, target_cell_type)
        result["target_selectivity"] = round(target_score, 4)
        liver_score = result["liver_selectivity"]
        if liver_score != 0:
            result["selectivity_ratio"] = round(target_score / liver_score, 4)

    return result
