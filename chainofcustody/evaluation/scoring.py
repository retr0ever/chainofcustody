"""Orchestrate the full 4-metric scoring pipeline over a parsed mRNA sequence."""

from chainofcustody.sequence import CAP5, POLY_A_LENGTH, mRNASequence
from chainofcustody.evaluation.structure import fold_sequence_bounded, fold_sequence, score_structure
from chainofcustody.evaluation.manufacturing import score_manufacturing
from chainofcustody.evaluation.stability import score_stability
from chainofcustody.evaluation.ribonn import score_ribonn


def score_parsed(
    parsed: mRNASequence,
    target: str | None = None,
    _ribonn_scores: dict | None = None,
    _fast_fold: bool = False,
    target_cell_type: str = "megakaryocytes",
) -> dict:
    """Run all 4 evaluation metrics on an already-parsed mRNA sequence.

    The full-sequence fold is computed **once** and shared between
    :func:`score_structure` (global MFE) and :func:`score_stability` (MFE/nt),
    avoiding two independent O(n³) ViennaRNA calls on the same sequence.

    Args:
        parsed: An ``mRNASequence`` whose regions are correctly delimited.
        target: Unused — reserved for future per-cell-type scoring.
        _ribonn_scores: Pre-computed RiboNN result dict (batch optimiser path).
        _fast_fold: When True, caps the global fold to 150 nt (covers the
            entire variable 5'UTR region).  Used by the batch optimiser to
            keep each fold at ~12 ms instead of ~4 s.  The final per-candidate
            report always uses the full fold (default False).

    Returns:
        Full report dict with keys: ``sequence_info``, ``structure_scores``,
        ``manufacturing_scores``, ``stability_scores``, ``ribonn_scores``,
        ``summary``.
    """
    # Fold ONCE, share between structure and stability to avoid duplicate work.
    seq = str(parsed)
    if _fast_fold:
        global_fold = fold_sequence_bounded(seq)   # caps at _GLOBAL_FOLD_CAP (150 nt)
    else:
        structure_str, mfe = fold_sequence(seq)
        global_fold = (structure_str, mfe)

    structure_scores = score_structure(parsed, _precomputed_global=global_fold)
    manufacturing_scores = score_manufacturing(parsed)
    stability_scores = score_stability(parsed, _precomputed_mfe=global_fold[1])
    ribonn_scores = _ribonn_scores if _ribonn_scores is not None else score_ribonn(parsed, target_cell_type=target_cell_type)

    mfg_violations = manufacturing_scores.get("total_violations", 0)

    summary = {
        "utr5_accessibility": structure_scores.get("utr5_accessibility", {}).get("status", "GREY"),
        "manufacturability": _traffic_light(-mfg_violations, (-3, 0), (-999, 0)),
        "stability": stability_scores.get("status", "GREY"),
        "specificity": ribonn_scores.get("status", "GREY"),
    }

    return {
        "sequence_info": {
            "total_length": len(parsed),
            "full_length": parsed.full_length,
            "cap5_length": len(CAP5),
            "utr5_length": len(parsed.utr5),
            "cds_length": len(parsed.cds),
            "utr3_length": len(parsed.utr3),
            "poly_a_length": POLY_A_LENGTH,
            "num_codons": len(parsed.codons),
        },
        "structure_scores": structure_scores,
        "manufacturing_scores": manufacturing_scores,
        "stability_scores": stability_scores,
        "ribonn_scores": ribonn_scores,
        "summary": summary,
    }


def _traffic_light(value: float | None, green_range: tuple, amber_range: tuple) -> str:
    if value is None:
        return "GREY"
    g_lo, g_hi = green_range
    a_lo, a_hi = amber_range
    if g_lo <= value <= g_hi:
        return "GREEN"
    elif a_lo <= value <= a_hi:
        return "AMBER"
    return "RED"
