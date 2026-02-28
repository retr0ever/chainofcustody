"""Orchestrate the full 3-metric scoring pipeline over a parsed mRNA sequence."""

from chainofcustody.sequence import mRNASequence
from chainofcustody.evaluation.structure import score_structure
from chainofcustody.evaluation.manufacturing import score_manufacturing
from chainofcustody.evaluation.stability import score_stability


def score_parsed(
    parsed: mRNASequence,
    target: str | None = None,
) -> dict:
    """Run all 3 evaluation metrics on an already-parsed mRNA sequence.

    Args:
        parsed: An ``mRNASequence`` whose regions are correctly delimited.
        target: Unused — reserved for future per-cell-type scoring.

    Returns:
        Full report dict with keys: ``sequence_info``, ``structure_scores``,
        ``manufacturing_scores``, ``stability_scores``, ``summary``.
    """
    structure_scores = score_structure(parsed)
    manufacturing_scores = score_manufacturing(parsed)
    stability_scores = score_stability(parsed)

    mfg_violations = manufacturing_scores.get("total_violations", 0)

    summary = {
        "utr5_accessibility": structure_scores.get("utr5_accessibility", {}).get("status", "GREY"),
        "manufacturability": _traffic_light(-mfg_violations, (-3, 0), (-999, 0)),
        "stability": stability_scores.get("status", "GREY"),
    }

    return {
        "sequence_info": {
            "total_length": len(parsed),
            "utr5_length": len(parsed.utr5),
            "cds_length": len(parsed.cds),
            "utr3_length": len(parsed.utr3),
            "num_codons": len(parsed.codons),
        },
        "structure_scores": structure_scores,
        "manufacturing_scores": manufacturing_scores,
        "stability_scores": stability_scores,
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
