"""Orchestrate the full 6-metric scoring pipeline over a parsed mRNA sequence."""

from chainofcustody.sequence import mRNASequence
from chainofcustody.evaluation.codons import score_codons
from chainofcustody.evaluation.mirna import score_mirna
from chainofcustody.evaluation.structure import score_structure
from chainofcustody.evaluation.manufacturing import score_manufacturing
from chainofcustody.evaluation.stability import score_stability


def score_parsed(
    parsed: mRNASequence,
    target: str | None = None,
) -> dict:
    """Run all 6 evaluation metrics on an already-parsed mRNA sequence.

    This is the canonical scoring entry point when the 5'UTR / CDS / 3'UTR
    boundaries are already known (e.g. inside the optimiser). It avoids the
    assemble-then-re-parse roundtrip required by ``score_sequence``.

    Args:
        parsed: An ``mRNASequence`` whose regions are correctly delimited.
        target: Optional target cell-type name for selectivity scoring.

    Returns:
        Full report dict with keys: ``sequence_info``, ``codon_scores``,
        ``mirna_scores``, ``structure_scores``, ``manufacturing_scores``,
        ``stability_scores``, ``summary``.
    """
    codon_scores = score_codons(parsed, target_cell_type=target)
    mirna_scores = score_mirna(parsed)

    mir122_positions = (
        mirna_scores.get("detargeting", {})
        .get("miR-122-5p", {})
        .get("positions", [])
    )
    structure_scores = score_structure(parsed, mirna_site_positions=mir122_positions)
    manufacturing_scores = score_manufacturing(parsed)
    stability_scores = score_stability(parsed)

    cai = codon_scores.get("cai", 0)
    gc = codon_scores.get("gc_content", {}).get("cds", 0)
    mir122_count = (
        mirna_scores.get("detargeting", {}).get("miR-122-5p", {}).get("utr3_sites", 0)
    )
    mfg_violations = manufacturing_scores.get("total_violations", 0)

    summary = {
        "codon_quality": _traffic_light(cai, (0.8, 1.0), (0.6, 1.0)),
        "gc_content": _traffic_light(gc, (40, 60), (30, 70)),
        "mir122_detargeting": _traffic_light(mir122_count, (3, 100), (1, 100)),
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
        "codon_scores": codon_scores,
        "mirna_scores": mirna_scores,
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
