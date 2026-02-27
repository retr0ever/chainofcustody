"""mRNA sequence evaluation pipeline."""

from .report import score_sequence, report_to_json
from .fitness import compute_fitness
from .evolve import evolve


def evaluate_candidate(
    seq: str,
    target: str | None = None,
    weights: dict[str, float] | None = None,
    label: str | None = None,
) -> dict:
    """Score a single candidate sequence. Returns report + fitness + suggestions."""
    report = score_sequence(seq, target=target)
    fitness = compute_fitness(report, weights=weights)
    return {"label": label or "candidate", "report": report, "fitness": fitness, "suggestions": fitness["suggestions"]}


def evaluate_batch(
    sequences: list[str] | list[dict],
    target: str | None = None,
    weights: dict[str, float] | None = None,
) -> list[dict]:
    """
    Score and rank multiple candidate sequences.

    Args:
        sequences: List of raw sequences, or list of {"seq": str, "label": str} dicts.
        target: Target cell type column name.
        weights: Custom metric weights.

    Returns: List of results sorted by fitness score (best first).
    """
    results = []
    for i, entry in enumerate(sequences):
        if isinstance(entry, dict):
            seq = entry["seq"]
            label = entry.get("label", f"candidate_{i + 1}")
        else:
            seq = entry
            label = f"candidate_{i + 1}"

        result = evaluate_candidate(seq, target=target, weights=weights, label=label)
        results.append(result)

    results.sort(key=lambda r: r["fitness"]["overall"], reverse=True)
    return results


__all__ = ["score_sequence", "compute_fitness", "evaluate_candidate", "evaluate_batch", "report_to_json", "evolve"]
