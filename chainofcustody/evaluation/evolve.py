"""Evolutionary optimisation loop for mRNA sequences."""

import random
from itertools import combinations
from typing import Callable

from .report import score_sequence
from .fitness import compute_fitness
from .mutations import synonymous_swap, crossover, insert_mir122_sites
from .parser import parse_sequence


def _score_population(
    population: list[str],
    generation: int,
    target: str | None = None,
    weights: dict[str, float] | None = None,
) -> list[dict]:
    """Score all candidates and return sorted results (best first)."""
    results = []
    for i, seq in enumerate(population):
        try:
            report = score_sequence(seq, target=target)
            fitness = compute_fitness(report, weights=weights)
            results.append({
                "label": f"gen{generation}_c{i}",
                "seq": seq,
                "report": report,
                "fitness": fitness,
            })
        except Exception:
            continue

    results.sort(key=lambda r: r["fitness"]["overall"], reverse=True)
    return results


def _compute_metric_ranges(results: list[dict]) -> dict:
    """Compute min/max/mean for each metric across a population."""
    if not results:
        return {}

    metrics = list(results[0]["fitness"]["scores"].keys())
    ranges = {}

    for metric in metrics:
        values = [r["fitness"]["scores"][metric]["value"] for r in results]
        ranges[metric] = {
            "min": round(min(values), 4),
            "max": round(max(values), 4),
            "mean": round(sum(values) / len(values), 4),
        }

    return ranges


def _breed(parents: list[dict], children_per_pair: int, mutation_rate: float) -> list[str]:
    """Generate children from parent pairs via crossover + mutation."""
    children = []
    parent_seqs = [p["seq"] for p in parents]

    # All unique pairs of parents
    pairs = list(combinations(range(len(parent_seqs)), 2))
    if not pairs:
        # Single parent — just mutate
        pairs = [(0, 0)]

    for a_idx, b_idx in pairs:
        for _ in range(children_per_pair):
            if a_idx == b_idx:
                # Self-pair: just mutate
                child = parent_seqs[a_idx]
            else:
                child, _ = crossover(parent_seqs[a_idx], parent_seqs[b_idx])

            # Mutate
            try:
                parsed = parse_sequence(child)
                child = synonymous_swap(child, parsed, rate=mutation_rate)

                # If no miR-122 sites, insert them
                parsed = parse_sequence(child)
                mir122_sites = 0
                report = parents[a_idx]["report"]
                mir122_data = report["mirna_scores"].get("detargeting", {}).get("miR-122-5p", {})
                mir122_sites = mir122_data.get("utr3_sites", 0)
                if mir122_sites == 0:
                    child = insert_mir122_sites(child, parsed, n_sites=3)

                children.append(child)
            except Exception:
                continue

    return children


def evolve(
    population: list[str],
    generations: int = 20,
    top_k: int = 3,
    children_per_pair: int = 2,
    mutation_rate: float = 0.05,
    target: str | None = None,
    weights: dict[str, float] | None = None,
    on_generation: Callable | None = None,
) -> dict:
    """
    Run evolutionary optimisation on mRNA sequences.

    Args:
        population: Initial candidate sequences (raw DNA strings).
        generations: Number of iterations.
        top_k: How many top candidates survive as parents.
        children_per_pair: Children produced per parent pair.
        mutation_rate: Fraction of codons to mutate per child.
        target: Target cell type column name.
        weights: Custom metric weights.
        on_generation: Callback(generation_num, generation_results) for progress.

    Returns:
        {
            "best": best candidate result,
            "history": [generation stats with metric_ranges],
            "final_population": [all final candidates, ranked],
        }
    """
    history = []
    best_overall = 0.0
    stale_count = 0
    current_pop = list(population)

    for gen in range(generations):
        # Score everyone
        results = _score_population(current_pop, gen, target=target, weights=weights)

        if not results:
            break

        gen_best = results[0]["fitness"]["overall"]
        gen_avg = sum(r["fitness"]["overall"] for r in results) / len(results)
        metric_ranges = _compute_metric_ranges(results)

        gen_stats = {
            "generation": gen,
            "population_size": len(results),
            "best_fitness": round(gen_best, 4),
            "avg_fitness": round(gen_avg, 4),
            "metric_ranges": metric_ranges,
        }
        history.append(gen_stats)

        if on_generation:
            on_generation(gen, results)

        # Early stopping
        if gen_best > best_overall:
            best_overall = gen_best
            stale_count = 0
        else:
            stale_count += 1

        if stale_count >= 5:
            break

        # Select parents (elitism)
        parents = results[:top_k]

        # Breed children
        children = _breed(parents, children_per_pair, mutation_rate)

        # New population = parents + children
        current_pop = [p["seq"] for p in parents] + children

    # Final scoring
    final_results = _score_population(current_pop, generations, target=target, weights=weights)

    return {
        "best": final_results[0] if final_results else None,
        "history": history,
        "final_population": final_results,
    }
