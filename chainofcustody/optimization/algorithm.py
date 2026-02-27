import os
from multiprocessing.pool import Pool

import numpy as np
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.core.callback import Callback
from pymoo.operators.crossover.ux import UniformCrossover
from pymoo.optimize import minimize
from pymoo.parallelization.starmap import StarmapParallelization
from pymoo.util.ref_dirs import get_reference_directions

from chainofcustody.evaluation.fitness import DEFAULT_WEIGHTS
from chainofcustody.optimization.operators import NucleotideMutation, NucleotideSampling
from chainofcustody.optimization.problem import METRIC_NAMES, N_OBJECTIVES, NUCLEOTIDES, SequenceProblem

_DEFAULT_WORKERS = os.cpu_count() or 1


class _ProgressCallback(Callback):
    """Advance a Rich progress bar by one step after each generation."""

    def __init__(self, progress, task) -> None:
        super().__init__()
        self._progress = progress
        self._task = task

    def notify(self, algorithm) -> None:
        self._progress.advance(self._task)


def build_algorithm(
    pop_size: int = 128,
    mutation_rate: float = 0.01,
) -> NSGA3:
    """Construct an NSGA3 instance for nucleotide sequence optimisation.

    Args:
        pop_size: Number of individuals in the population.
        mutation_rate: Per-position probability of a point mutation.
    """
    ref_dirs = get_reference_directions("das-dennis", N_OBJECTIVES, n_partitions=4)

    return NSGA3(
        ref_dirs=ref_dirs,
        pop_size=pop_size,
        sampling=NucleotideSampling(),
        crossover=UniformCrossover(),
        mutation=NucleotideMutation(mutation_rate=mutation_rate),
        eliminate_duplicates=True,
    )


def _build_history(result) -> list[dict]:
    """Extract per-generation population snapshots from a pymoo result."""
    records = []
    for gen_state in result.history:
        gen = gen_state.n_gen
        X = gen_state.pop.get("X")
        F = gen_state.pop.get("F")
        if X is None or F is None:
            continue
        for x_row, f_row in zip(X, F):
            seq = "".join(NUCLEOTIDES[x_row])
            scores = {m: round(1.0 - float(f_val), 4) for m, f_val in zip(METRIC_NAMES, f_row)}
            overall = round(sum(scores[m] * DEFAULT_WEIGHTS.get(m, 0) for m in METRIC_NAMES), 4)
            records.append({"generation": gen, "sequence": seq, **scores, "overall": overall})
    return records


def run(
    seq_len: int = 100,
    pop_size: int = 128,
    n_gen: int = 50,
    mutation_rate: float = 0.01,
    seed: int | None = None,
    verbose: bool = False,
    n_workers: int | None = None,
    progress=None,
    progress_task=None,
) -> tuple[np.ndarray, np.ndarray, list[dict]]:
    """Run NSGA3 on the sequence optimisation problem.

    Args:
        seq_len: Length of each candidate nucleotide sequence.
        pop_size: Population size.
        n_gen: Number of generations to evolve.
        mutation_rate: Per-position point-mutation probability.
        seed: Random seed for reproducibility.
        verbose: Print per-generation progress.
        n_workers: Number of worker processes for parallel fitness evaluation.
            Defaults to the number of logical CPU cores. Pass ``1`` to run
            single-threaded (no pool overhead).
        progress: Optional Rich Progress instance for a live progress bar.
        progress_task: Task ID returned by ``progress.add_task``.

    Returns:
        A tuple ``(X, F, history)`` where ``X`` is the integer-encoded
        Pareto-front population, ``F`` are the corresponding objective values,
        and ``history`` is a list of per-generation population records suitable
        for CSV export.
    """
    workers = n_workers if n_workers is not None else _DEFAULT_WORKERS

    algorithm = build_algorithm(pop_size=pop_size, mutation_rate=mutation_rate)

    minimize_kwargs = dict(
        termination=("n_gen", n_gen),
        seed=seed,
        verbose=verbose,
        save_history=True,
    )

    if progress is not None:
        minimize_kwargs["callback"] = _ProgressCallback(progress, progress_task)

    if workers == 1:
        problem = SequenceProblem(seq_len=seq_len)
        result = minimize(problem, algorithm, **minimize_kwargs)
    else:
        with Pool(workers) as pool:
            runner = StarmapParallelization(pool.starmap)
            problem = SequenceProblem(seq_len=seq_len, elementwise_runner=runner)
            result = minimize(problem, algorithm, **minimize_kwargs)

    return result.X, result.F, _build_history(result)
