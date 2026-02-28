import numpy as np
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.core.callback import Callback
from pymoo.operators.crossover.ux import UniformCrossover
from pymoo.optimize import minimize
from pymoo.util.ref_dirs import get_reference_directions

from chainofcustody.evaluation.fitness import DEFAULT_WEIGHTS
from chainofcustody.optimization.operators import NucleotideMutation, NucleotideSampling
from chainofcustody.sequence import KOZAK
from chainofcustody.optimization.problem import METRIC_NAMES, N_OBJECTIVES, NUCLEOTIDES, SequenceProblem, assemble_mrna


class _ProgressCallback(Callback):
    """Advance a Rich progress bar by one step after each generation."""

    def __init__(self, progress, task, n_gen: int) -> None:
        super().__init__()
        self._progress = progress
        self._task = task
        self._n_gen = n_gen

    def notify(self, algorithm) -> None:
        gen = algorithm.n_gen
        self._progress.update(
            self._task,
            advance=1,
            description=f"Evolving 5'UTR  (gen {gen}/{self._n_gen})",
        )


def build_algorithm(
    pop_size: int = 128,
    mutation_rate: float = 0.05,
    initial_length: int | None = None,
    max_length_delta: int = 50,
) -> NSGA3:
    """Construct an NSGA3 instance for nucleotide sequence optimisation.

    Args:
        pop_size: Number of individuals in the population.
        mutation_rate: Per-position probability of a point mutation.
        initial_length: Seed the population around this 5'UTR length (None = uniform).
        max_length_delta: Maximum nt change to the length variable per mutation event.
    """
    # n_partitions=3 → 84 Das-Dennis reference points for 6 objectives,
    # which fits within the default pop_size=128 (NSGA-III requires pop_size ≥ n_ref_points).
    ref_dirs = get_reference_directions("das-dennis", N_OBJECTIVES, n_partitions=3)

    return NSGA3(
        ref_dirs=ref_dirs,
        pop_size=pop_size,
        sampling=NucleotideSampling(initial_length=initial_length),
        crossover=UniformCrossover(),
        mutation=NucleotideMutation(mutation_rate=mutation_rate, max_length_delta=max_length_delta),
        eliminate_duplicates=True,
    )


def _build_history(result, cds: str, utr3: str) -> list[dict]:
    """Extract per-generation population snapshots from a pymoo result."""
    records = []
    for gen_state in result.history:
        gen = gen_state.n_gen
        X = gen_state.pop.get("X")
        F = gen_state.pop.get("F")
        if X is None or F is None:
            continue
        for x_row, f_row in zip(X, F):
            utr5_len = int(x_row[0])
            utr5 = "".join(NUCLEOTIDES[x_row[1:utr5_len + 1]])
            seq = assemble_mrna(utr5, cds, utr3)
            scores = {m: round(1.0 - float(f_val), 4) for m, f_val in zip(METRIC_NAMES, f_row)}
            overall = round(sum(scores[m] * DEFAULT_WEIGHTS.get(m, 0) for m in METRIC_NAMES), 4)
            records.append({"generation": gen, "sequence": seq, **scores, "overall": overall})
    return records


def run(
    utr5_min: int = 20,
    utr5_max: int = 1000,
    cds: str = "",
    utr3: str = "",
    pop_size: int = 128,
    n_gen: int = 50,
    mutation_rate: float = 0.05,
    seed: int | None = None,
    verbose: bool = False,
    n_workers: int | None = None,
    progress=None,
    progress_task=None,
    target_cell_type: str = "megakaryocytes",
    initial_length: int | None = 200,
    max_length_delta: int = 50,
) -> tuple[np.ndarray, np.ndarray, list[dict]]:
    """Run NSGA3 on the sequence optimisation problem.

    The genetic algorithm evolves the 5'UTR region — both its sequence and its
    length (which varies between utr5_min and utr5_max). The CDS and 3'UTR are
    fixed and concatenated to the evolved 5'UTR before each evaluation.

    RiboNN inference is batched across the whole population each generation,
    keeping GPU utilisation high. The ``n_workers`` argument is accepted for
    API compatibility but ignored — multiprocessing is not used because the
    GPU-batched evaluator is already faster and CUDA cannot be forked safely.

    Args:
        utr5_min: Minimum allowed 5'UTR length.
        utr5_max: Maximum allowed 5'UTR length.
        cds: Fixed CDS sequence (RNA, uppercase).
        utr3: Fixed 3'UTR sequence (RNA, uppercase).
        pop_size: Population size.
        n_gen: Number of generations to evolve.
        mutation_rate: Per-position point-mutation probability.
        seed: Random seed for reproducibility.
        verbose: Print per-generation progress.
        n_workers: Ignored. Kept for API compatibility.
        progress: Optional Rich Progress instance for a live progress bar.
        progress_task: Task ID returned by ``progress.add_task``.
        initial_length: Seed population 5'UTR lengths around this value (None = uniform).
        max_length_delta: Maximum nt change to the length variable per mutation event.

    Returns:
        A tuple ``(X, F, history)`` where ``X`` is the integer-encoded
        Pareto-front population, ``F`` are the corresponding objective values,
        and ``history`` is a list of per-generation population records
        (full assembled sequences) suitable for CSV export.
    """
    algorithm = build_algorithm(
        pop_size=pop_size,
        mutation_rate=mutation_rate,
        initial_length=initial_length,
        max_length_delta=max_length_delta,
    )
    problem = SequenceProblem(utr5_min=utr5_min, utr5_max=utr5_max, cds=cds, utr3=utr3, target_cell_type=target_cell_type)

    from chainofcustody.progress import update_status  # noqa: PLC0415
    update_status("loading RiboNN models into GPU…")
    from chainofcustody.evaluation.ribonn import get_predictor  # noqa: PLC0415
    get_predictor()
    update_status("models ready")

    minimize_kwargs = dict(
        termination=("n_gen", n_gen),
        seed=seed,
        verbose=verbose,
        save_history=True,
    )

    if progress is not None:
        minimize_kwargs["callback"] = _ProgressCallback(progress, progress_task, n_gen)

    result = minimize(problem, algorithm, **minimize_kwargs)
    return result.X, result.F, _build_history(result, cds, utr3)
