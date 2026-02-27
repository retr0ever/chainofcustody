import numpy as np
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.operators.crossover.ux import UniformCrossover
from pymoo.optimize import minimize
from pymoo.util.ref_dirs import get_reference_directions

from chainofcustody.optimization.operators import NucleotideMutation, NucleotideSampling
from chainofcustody.optimization.problem import N_OBJECTIVES, SequenceProblem


def build_algorithm(
    pop_size: int = 100,
    mutation_rate: float = 0.01,
) -> NSGA3:
    """Construct an NSGA3 instance for nucleotide sequence optimisation.

    Args:
        pop_size: Number of individuals in the population.
        mutation_rate: Per-position probability of a point mutation.
    """
    ref_dirs = get_reference_directions("das-dennis", N_OBJECTIVES, n_partitions=12)

    return NSGA3(
        ref_dirs=ref_dirs,
        pop_size=pop_size,
        sampling=NucleotideSampling(),
        crossover=UniformCrossover(),
        mutation=NucleotideMutation(mutation_rate=mutation_rate),
        eliminate_duplicates=True,
    )


def run(
    seq_len: int = 100,
    pop_size: int = 100,
    n_gen: int = 50,
    mutation_rate: float = 0.01,
    seed: int | None = None,
    verbose: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Run NSGA3 on the sequence optimisation problem.

    Args:
        seq_len: Length of each candidate nucleotide sequence.
        pop_size: Population size.
        n_gen: Number of generations to evolve.
        mutation_rate: Per-position point-mutation probability.
        seed: Random seed for reproducibility.
        verbose: Print per-generation progress.

    Returns:
        A tuple ``(X, F)`` where ``X`` is the integer-encoded Pareto-front
        population and ``F`` are the corresponding objective values.
    """
    problem = SequenceProblem(seq_len=seq_len)
    algorithm = build_algorithm(pop_size=pop_size, mutation_rate=mutation_rate)

    result = minimize(
        problem,
        algorithm,
        termination=("n_gen", n_gen),
        seed=seed,
        verbose=verbose,
    )

    return result.X, result.F
