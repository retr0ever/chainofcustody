import numpy as np
from pymoo.core.mutation import Mutation
from pymoo.core.sampling import Sampling

from chainofcustody.optimization.problem import N_NUCLEOTIDES


class NucleotideSampling(Sampling):
    """Generate a random initial population of nucleotide sequences."""

    def _do(self, problem, n_samples: int, **kwargs) -> np.ndarray:
        return np.random.randint(0, N_NUCLEOTIDES, size=(n_samples, problem.n_var))


class NucleotideMutation(Mutation):
    """Per-position point mutation on nucleotide sequences.

    Each position is independently replaced with a uniformly drawn nucleotide
    (possibly the same one) with probability `mutation_rate`.
    """

    def __init__(self, mutation_rate: float = 0.01) -> None:
        super().__init__()
        if not 0.0 <= mutation_rate <= 1.0:
            raise ValueError(f"mutation_rate must be in [0, 1], got {mutation_rate}")
        self.mutation_rate = mutation_rate

    def _do(self, problem, X: np.ndarray, **kwargs) -> np.ndarray:
        mutated = X.copy()
        mask = np.random.random(X.shape) < self.mutation_rate
        mutated[mask] = np.random.randint(0, N_NUCLEOTIDES, size=mask.sum())
        return mutated
