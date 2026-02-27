import numpy as np
from pymoo.core.mutation import Mutation
from pymoo.core.sampling import Sampling

from chainofcustody.optimization.problem import N_NUCLEOTIDES, NUCLEOTIDES

_NUCLEOTIDE_INDEX = {nt: i for i, nt in enumerate(NUCLEOTIDES)}

SEED_SEQUENCE = "GAGTAGTCCCTTCGCAAGCCCTCATTTCACCAGGCCCCCGGCTTGGGGCGCCTTCCTTCCCC"


def _encode(seq: str) -> np.ndarray:
    """Encode a nucleotide string to an integer array (A=0, C=1, G=2, U=3).

    DNA thymine (T) is mapped to U so DNA seed sequences are accepted.
    """
    return np.array([_NUCLEOTIDE_INDEX[nt.upper().replace("T", "U")] for nt in seq])


class NucleotideSampling(Sampling):
    """Seed the initial population from SEED_SEQUENCE with random mutations.

    The seed is tiled/trimmed to match the problem's sequence length, then each
    individual in the population is independently mutated at every position with
    probability ``seed_mutation_rate`` so the optimizer starts with diversity
    around the known good sequence rather than a purely random population.
    """

    def __init__(self, seed_mutation_rate: float = 0.1) -> None:
        super().__init__()
        self.seed_mutation_rate = seed_mutation_rate

    def _do(self, problem, n_samples: int, **kwargs) -> np.ndarray:
        seq_len = problem.n_var
        seed = _encode(SEED_SEQUENCE)

        # Tile or trim the seed to the required length
        if len(seed) < seq_len:
            repeats = -(-seq_len // len(seed))  # ceiling division
            seed = np.tile(seed, repeats)
        seed = seed[:seq_len]

        population = np.tile(seed, (n_samples, 1))
        mask = np.random.random(population.shape) < self.seed_mutation_rate
        population[mask] = np.random.randint(0, N_NUCLEOTIDES, size=mask.sum())
        return population


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
