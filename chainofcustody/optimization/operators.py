import numpy as np
from pymoo.core.mutation import Mutation
from pymoo.core.sampling import Sampling

from chainofcustody.optimization.problem import N_NUCLEOTIDES, NUCLEOTIDES

_NUCLEOTIDE_INDEX = {nt: i for i, nt in enumerate(NUCLEOTIDES)}

# Shared UTR sequence used as seed for 5'UTR evolution and as the fixed 3'UTR.
UTR_SEED = "GAGUAGUCCCUUCGCAAGCCCUCAUUUCACCAGGCCCCCGGCUUGGGGCGCCUUCCUUCCCC"


def _encode(seq: str) -> np.ndarray:
    """Encode a nucleotide string to an integer array (A=0, C=1, G=2, U=3).

    DNA thymine (T) is mapped to U so external DNA sequences are accepted.
    """
    return np.array([_NUCLEOTIDE_INDEX[nt.upper().replace("T", "U")] for nt in seq])


class NucleotideSampling(Sampling):
    """Seed the initial population from UTR_SEED with random mutations.

    Chromosome layout mirrors SequenceProblem:
        column 0          — 5'UTR length (seeded to len(UTR_SEED), clamped to [utr5_min, utr5_max])
        columns 1..utr5_max — nucleotide positions (seeded from UTR_SEED, tiled/trimmed as needed)

    Each individual is independently mutated at every nucleotide position with
    probability ``seed_mutation_rate`` so the optimizer starts with diversity
    around the known good sequence.
    """

    def __init__(self, seed_mutation_rate: float = 0.1) -> None:
        super().__init__()
        self.seed_mutation_rate = seed_mutation_rate

    def _do(self, problem, n_samples: int, **kwargs) -> np.ndarray:
        utr5_min = int(problem.xl[0])
        utr5_max = int(problem.xu[0])

        # Seed length: len(UTR_SEED) clamped to the allowed range
        seed_len = int(np.clip(len(UTR_SEED), utr5_min, utr5_max))

        # Seed nucleotides tiled/trimmed to fill all utr5_max positions
        seed_nts = _encode(UTR_SEED)
        if len(seed_nts) < utr5_max:
            repeats = -(-utr5_max // len(seed_nts))  # ceiling division
            seed_nts = np.tile(seed_nts, repeats)
        seed_nts = seed_nts[:utr5_max]

        # Build chromosome: [length, nt_1, nt_2, ..., nt_utr5_max]
        seed_row = np.concatenate([[seed_len], seed_nts])
        population = np.tile(seed_row, (n_samples, 1))

        # Mutate only the nucleotide columns (skip column 0 = length)
        nt_cols = population[:, 1:]
        mask = np.random.random(nt_cols.shape) < self.seed_mutation_rate
        nt_cols[mask] = np.random.randint(0, N_NUCLEOTIDES, size=mask.sum())
        population[:, 1:] = nt_cols

        return population


class NucleotideMutation(Mutation):
    """Per-position point mutation on nucleotide sequences.

    Each position is independently replaced with a uniformly drawn value within
    that variable's bounds (from problem.xl / problem.xu) with probability
    `mutation_rate`. This correctly handles both the length variable (column 0,
    bounded [utr5_min, utr5_max]) and the nucleotide variables (columns 1+,
    bounded [0, 3]).
    """

    def __init__(self, mutation_rate: float = 0.01) -> None:
        super().__init__()
        if not 0.0 <= mutation_rate <= 1.0:
            raise ValueError(f"mutation_rate must be in [0, 1], got {mutation_rate}")
        self.mutation_rate = mutation_rate

    def _do(self, problem, X: np.ndarray, **kwargs) -> np.ndarray:
        mutated = X.copy()
        xl = problem.xl.astype(int)
        xu = problem.xu.astype(int)
        for j in range(X.shape[1]):
            col_mask = np.random.random(X.shape[0]) < self.mutation_rate
            if col_mask.any():
                mutated[col_mask, j] = np.random.randint(xl[j], xu[j] + 1, size=col_mask.sum())
        return mutated
