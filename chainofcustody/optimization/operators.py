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
    """Initialise the population with uniformly random 5'UTR sequences.

    Chromosome layout mirrors SequenceProblem:
        column 0            — 5'UTR length, centred on *initial_length* when provided
        columns 1..utr5_max — nucleotide positions, each drawn uniformly from {A, C, G, U}

    When *initial_length* is set, lengths are drawn from a Gaussian centred at that
    value (σ = 10% of initial_length) and clamped to [utr5_min, utr5_max].  This
    seeds the population at a sensible starting length while preserving diversity.
    Without *initial_length*, lengths are drawn uniformly across the full range.
    """

    def __init__(self, initial_length: int | None = None) -> None:
        super().__init__()
        self.initial_length = initial_length

    def _do(self, problem, n_samples: int, **kwargs) -> np.ndarray:
        utr5_min = int(problem.xl[0])
        utr5_max = int(problem.xu[0])
        n_var = utr5_max + 1

        population = np.zeros((n_samples, n_var), dtype=int)

        if self.initial_length is not None:
            init_len = int(np.clip(self.initial_length, utr5_min, utr5_max))
            sigma = max(1, int(init_len * 0.10))
            lengths = np.round(np.random.normal(init_len, sigma, n_samples)).astype(int)
            population[:, 0] = np.clip(lengths, utr5_min, utr5_max)
        else:
            population[:, 0] = np.random.randint(utr5_min, utr5_max + 1, size=n_samples)

        population[:, 1:] = np.random.randint(0, N_NUCLEOTIDES, size=(n_samples, utr5_max))
        return population


class NucleotideMutation(Mutation):
    """Per-position point mutation on nucleotide sequences.

    Nucleotide columns (1+) are mutated with *mutation_rate* probability by
    replacement with a uniformly drawn nucleotide.

    The length column (0) is mutated with the same probability, but using a
    bounded random-walk step in [-max_length_delta, +max_length_delta] rather
    than a jump to a completely random value.  This prevents disruptive length
    changes while still allowing the population to explore the full length range
    over many generations.
    """

    def __init__(self, mutation_rate: float = 0.01, max_length_delta: int = 50) -> None:
        super().__init__()
        if not 0.0 <= mutation_rate <= 1.0:
            raise ValueError(f"mutation_rate must be in [0, 1], got {mutation_rate}")
        self.mutation_rate = mutation_rate
        self.max_length_delta = max_length_delta

    def _do(self, problem, X: np.ndarray, **kwargs) -> np.ndarray:
        mutated = X.copy()
        xl = problem.xl.astype(int)
        xu = problem.xu.astype(int)

        # ── Nucleotide mutations (columns 1+): random replacement ──────────────
        nt_mask = np.random.random(X[:, 1:].shape) < self.mutation_rate
        if nt_mask.any():
            lo = xl[1:][np.newaxis, :]
            hi = xu[1:][np.newaxis, :]
            span = (hi - lo + 1).astype(np.float64)
            replacements = lo + (np.random.random(X[:, 1:].shape) * span).astype(int)
            mutated[:, 1:] = np.where(nt_mask, replacements, mutated[:, 1:])

        # ── Length mutation (column 0): bounded random walk ─────────────────────
        len_mask = np.random.random(X.shape[0]) < self.mutation_rate
        if len_mask.any():
            delta = np.random.randint(
                -self.max_length_delta, self.max_length_delta + 1, size=X.shape[0]
            )
            new_lengths = np.clip(mutated[:, 0] + delta, xl[0], xu[0])
            mutated[:, 0] = np.where(len_mask, new_lengths, mutated[:, 0])

        return mutated
