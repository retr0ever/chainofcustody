import numpy as np
from pymoo.core.problem import Problem

# Nucleotide encoding: 0=A, 1=C, 2=G, 3=T
NUCLEOTIDES = np.array(["A", "C", "G", "T"])
N_NUCLEOTIDES = len(NUCLEOTIDES)
N_OBJECTIVES = 3


class SequenceProblem(Problem):
    """Multi-objective sequence optimisation problem.

    Variables are nucleotide positions encoded as integers 0–3 (A/C/G/T).
    Objectives (all minimised by NSGA3):
      F[0] — fraction of A (encoding 0)
      F[1] — fraction of C (encoding 1)
      F[2] — fraction of T (encoding 3)
    """

    def __init__(self, seq_len: int = 100) -> None:
        super().__init__(
            n_var=seq_len,
            n_obj=N_OBJECTIVES,
            xl=0,
            xu=N_NUCLEOTIDES - 1,
            vtype=int,
        )

    def _evaluate(self, X: np.ndarray, out: dict, *args, **kwargs) -> None:
        n_pos = X.shape[1]
        # Fraction of each nucleotide per sequence (minimised by NSGA3).
        # Encoding: 0=A, 1=C, 2=G, 3=T
        frac_a = np.sum(X == 0, axis=1) / n_pos
        frac_c = np.sum(X == 1, axis=1) / n_pos
        frac_t = np.sum(X == 3, axis=1) / n_pos
        out["F"] = np.column_stack([frac_a, frac_c, frac_t])

    def decode(self, X: np.ndarray) -> list[str]:
        """Convert integer-encoded population rows to nucleotide strings."""
        return ["".join(NUCLEOTIDES[row]) for row in X]
