import numpy as np
from pymoo.core.problem import Problem

# Nucleotide encoding: 0=A, 1=C, 2=G, 3=U
NUCLEOTIDES = np.array(["A", "C", "G", "U"])
N_NUCLEOTIDES = len(NUCLEOTIDES)
N_OBJECTIVES = 3


class SequenceProblem(Problem):
    """Multi-objective sequence optimisation problem.

    Variables are nucleotide positions encoded as integers 0–3 (A/C/G/U).
    Objectives (all minimised by NSGA3):
      F[0] — fraction of A (encoding 0)
      F[1] — fraction of C (encoding 1)
      F[2] — fraction of U (encoding 3)
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
        sequences = self.decode(X)
        frac_a = np.array([s.count("A") / len(s) for s in sequences])
        frac_c = np.array([s.count("C") / len(s) for s in sequences])
        frac_u = np.array([s.count("U") / len(s) for s in sequences])
        out["F"] = np.column_stack([frac_a, frac_c, frac_u])

    def decode(self, X: np.ndarray) -> list[str]:
        """Convert integer-encoded population rows to nucleotide strings."""
        return ["".join(NUCLEOTIDES[row]) for row in X]
