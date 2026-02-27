import numpy as np
from pymoo.core.problem import ElementwiseProblem

from chainofcustody.evaluation.report import score_sequence
from chainofcustody.evaluation.fitness import compute_fitness

# Nucleotide encoding: 0=A, 1=C, 2=G, 3=U
NUCLEOTIDES = np.array(["A", "C", "G", "U"])
N_NUCLEOTIDES = len(NUCLEOTIDES)

# One objective per fitness metric + overall
METRIC_NAMES = [
    "codon_quality",
    "gc_content",
    "mir122_detargeting",
    "utr5_accessibility",
    "manufacturability",
    "stability",
]
N_OBJECTIVES = len(METRIC_NAMES)


class SequenceProblem(ElementwiseProblem):
    """Multi-objective sequence optimisation problem.

    Variables are nucleotide positions encoded as integers 0–3 (A/C/G/U).
    Objectives: minimise (1 - metric_score) for each of the 6 evaluation metrics.
    Lower = better (pymoo minimises).

    Inherits from ElementwiseProblem so that individual solutions are dispatched
    one-at-a-time, enabling transparent parallelisation via an elementwise_runner
    (e.g. StarmapParallelization backed by a process or thread pool).
    """

    def __init__(self, seq_len: int = 100, **kwargs) -> None:
        super().__init__(
            n_var=seq_len,
            n_obj=N_OBJECTIVES,
            xl=0,
            xu=N_NUCLEOTIDES - 1,
            vtype=int,
            **kwargs,
        )

    def _evaluate(self, x: np.ndarray, out: dict, *args, **kwargs) -> None:
        """Evaluate a single solution vector ``x``."""
        seq = "".join(NUCLEOTIDES[x])
        try:
            report = score_sequence(seq)
            fitness = compute_fitness(report)
            out["F"] = np.array([1.0 - fitness["scores"][m]["value"] for m in METRIC_NAMES])
        except Exception:
            out["F"] = np.ones(N_OBJECTIVES)

    def decode(self, X: np.ndarray) -> list[str]:
        """Convert integer-encoded population rows to nucleotide strings."""
        return ["".join(NUCLEOTIDES[row]) for row in X]
