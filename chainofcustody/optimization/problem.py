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

    Only the 5'UTR is evolved; the CDS and 3'UTR are fixed at construction time.
    Decision variables are nucleotide positions in the 5'UTR encoded as integers
    0–3 (A/C/G/U). The full sequence evaluated by the scorer is:

        5'UTR (evolved)  +  CDS (fixed)  +  3'UTR (fixed)

    Objectives: minimise (1 - metric_score) for each of the 6 evaluation metrics.
    Lower = better (pymoo minimises).

    Inherits from ElementwiseProblem so that individual solutions are dispatched
    one-at-a-time, enabling transparent parallelisation via an elementwise_runner
    (e.g. StarmapParallelization backed by a process or thread pool).
    """

    def __init__(self, utr5_len: int, cds: str, utr3: str, **kwargs) -> None:
        self.cds = cds
        self.utr3 = utr3
        super().__init__(
            n_var=utr5_len,
            n_obj=N_OBJECTIVES,
            xl=0,
            xu=N_NUCLEOTIDES - 1,
            vtype=int,
            **kwargs,
        )

    def _evaluate(self, x: np.ndarray, out: dict, *args, **kwargs) -> None:
        """Evaluate a single solution vector ``x`` (the 5'UTR encoding)."""
        utr5 = "".join(NUCLEOTIDES[x])
        full_seq = utr5 + self.cds + self.utr3
        try:
            report = score_sequence(
                full_seq,
                utr5_end=len(utr5),
                cds_end=len(utr5) + len(self.cds),
            )
            fitness = compute_fitness(report)
            out["F"] = np.array([1.0 - fitness["scores"][m]["value"] for m in METRIC_NAMES])
        except Exception:
            out["F"] = np.ones(N_OBJECTIVES)

    def decode(self, X: np.ndarray) -> list[str]:
        """Convert integer-encoded 5'UTR rows to full assembled sequences."""
        return ["".join(NUCLEOTIDES[row]) + self.cds + self.utr3 for row in X]
