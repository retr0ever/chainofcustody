import logging
import numpy as np
from pymoo.core.problem import ElementwiseProblem

from chainofcustody.sequence import KOZAK, mRNASequence
from chainofcustody.evaluation.scoring import score_parsed
from chainofcustody.evaluation.fitness import compute_fitness

logger = logging.getLogger(__name__)

# Nucleotide encoding: 0=A, 1=C, 2=G, 3=U
NUCLEOTIDES = np.array(["A", "C", "G", "U"])
N_NUCLEOTIDES = len(NUCLEOTIDES)

# One objective per fitness metric + overall
METRIC_NAMES = [
    "utr5_accessibility",
    "manufacturability",
    "stability",
]
N_OBJECTIVES = len(METRIC_NAMES)


def assemble_mrna(utr5: str, cds: str, utr3: str) -> str:
    """Assemble a full mRNA sequence from its three regions.

    Inserts the Kozak consensus (``GCCACC``) between the 5'UTR and CDS.

    Args:
        utr5: 5'UTR sequence (RNA, variable length).
        cds:  Coding sequence (RNA, start codon through stop codon).
        utr3: 3'UTR sequence (RNA).

    Returns:
        Full mRNA: ``5'UTR + KOZAK + CDS + 3'UTR``.
    """
    return utr5 + KOZAK + cds + utr3


class SequenceProblem(ElementwiseProblem):
    """Multi-objective sequence optimisation problem.

    Only the 5'UTR is evolved; its length is also a decision variable. The CDS
    and 3'UTR are fixed at construction time. The full evaluated sequence is:

        5'UTR (evolved, variable length)  +  KOZAK  +  CDS (fixed)  +  3'UTR (fixed)

    The Kozak consensus (``GCCACC``) is always inserted between the evolved
    5'UTR and the CDS start codon. It is not part of the evolved chromosome.

    Chromosome layout (n_var = utr5_max + 1):
        x[0]          — 5'UTR length  ∈ [utr5_min, utr5_max]
        x[1 : x[0]+1] — active 5'UTR nucleotides (0=A, 1=C, 2=G, 3=U)
        x[x[0]+1 :]   — inactive padding (ignored during evaluation)

    Objectives: minimise (1 - metric_score) for each of the 6 evaluation
    metrics. Lower = better (pymoo minimises).

    Inherits from ElementwiseProblem so that individual solutions are dispatched
    one-at-a-time, enabling transparent parallelisation via an elementwise_runner
    (e.g. StarmapParallelization backed by a process or thread pool).
    """

    def __init__(self, utr5_min: int, utr5_max: int, cds: str, utr3: str, **kwargs) -> None:
        if utr5_min < 0 or utr5_max < utr5_min:
            raise ValueError(f"utr5_min/utr5_max must satisfy 0 ≤ utr5_min ≤ utr5_max, got {utr5_min}/{utr5_max}")
        self.utr5_min = utr5_min
        self.utr5_max = utr5_max
        self.cds = cds
        self.utr3 = utr3

        # Variable 0 encodes length; variables 1..utr5_max encode nucleotides.
        xl = np.array([utr5_min] + [0] * utr5_max)
        xu = np.array([utr5_max] + [N_NUCLEOTIDES - 1] * utr5_max)

        super().__init__(
            n_var=utr5_max + 1,
            n_obj=N_OBJECTIVES,
            xl=xl,
            xu=xu,
            vtype=int,
            **kwargs,
        )

    def _evaluate(self, x: np.ndarray, out: dict, *args, **kwargs) -> None:
        """Evaluate a single solution vector ``x``."""
        utr5_len = int(x[0])
        utr5 = "".join(NUCLEOTIDES[x[1:utr5_len + 1]])
        parsed = mRNASequence(
            utr5=utr5 + KOZAK,
            cds=self.cds,
            utr3=self.utr3,
        )
        try:
            report = score_parsed(parsed)
            fitness = compute_fitness(report)
            out["F"] = np.array([1.0 - fitness["scores"][m]["value"] for m in METRIC_NAMES])
        except Exception as exc:
            logger.warning("Scoring failed for sequence %r…: %s", str(parsed)[:30], exc)
            out["F"] = np.ones(N_OBJECTIVES)

    def decode(self, X: np.ndarray) -> list[str]:
        """Convert integer-encoded rows to full assembled sequences."""
        result = []
        for row in X:
            utr5_len = int(row[0])
            utr5 = "".join(NUCLEOTIDES[row[1:utr5_len + 1]])
            result.append(assemble_mrna(utr5, self.cds, self.utr3))
        return result
