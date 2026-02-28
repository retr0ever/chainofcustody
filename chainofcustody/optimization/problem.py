import logging
import os
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from pymoo.core.problem import Problem

from chainofcustody.sequence import KOZAK, mRNASequence
from chainofcustody.evaluation.scoring import score_parsed
from chainofcustody.evaluation.ribonn import score_ribonn_batch
from chainofcustody.evaluation.fitness import compute_fitness
from chainofcustody.progress import update_status, update_best_score

_CPU_WORKERS = os.cpu_count() or 1

logger = logging.getLogger(__name__)

# Nucleotide encoding: 0=A, 1=C, 2=G, 3=U
NUCLEOTIDES = np.array(["A", "C", "G", "U"])
N_NUCLEOTIDES = len(NUCLEOTIDES)

# One objective per fitness metric
METRIC_NAMES = [
    "utr5_accessibility",
    "manufacturability",
    "stability",
    "specificity",
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


class SequenceProblem(Problem):
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

    Objectives: minimise (1 - metric_score) for each of the 4 evaluation
    metrics. Lower = better (pymoo minimises).

    Inherits from ``Problem`` (vectorised) so that the whole population is
    evaluated in a single call. RiboNN inference is batched across the entire
    population, keeping GPU utilisation high. ViennaRNA folding runs in a
    ThreadPoolExecutor (it releases the GIL, so threads scale well).
    """

    def __init__(self, utr5_min: int, utr5_max: int, cds: str, utr3: str, target_cell_type: str = "megakaryocytes", **kwargs) -> None:
        if utr5_min < 0 or utr5_max < utr5_min:
            raise ValueError(
                f"utr5_min/utr5_max must satisfy 0 ≤ utr5_min ≤ utr5_max, "
                f"got {utr5_min}/{utr5_max}"
            )
        self.utr5_min = utr5_min
        self.utr5_max = utr5_max
        self.cds = cds
        self.utr3 = utr3
        self.target_cell_type = target_cell_type
        self._gen = 0  # incremented on each _evaluate call

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

    def _evaluate(self, X: np.ndarray, out: dict, *args, **kwargs) -> None:
        """Evaluate the entire population matrix ``X`` (shape: pop_size × n_var).

        Strategy:
        1. Batch RiboNN GPU inference across all sequences simultaneously.
        2. CPU-bound scoring (ViennaRNA folding, manufacturing, stability) is
           parallelised across all cores with a ThreadPoolExecutor (ViennaRNA
           releases the GIL so threads scale well).
        """
        self._gen += 1
        gen_tag = f"gen {self._gen}"
        n = len(X)
        F = np.ones((n, N_OBJECTIVES))

        # Decode all chromosomes into mRNASequence objects
        parsed_list: list[mRNASequence] = []
        for row in X:
            utr5_len = int(row[0])
            utr5 = "".join(NUCLEOTIDES[row[1:utr5_len + 1]])
            parsed_list.append(
                mRNASequence(utr5=utr5 + KOZAK, cds=self.cds, utr3=self.utr3)
            )

        # --- GPU: RiboNN batch inference ---
        update_status(f"{gen_tag}  RiboNN GPU inference ({n} seqs)")
        try:
            ribonn_results = score_ribonn_batch(parsed_list, target_cell_type=self.target_cell_type)
        except Exception as exc:
            logger.warning("Batch RiboNN scoring failed: %s", exc)
            ribonn_results = [None] * n

        # --- CPU: parallel ViennaRNA folding + manufacturing + stability ---
        update_status(f"{gen_tag}  CPU scoring ({n} seqs, {_CPU_WORKERS} threads)")

        def _score_one(args: tuple[int, mRNASequence, dict | None]) -> tuple[int, np.ndarray]:
            idx, parsed, ribonn_scores = args
            try:
                report = score_parsed(parsed, _ribonn_scores=ribonn_scores, _fast_fold=True, target_cell_type=self.target_cell_type)
                fitness = compute_fitness(report)
                f_row = np.array([1.0 - fitness["scores"][m]["value"] for m in METRIC_NAMES])
            except Exception as exc:
                logger.warning(
                    "Scoring failed for sequence %r…: %s", str(parsed)[:30], exc
                )
                f_row = np.ones(N_OBJECTIVES)
            return idx, f_row

        work = list(zip(range(n), parsed_list, ribonn_results))
        with ThreadPoolExecutor(max_workers=_CPU_WORKERS) as pool:
            for idx, f_row in pool.map(_score_one, work):
                F[idx] = f_row

        update_status(f"{gen_tag}  done")
        out["F"] = F

        # Broadcast the best weighted overall score in this generation
        from chainofcustody.evaluation.fitness import DEFAULT_WEIGHTS  # noqa: PLC0415
        weights = np.array([DEFAULT_WEIGHTS.get(m, 0) for m in METRIC_NAMES])
        overall_scores = 1.0 - F @ weights  # shape (n,); higher = better
        update_best_score(float(overall_scores.max()))

    def decode(self, X: np.ndarray) -> list[str]:
        """Convert integer-encoded rows to full assembled sequences."""
        result = []
        for row in X:
            utr5_len = int(row[0])
            utr5 = "".join(NUCLEOTIDES[row[1:utr5_len + 1]])
            result.append(assemble_mrna(utr5, self.cds, self.utr3))
        return result
