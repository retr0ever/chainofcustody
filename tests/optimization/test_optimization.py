import numpy as np
import pytest

from chainofcustody.optimization import (
    KOZAK,
    METRIC_NAMES,
    NucleotideMutation,
    NucleotideSampling,
    SequenceProblem,
    build_algorithm,
    run,
)
from chainofcustody.three_prime.generate_utr3 import generate_mrna_sponge_utr

N_METRICS = len(METRIC_NAMES)

# Minimal fixed sequence parts used across tests.
_CDS = "AUGCCCAAGUAA"   # AUG + 1 sense codon + stop (12 nt, divisible by 3)
# Derive the 3'UTR from the same pure function used at runtime, so tests
# exercise the actual three_prime code without touching the expression DB.
_UTR3 = generate_mrna_sponge_utr(["UAGCUUAUCAGACUGAUGUUGA"], num_sites=2)["full_utr"]
_UTR5_MIN = 4
_UTR5_MAX = 20

_NULL_RIBONN = {
    "mean_te": 1.0,
    "per_tissue": None,
    "status": "AMBER",
    "message": "mocked",
}


@pytest.fixture(autouse=True)
def mock_ribonn(mocker):
    """Prevent RiboNN from loading GPU models in optimization unit tests."""
    mocker.patch(
        "chainofcustody.evaluation.ribonn.score_ribonn_batch",
        side_effect=lambda seqs: [_NULL_RIBONN] * len(seqs),
    )
    mocker.patch(
        "chainofcustody.evaluation.ribonn.score_ribonn",
        return_value=_NULL_RIBONN,
    )


def _problem(utr5_min: int = _UTR5_MIN, utr5_max: int = _UTR5_MAX) -> SequenceProblem:
    return SequenceProblem(utr5_min=utr5_min, utr5_max=utr5_max, cds=_CDS, utr3=_UTR3)


# ── Problem ──────────────────────────────────────────────────────────────────

def test_problem_dimensions():
    problem = _problem(utr5_min=10, utr5_max=50)
    # n_var = utr5_max + 1 (one length variable + utr5_max nucleotide slots)
    assert problem.n_var == 51
    assert problem.n_obj == N_METRICS


def test_problem_bounds():
    problem = _problem(utr5_min=5, utr5_max=20)
    assert int(problem.xl[0]) == 5   # length lower bound
    assert int(problem.xu[0]) == 20  # length upper bound
    assert int(problem.xl[1]) == 0   # nucleotide lower bound
    assert int(problem.xu[1]) == 3   # nucleotide upper bound


def test_problem_evaluate_shape():
    """Vectorized _evaluate receives a full population matrix; drive via evaluate()."""
    problem = _problem()
    # Column 0 = length in [utr5_min, utr5_max], columns 1+ = nucleotides
    X = np.column_stack([
        np.random.randint(_UTR5_MIN, _UTR5_MAX + 1, size=10),
        np.random.randint(0, 4, size=(10, _UTR5_MAX)),
    ])
    result = problem.evaluate(X)
    assert result.shape == (10, N_METRICS)
    assert np.all((result >= 0) & (result <= 1))


def test_problem_decode():
    problem = _problem(utr5_min=4, utr5_max=4)
    # x[0]=4 (length), x[1:5]=[0,1,2,3] → ACGU
    X = np.array([[4, 0, 1, 2, 3]])
    decoded = problem.decode(X)
    assert decoded == ["ACGU" + KOZAK + _CDS + _UTR3]


def test_problem_decode_uses_length():
    """decode respects x[0] and ignores trailing padding."""
    problem = _problem(utr5_min=2, utr5_max=6)
    # length=2, active=[A,C], padding=[G,U,A,C] (ignored)
    X = np.array([[2, 0, 1, 2, 3, 0, 1]])
    decoded = problem.decode(X)
    assert decoded == ["AC" + KOZAK + _CDS + _UTR3]


# ── Sampling ─────────────────────────────────────────────────────────────────

def test_sampling_shape():
    problem = _problem()
    sampling = NucleotideSampling()
    X = sampling._do(problem, n_samples=20)
    assert X.shape == (20, _UTR5_MAX + 1)
    # Length column within bounds
    assert np.all(X[:, 0] >= _UTR5_MIN) and np.all(X[:, 0] <= _UTR5_MAX)
    # Nucleotide columns in [0, 3]
    assert X[:, 1:].min() >= 0 and X[:, 1:].max() <= 3


# ── Mutation ─────────────────────────────────────────────────────────────────

def test_mutation_preserves_shape():
    problem = _problem()
    mutation = NucleotideMutation(mutation_rate=0.5)
    X = np.column_stack([
        np.full(10, _UTR5_MIN),
        np.zeros((10, _UTR5_MAX), dtype=int),
    ])
    X_mut = mutation._do(problem, X)
    assert X_mut.shape == X.shape
    # Length column stays within bounds
    assert np.all(X_mut[:, 0] >= _UTR5_MIN) and np.all(X_mut[:, 0] <= _UTR5_MAX)
    # Nucleotide columns in [0, 3]
    assert X_mut[:, 1:].min() >= 0 and X_mut[:, 1:].max() <= 3


def test_mutation_rate_zero_is_identity():
    problem = _problem()
    mutation = NucleotideMutation(mutation_rate=0.0)
    X = np.column_stack([
        np.full(5, _UTR5_MIN),
        np.random.randint(0, 4, size=(5, _UTR5_MAX)),
    ])
    assert np.array_equal(mutation._do(problem, X), X)


def test_mutation_invalid_rate():
    with pytest.raises(ValueError):
        NucleotideMutation(mutation_rate=1.5)


# ── End-to-end ────────────────────────────────────────────────────────────────

def test_run_returns_pareto_front():
    X, F, history = run(utr5_min=4, utr5_max=20, cds=_CDS, utr3=_UTR3, pop_size=128, n_gen=3, seed=42)
    assert X.ndim == 2
    assert X.shape[1] == 21  # utr5_max + 1
    assert F.shape[1] == N_METRICS
    # Length variable in bounds
    assert np.all(X[:, 0] >= 4) and np.all(X[:, 0] <= 20)
    assert len(history) > 0
    assert {"generation", "sequence", "overall"} <= history[0].keys()
    # History sequences are full assembled sequences (length varies per individual)
    seq = history[0]["sequence"]
    assert len(seq) >= 4 + len(KOZAK) + len(_CDS) + len(_UTR3)
    assert len(seq) <= 20 + len(KOZAK) + len(_CDS) + len(_UTR3)
