import numpy as np
import pytest

from chainofcustody.optimization import (
    METRIC_NAMES,
    NucleotideMutation,
    NucleotideSampling,
    SequenceProblem,
    build_algorithm,
    run,
)

N_METRICS = len(METRIC_NAMES)

# Minimal fixed sequence parts used across tests.
_CDS = "AUGCCCAAGUAA"   # AUG + 1 sense codon + stop (12 nt, divisible by 3)
_UTR3 = "GAGTAGUCCC"


def _problem(utr5_len: int = 20) -> SequenceProblem:
    return SequenceProblem(utr5_len=utr5_len, cds=_CDS, utr3=_UTR3)


# ── Problem ──────────────────────────────────────────────────────────────────

def test_problem_dimensions():
    problem = _problem(utr5_len=50)
    assert problem.n_var == 50
    assert problem.n_obj == N_METRICS


def test_problem_evaluate_shape():
    """ElementwiseProblem._evaluate receives a single row; drive via evaluate()."""
    problem = _problem(utr5_len=20)
    X = np.random.randint(0, 4, size=(10, 20))
    result = problem.evaluate(X)
    assert result.shape == (10, N_METRICS)
    assert np.all((result >= 0) & (result <= 1))


def test_problem_evaluate_returns_valid_scores():
    problem = _problem(utr5_len=20)
    X = np.random.randint(0, 4, size=(3, 20))
    result = problem.evaluate(X)
    assert np.all(result >= 0)
    assert np.all(result <= 1)


def test_problem_decode():
    problem = _problem(utr5_len=4)
    X = np.array([[0, 1, 2, 3]])  # A C G U
    decoded = problem.decode(X)
    assert decoded == ["ACGU" + _CDS + _UTR3]


# ── Sampling ─────────────────────────────────────────────────────────────────

def test_sampling_shape():
    problem = _problem(utr5_len=30)
    sampling = NucleotideSampling()
    X = sampling._do(problem, n_samples=20)
    assert X.shape == (20, 30)
    assert X.min() >= 0 and X.max() <= 3


# ── Mutation ─────────────────────────────────────────────────────────────────

def test_mutation_preserves_shape():
    problem = _problem(utr5_len=30)
    mutation = NucleotideMutation(mutation_rate=0.1)
    X = np.zeros((10, 30), dtype=int)
    X_mut = mutation._do(problem, X)
    assert X_mut.shape == X.shape
    assert X_mut.min() >= 0 and X_mut.max() <= 3


def test_mutation_rate_zero_is_identity():
    problem = _problem(utr5_len=30)
    mutation = NucleotideMutation(mutation_rate=0.0)
    X = np.random.randint(0, 4, size=(5, 30))
    assert np.array_equal(mutation._do(problem, X), X)


def test_mutation_invalid_rate():
    with pytest.raises(ValueError):
        NucleotideMutation(mutation_rate=1.5)


# ── End-to-end ────────────────────────────────────────────────────────────────

def test_run_returns_pareto_front():
    X, F, history = run(utr5_len=20, cds=_CDS, utr3=_UTR3, pop_size=128, n_gen=3, seed=42)
    assert X.ndim == 2
    assert X.shape[1] == 20
    assert F.shape[1] == N_METRICS
    assert X.min() >= 0 and X.max() <= 3
    assert len(history) > 0
    assert {"generation", "sequence", "overall"} <= history[0].keys()
    # History sequences should be the full assembled sequence
    assert len(history[0]["sequence"]) == 20 + len(_CDS) + len(_UTR3)
