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


# ── Problem ──────────────────────────────────────────────────────────────────

def test_problem_dimensions():
    problem = SequenceProblem(seq_len=50)
    assert problem.n_var == 50
    assert problem.n_obj == N_METRICS


def test_problem_evaluate_shape():
    """ElementwiseProblem._evaluate receives a single row; drive via evaluate()."""
    problem = SequenceProblem(seq_len=20)
    X = np.random.randint(0, 4, size=(10, 20))
    result = problem.evaluate(X)
    assert result.shape == (10, N_METRICS)
    assert np.all((result >= 0) & (result <= 1))


def test_problem_evaluate_returns_valid_scores():
    problem = SequenceProblem(seq_len=90)
    X = np.random.randint(0, 4, size=(3, 90))
    result = problem.evaluate(X)
    assert np.all(result >= 0)
    assert np.all(result <= 1)


def test_problem_decode():
    problem = SequenceProblem(seq_len=4)
    X = np.array([[0, 1, 2, 3]])  # A C G U
    assert problem.decode(X) == ["ACGU"]


# ── Sampling ─────────────────────────────────────────────────────────────────

def test_sampling_shape():
    problem = SequenceProblem(seq_len=30)
    sampling = NucleotideSampling()
    X = sampling._do(problem, n_samples=20)
    assert X.shape == (20, 30)
    assert X.min() >= 0 and X.max() <= 3


# ── Mutation ─────────────────────────────────────────────────────────────────

def test_mutation_preserves_shape():
    problem = SequenceProblem(seq_len=30)
    mutation = NucleotideMutation(mutation_rate=0.1)
    X = np.zeros((10, 30), dtype=int)
    X_mut = mutation._do(problem, X)
    assert X_mut.shape == X.shape
    assert X_mut.min() >= 0 and X_mut.max() <= 3


def test_mutation_rate_zero_is_identity():
    problem = SequenceProblem(seq_len=30)
    mutation = NucleotideMutation(mutation_rate=0.0)
    X = np.random.randint(0, 4, size=(5, 30))
    assert np.array_equal(mutation._do(problem, X), X)


def test_mutation_invalid_rate():
    with pytest.raises(ValueError):
        NucleotideMutation(mutation_rate=1.5)


# ── End-to-end ────────────────────────────────────────────────────────────────

def test_run_returns_pareto_front():
    X, F = run(seq_len=20, pop_size=20, n_gen=3, seed=42)
    assert X.ndim == 2
    assert X.shape[1] == 20
    assert F.shape[1] == N_METRICS
    assert X.min() >= 0 and X.max() <= 3
