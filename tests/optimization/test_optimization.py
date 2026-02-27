import numpy as np
import pytest

from chainofcustody.optimization import (
    NucleotideMutation,
    NucleotideSampling,
    SequenceProblem,
    build_algorithm,
    run,
)


# ── Problem ──────────────────────────────────────────────────────────────────

def test_problem_dimensions():
    problem = SequenceProblem(seq_len=50)
    assert problem.n_var == 50
    assert problem.n_obj == 3


def test_problem_evaluate_shape():
    problem = SequenceProblem(seq_len=20)
    X = np.random.randint(0, 4, size=(10, 20))
    out = {}
    problem._evaluate(X, out)
    assert out["F"].shape == (10, 3)
    assert np.all((out["F"] >= 0) & (out["F"] <= 1))


def test_problem_evaluate_fractions():
    problem = SequenceProblem(seq_len=4)
    # All-A sequence: frac_a=1, frac_c=0, frac_t=0
    X = np.array([[0, 0, 0, 0]])
    out = {}
    problem._evaluate(X, out)
    np.testing.assert_array_almost_equal(out["F"], [[1.0, 0.0, 0.0]])

    # One of each A, C, G, T: each fraction = 0.25
    X = np.array([[0, 1, 2, 3]])
    problem._evaluate(X, out)
    np.testing.assert_array_almost_equal(out["F"], [[0.25, 0.25, 0.25]])


def test_problem_decode():
    problem = SequenceProblem(seq_len=4)
    X = np.array([[0, 1, 2, 3]])  # A C G T
    assert problem.decode(X) == ["ACGT"]


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
    assert F.shape[1] == 3
    assert X.min() >= 0 and X.max() <= 3
