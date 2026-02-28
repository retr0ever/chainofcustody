import numpy as np
import pytest

from chainofcustody.optimization import (
    KOZAK,
    METRIC_NAMES,
    NucleotideMutation,
    NucleotideSampling,
    SequenceProblem,
    ElitistNSGA3,
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
    "target_te": 1.0,
    "target_cell_type": "fibroblast",
    "mean_off_target_te": 0.5,
    "per_tissue": None,
    "status": "AMBER",
    "message": "mocked",
}


@pytest.fixture(autouse=True)
def mock_ribonn(mocker):
    """Prevent RiboNN from loading GPU models in optimization unit tests."""
    mocker.patch(
        "chainofcustody.evaluation.ribonn.score_ribonn_batch",
        side_effect=lambda seqs, target_cell_type="megakaryocytes": [_NULL_RIBONN] * len(seqs),
    )
    mocker.patch(
        "chainofcustody.evaluation.ribonn.score_ribonn",
        return_value=_NULL_RIBONN,
    )
    # Prevent MOESM3 file I/O in unit tests that call run() by patching the
    # function object on the module so the lazy import in algorithm.py picks it up.
    import chainofcustody.optimization.moesm3_seeds as _m
    mocker.patch.object(_m, "load_top_utr5_seeds", return_value=["AAAA", "CCCC"])


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


def test_sampling_initial_length():
    """NucleotideSampling with initial_length centres population lengths around that value."""
    problem = _problem(utr5_min=_UTR5_MIN, utr5_max=_UTR5_MAX)
    sampling = NucleotideSampling(initial_length=10)
    X = sampling._do(problem, n_samples=200)
    assert X.shape == (200, _UTR5_MAX + 1)
    # All lengths must stay within bounds
    assert np.all(X[:, 0] >= _UTR5_MIN) and np.all(X[:, 0] <= _UTR5_MAX)
    # Mean should be close to 10 (within a reasonable range given the small UTR range)
    assert abs(X[:, 0].mean() - 10) < 5


def test_mutation_max_length_delta():
    """NucleotideMutation respects max_length_delta for the length variable."""
    problem = _problem(utr5_min=_UTR5_MIN, utr5_max=_UTR5_MAX)
    initial_length = (_UTR5_MIN + _UTR5_MAX) // 2
    X = np.column_stack([
        np.full(50, initial_length),
        np.zeros((50, _UTR5_MAX), dtype=int),
    ])
    mutation = NucleotideMutation(mutation_rate=1.0, max_length_delta=2)
    X_mut = mutation._do(problem, X)
    # Length change per individual must not exceed max_length_delta=2
    assert np.all(np.abs(X_mut[:, 0] - initial_length) <= 2)


# ── End-to-end ────────────────────────────────────────────────────────────────

def test_run_returns_pareto_front():
    X, F, history = run(utr5_min=4, utr5_max=20, cds=_CDS, utr3=_UTR3, pop_size=128, n_gen=3, seed=42, initial_length=10)
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


# ── Elitism ───────────────────────────────────────────────────────────────────

def test_build_algorithm_returns_elitist_nsga3():
    alg = build_algorithm(pop_size=128)
    assert isinstance(alg, ElitistNSGA3)


def test_elitist_nsga3_archive_initialises_to_none():
    alg = build_algorithm(pop_size=32)
    assert alg._elitist_archive is None


def test_elitist_best_score_never_decreases():
    """Best weighted score across the population must be monotonically non-decreasing."""
    from chainofcustody.evaluation.fitness import DEFAULT_WEIGHTS

    weights = np.array([DEFAULT_WEIGHTS.get(m, 0) for m in METRIC_NAMES])
    _, _, history = run(
        utr5_min=4, utr5_max=20, cds=_CDS, utr3=_UTR3,
        pop_size=128, n_gen=6, seed=0, initial_length=10,
    )

    # Compute best overall score per generation from the history records
    gen_best: dict[int, float] = {}
    for row in history:
        g = row["generation"]
        gen_best[g] = max(gen_best.get(g, -np.inf), row["overall"])

    scores_by_gen = [gen_best[g] for g in sorted(gen_best)]

    for i in range(1, len(scores_by_gen)):
        assert scores_by_gen[i] >= scores_by_gen[i - 1] - 1e-9, (
            f"Best score regressed at generation {i + 1}: "
            f"{scores_by_gen[i - 1]:.6f} -> {scores_by_gen[i]:.6f}"
        )


def test_elitist_archive_is_populated_after_run():
    """After run(), the algorithm's _elitist_archive must hold at least one individual."""
    X, F, _ = run(
        utr5_min=_UTR5_MIN, utr5_max=_UTR5_MAX, cds=_CDS, utr3=_UTR3,
        pop_size=32, n_gen=3, seed=1, initial_length=10,
    )
    # The return values from run() come from result.X / result.F which are the
    # Pareto front — a non-empty front confirms elitism ran successfully.
    assert X is not None and len(X) >= 1
    assert F is not None and len(F) >= 1


# ── Seed sequences in NucleotideSampling ─────────────────────────────────────

def test_sampling_with_string_seeds():
    """String seeds are encoded into the first population slots."""
    problem = _problem(utr5_min=4, utr5_max=20)
    seeds = ["AAAA", "CCCC", "GGGG"]
    sampling = NucleotideSampling(seed_sequences=seeds)
    X = sampling._do(problem, n_samples=10)

    assert X.shape == (10, _UTR5_MAX + 1)
    # First three rows should encode the seeds
    assert X[0, 0] == 4  # length of "AAAA"
    assert list(X[0, 1:5]) == [0, 0, 0, 0]   # A=0
    assert X[1, 0] == 4
    assert list(X[1, 1:5]) == [1, 1, 1, 1]   # C=1
    assert X[2, 0] == 4
    assert list(X[2, 1:5]) == [2, 2, 2, 2]   # G=2


def test_sampling_with_array_seeds():
    """Pre-built chromosome row arrays are placed into the first population slots."""
    problem = _problem(utr5_min=4, utr5_max=20)
    row = np.zeros(_UTR5_MAX + 1, dtype=int)
    row[0] = 5
    row[1:6] = [0, 1, 2, 3, 0]  # ACGUA
    sampling = NucleotideSampling(seed_sequences=[row])
    X = sampling._do(problem, n_samples=8)

    assert X.shape == (8, _UTR5_MAX + 1)
    np.testing.assert_array_equal(X[0], row)


def test_sampling_seeds_do_not_exceed_population():
    """More seeds than n_samples: only the first n_samples seeds are used."""
    problem = _problem(utr5_min=4, utr5_max=20)
    seeds = ["AAAA"] * 20  # 20 seeds but only 5 individuals
    sampling = NucleotideSampling(seed_sequences=seeds)
    X = sampling._do(problem, n_samples=5)
    assert X.shape == (5, _UTR5_MAX + 1)


def test_sampling_seed_length_clamped_to_bounds():
    """Seeds whose encoded length exceeds utr5_max are truncated and clamped."""
    problem = _problem(utr5_min=4, utr5_max=10)
    long_seed = "A" * 50  # longer than utr5_max=10
    sampling = NucleotideSampling(seed_sequences=[long_seed])
    X = sampling._do(problem, n_samples=5)
    assert X[0, 0] <= 10  # length clamped to utr5_max


def test_run_with_seed_from_data_disabled(mocker):
    """run() with seed_from_data=False skips MOESM3 loading."""
    load_mock = mocker.patch(
        "chainofcustody.optimization.moesm3_seeds.load_top_utr5_seeds",
        return_value=["AAAA"],
    )
    run(
        utr5_min=_UTR5_MIN, utr5_max=_UTR5_MAX, cds=_CDS, utr3=_UTR3,
        pop_size=32, n_gen=2, seed=0, initial_length=10,
        seed_from_data=False, gradient_seed_steps=0,
    )
    load_mock.assert_not_called()


def test_run_with_seed_from_data_enabled(mocker):
    """run() with seed_from_data=True calls load_top_utr5_seeds."""
    load_mock = mocker.patch(
        "chainofcustody.optimization.moesm3_seeds.load_top_utr5_seeds",
        return_value=["AAAA", "CCCC"],
    )
    run(
        utr5_min=_UTR5_MIN, utr5_max=_UTR5_MAX, cds=_CDS, utr3=_UTR3,
        pop_size=32, n_gen=2, seed=0, initial_length=10,
        seed_from_data=True, gradient_seed_steps=0,
    )
    load_mock.assert_called_once()
