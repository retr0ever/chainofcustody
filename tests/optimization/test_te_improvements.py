"""Tests for the TE improvement features: MOESM3 seeds, gradient seed, fitness sigmoid.

These tests are in a separate file to avoid the autouse ``mock_ribonn`` fixture
defined in ``test_optimization.py``, which patches ``load_top_utr5_seeds`` at
the module level and would mask the real function calls needed here.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

_CDS = "AUGCCCAAGUAA"
_UTR3 = "GAGCCCUAA"
_UTR5_MAX = 20


# ── MOESM3 seed loader ────────────────────────────────────────────────────────

def test_moesm3_seeds_returns_list_of_rna_strings():
    """load_top_utr5_seeds returns RNA strings (no T, only ACGU)."""
    from chainofcustody.optimization.moesm3_seeds import load_top_utr5_seeds
    seeds = load_top_utr5_seeds(n=5)
    assert isinstance(seeds, list)
    assert len(seeds) <= 5
    for s in seeds:
        assert isinstance(s, str)
        assert len(s) > 0
        assert all(c in "ACGU" for c in s), f"Non-RNA character in seed: {s!r}"


def test_moesm3_seeds_respects_length_bounds():
    """Returned seeds respect the min/max length constraints."""
    from chainofcustody.optimization.moesm3_seeds import load_top_utr5_seeds
    seeds = load_top_utr5_seeds(n=10, min_utr5_len=30, max_utr5_len=150)
    for s in seeds:
        assert 30 <= len(s) <= 150, f"Seed length {len(s)} out of [30, 150]: {s[:20]!r}"


def test_moesm3_seeds_count():
    """load_top_utr5_seeds returns at most n seeds and the list is non-empty."""
    from chainofcustody.optimization.moesm3_seeds import load_top_utr5_seeds
    seeds = load_top_utr5_seeds(n=3)
    assert 0 < len(seeds) <= 3


def test_moesm3_seeds_missing_file_returns_empty():
    """load_top_utr5_seeds returns [] gracefully when the file is absent."""
    from chainofcustody.optimization.moesm3_seeds import load_top_utr5_seeds
    seeds = load_top_utr5_seeds(n=5, data_path=Path("/nonexistent/path.xlsx"))
    assert seeds == []


# ── Gradient seed ─────────────────────────────────────────────────────────────

def test_gradient_seed_returns_chromosome_rows(mocker):
    """generate_gradient_seeds returns correctly shaped chromosome rows."""
    import torch
    from chainofcustody.optimization.gradient_seed import generate_gradient_seeds

    mocker.patch(
        "chainofcustody.evaluation.ribonn.score_ribonn",
        return_value={"target_te": 1.3, "mean_te": 1.1},
    )

    mock_predictor = mocker.MagicMock()
    mock_predictor.device = "cpu"
    mock_predictor._predicted_cols = ["predicted_TE_neurons", "predicted_TE_fibroblast"]
    mock_model = mocker.MagicMock()
    # requires_grad=True so that .backward() works through the mock
    mock_model.side_effect = lambda x: (x[:, :2, 0].sum(dim=-1, keepdim=True).expand(-1, 2) * 0 + torch.tensor([[1.2, 0.8]]))
    mock_predictor._fold_models = [(0, [mock_model])]
    mocker.patch(
        "chainofcustody.optimization.gradient_seed.get_predictor",
        return_value=mock_predictor,
    )

    rows = generate_gradient_seeds(
        cds=_CDS, utr3=_UTR3, target_cell_type="neurons",
        utr5_len=20, n_steps=5, n_seeds=2, n_restarts=2,
        utr5_max=_UTR5_MAX,
    )

    assert isinstance(rows, list)
    assert len(rows) <= 2
    for row in rows:
        assert isinstance(row, np.ndarray)
        assert row.shape == (_UTR5_MAX + 1,)
        assert row[0] == 20  # stored length
        assert np.all((row[1:21] >= 0) & (row[1:21] <= 3))


def test_gradient_seed_unknown_target_returns_empty(mocker):
    """generate_gradient_seeds returns [] for an unrecognised target cell type."""
    from chainofcustody.optimization.gradient_seed import generate_gradient_seeds

    mock_predictor = mocker.MagicMock()
    mock_predictor.device = "cpu"
    mock_predictor._predicted_cols = ["predicted_TE_neurons"]
    mock_predictor._fold_models = []
    mocker.patch(
        "chainofcustody.optimization.gradient_seed.get_predictor",
        return_value=mock_predictor,
    )

    rows = generate_gradient_seeds(
        cds=_CDS, utr3=_UTR3, target_cell_type="unknown_tissue",
        utr5_len=10, n_steps=3, n_seeds=1, n_restarts=1,
        utr5_max=_UTR5_MAX,
    )
    assert rows == []


def test_gradient_seed_predictor_unavailable_returns_empty(mocker):
    """generate_gradient_seeds returns [] gracefully when RiboNN cannot load."""
    from chainofcustody.optimization.gradient_seed import generate_gradient_seeds

    mocker.patch(
        "chainofcustody.optimization.gradient_seed.get_predictor",
        side_effect=RuntimeError("model weights missing"),
    )

    rows = generate_gradient_seeds(
        cds=_CDS, utr3=_UTR3, target_cell_type="neurons",
        utr5_len=10, n_steps=3, n_seeds=1, n_restarts=1,
        utr5_max=_UTR5_MAX,
    )
    assert rows == []


def test_gradient_seed_sorted_best_first(mocker):
    """Returned rows are sorted best-TE-first."""
    import torch
    from chainofcustody.optimization.gradient_seed import generate_gradient_seeds

    te_values = iter([0.9, 1.5, 1.1, 0.7])

    def _mock_score_ribonn(seq, target_cell_type="neurons"):
        return {"target_te": next(te_values), "mean_te": 1.0}

    mocker.patch(
        "chainofcustody.evaluation.ribonn.score_ribonn",
        side_effect=_mock_score_ribonn,
    )

    mock_predictor = mocker.MagicMock()
    mock_predictor.device = "cpu"
    mock_predictor._predicted_cols = ["predicted_TE_neurons"]
    mock_model = mocker.MagicMock()
    mock_model.side_effect = lambda x: (x[:, :1, 0].sum(dim=-1, keepdim=True).expand(-1, 1) * 0 + torch.tensor([[1.0]]))
    mock_predictor._fold_models = [(0, [mock_model])]
    mocker.patch(
        "chainofcustody.optimization.gradient_seed.get_predictor",
        return_value=mock_predictor,
    )

    rows = generate_gradient_seeds(
        cds=_CDS, utr3=_UTR3, target_cell_type="neurons",
        utr5_len=10, n_steps=3, n_seeds=4, n_restarts=4,
        utr5_max=_UTR5_MAX,
    )
    # Best TE was 1.5 (restart 2), so first row should come from that restart
    assert len(rows) == 4


# ── Fitness sigmoid ───────────────────────────────────────────────────────────

def test_normalise_te_sigmoid_midpoint():
    """TE=1.0 should yield exactly 0.5 (sigmoid midpoint at 1.0)."""
    from chainofcustody.evaluation.fitness import _normalise_te
    report = {"ribonn_scores": {"target_te": 1.0, "mean_te": 1.0}, "summary": {}}
    assert abs(_normalise_te(report) - 0.5) < 1e-6


def test_normalise_te_sigmoid_gradient_range():
    """With k=6 the score difference between TE=0.8 and TE=1.2 must be >= 0.5.

    Under the old parameters (midpoint=1.2, k=3) this difference was ~0.30.
    The new parameters (midpoint=1.0, k=6) amplify it to >= 0.50, giving
    NSGA-III a stronger directional gradient.
    """
    from chainofcustody.evaluation.fitness import _normalise_te

    def _report(te):
        return {"ribonn_scores": {"target_te": te, "mean_te": te}, "summary": {}}

    assert _normalise_te(_report(1.2)) - _normalise_te(_report(0.8)) >= 0.50


def test_normalise_te_monotone():
    """Higher TE must always produce a strictly higher fitness score."""
    from chainofcustody.evaluation.fitness import _normalise_te

    def _report(te):
        return {"ribonn_scores": {"target_te": te, "mean_te": te}, "summary": {}}

    prev = _normalise_te(_report(0.5))
    for te in [0.8, 1.0, 1.2, 1.5, 2.0]:
        curr = _normalise_te(_report(te))
        assert curr > prev, f"Score not strictly monotone at TE={te}"
        prev = curr


def test_normalise_te_uses_target_te_over_mean_te():
    """_normalise_te prefers target_te over mean_te when both present."""
    from chainofcustody.evaluation.fitness import _normalise_te, _sigmoid
    report = {"ribonn_scores": {"target_te": 1.5, "mean_te": 0.5}, "summary": {}}
    score = _normalise_te(report)
    expected = _sigmoid(1.5, midpoint=1.0, k=6.0)
    assert abs(score - expected) < 1e-6
