"""Tests for the cell_type_seed_map <-> RiboNN name mapping."""
import pytest

from chainofcustody.three_prime.cell_type_map import (
    RIBONN_TO_SEED_MAP,
    SEED_MAP_TO_RIBONN,
    get_valid_target_cell_types,
    ribonn_to_seed_map,
    seed_map_to_ribonn,
)


def test_fibroblast_forward():
    assert seed_map_to_ribonn("Fibroblast") == "fibroblast"


def test_fibroblast_reverse():
    assert ribonn_to_seed_map("fibroblast") == "Fibroblast"


def test_neuron_forward():
    assert seed_map_to_ribonn("Neuron") == "neurons"


def test_neuron_reverse():
    assert ribonn_to_seed_map("neurons") == "Neuron"


def test_round_trip():
    for seed_name, ribonn_name in SEED_MAP_TO_RIBONN.items():
        assert seed_map_to_ribonn(seed_name) == ribonn_name
        assert ribonn_to_seed_map(ribonn_name) == seed_name


def test_unknown_seed_map_raises():
    with pytest.raises(KeyError, match="NotACell"):
        seed_map_to_ribonn("NotACell")


def test_unknown_ribonn_raises():
    with pytest.raises(KeyError, match="not_a_tissue"):
        ribonn_to_seed_map("not_a_tissue")


def test_get_valid_target_cell_types_returns_list():
    valid = get_valid_target_cell_types()
    assert isinstance(valid, list)
    assert len(valid) >= 2
    assert "Fibroblast" in valid
    assert "Neuron" in valid


def test_mappings_are_consistent():
    """Forward and reverse dicts must be exact inverses."""
    assert len(SEED_MAP_TO_RIBONN) == len(RIBONN_TO_SEED_MAP)
    for k, v in SEED_MAP_TO_RIBONN.items():
        assert RIBONN_TO_SEED_MAP[v] == k
