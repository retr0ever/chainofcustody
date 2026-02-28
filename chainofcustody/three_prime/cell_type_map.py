"""Mapping between cell_type_seed_map.csv names and RiboNN tissue column names.

The expression database (``cell_type_seed_map.csv``) uses title-case names
(e.g. ``"Fibroblast"``), while RiboNN output columns use their own naming
convention (e.g. ``"fibroblast"``).  This module provides a single source of
truth for that mapping and helpers to convert between them.

Only cell types that appear in **both** databases are included.  Cell types
present in only one source cannot be used end-to-end and are therefore omitted.
"""

from __future__ import annotations

# Mapping: seed_map cell type name → RiboNN tissue column name.
# Extend this dict as new overlapping cell types are identified.
SEED_MAP_TO_RIBONN: dict[str, str] = {
    "Fibroblast": "fibroblast",
    "Neuron": "neurons",
}

# Reverse mapping for look-ups in the other direction.
RIBONN_TO_SEED_MAP: dict[str, str] = {v: k for k, v in SEED_MAP_TO_RIBONN.items()}


def seed_map_to_ribonn(cell_type: str) -> str:
    """Convert a ``cell_type_seed_map.csv`` cell type name to its RiboNN column name.

    Raises
    ------
    KeyError
        If *cell_type* is not in the cross-database mapping.  Use
        :data:`SEED_MAP_TO_RIBONN` to see available names.
    """
    try:
        return SEED_MAP_TO_RIBONN[cell_type]
    except KeyError:
        available = ", ".join(sorted(SEED_MAP_TO_RIBONN))
        raise KeyError(
            f"Cell type {cell_type!r} has no RiboNN mapping. "
            f"Available: {available}"
        ) from None


def ribonn_to_seed_map(cell_type: str) -> str:
    """Convert a RiboNN tissue column name to its ``cell_type_seed_map.csv`` name.

    Raises
    ------
    KeyError
        If *cell_type* is not in the cross-database mapping.
    """
    try:
        return RIBONN_TO_SEED_MAP[cell_type]
    except KeyError:
        available = ", ".join(sorted(RIBONN_TO_SEED_MAP))
        raise KeyError(
            f"RiboNN tissue {cell_type!r} has no seed-map mapping. "
            f"Available: {available}"
        ) from None


def get_valid_target_cell_types() -> list[str]:
    """Return all cell type names (seed_map format) that can be used as targets."""
    return sorted(SEED_MAP_TO_RIBONN.keys())
