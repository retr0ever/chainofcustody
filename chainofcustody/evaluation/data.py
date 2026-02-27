"""Load MOESM3_ESM.xlsx and compute per-codon translation efficiency weights by cell type."""

from pathlib import Path
from functools import lru_cache

import pandas as pd
import numpy as np

DATA_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "MOESM3_ESM.xlsx"

# All 64 codons
CODONS = [
    f"{a}{b}{c}"
    for a in "ATGC"
    for b in "ATGC"
    for c in "ATGC"
]

# Liver cell line columns in the dataset
LIVER_COLUMNS = ["TE_HepG2", "TE_Huh7", "TE_Huh-7.5"]


@lru_cache(maxsize=1)
def load_dataset(path: str | Path | None = None) -> pd.DataFrame:
    """Load the xlsx dataset. Cached so it's only read once."""
    path = Path(path) if path else DATA_PATH
    if not path.exists():
        raise FileNotFoundError(f"Dataset not found at {path}")
    df = pd.read_excel(path, sheet_name="Human", engine="openpyxl")
    return df


def get_te_columns(df: pd.DataFrame) -> list[str]:
    """Return all TE_ column names from the dataset."""
    return [c for c in df.columns if c.startswith("TE_")]


def _extract_codons_from_cds(tx_seq: str, utr5_size: int, cds_size: int) -> list[str]:
    """Extract codons from a transcript sequence given UTR5 and CDS sizes."""
    cds = tx_seq[utr5_size:utr5_size + cds_size]
    if len(cds) < 3:
        return []
    return [cds[i:i+3] for i in range(0, len(cds) - len(cds) % 3, 3)]


def compute_codon_weights(cell_type_columns: list[str], path: str | Path | None = None) -> dict[str, float]:
    """
    Compute per-codon average TE for given cell type column(s).

    For each codon, we:
    1. Find all genes in the dataset that use that codon
    2. For each usage, record the gene's TE in the specified cell type(s)
    3. Average across all usages

    Returns: dict mapping codon (e.g. "ATG") -> average TE (float)
    """
    df = load_dataset(path)

    # Average across the specified cell type columns
    te_values = df[cell_type_columns].mean(axis=1)

    # Accumulate codon -> list of TE values
    codon_te: dict[str, list[float]] = {c: [] for c in CODONS}

    for _, row in df.iterrows():
        te = te_values[row.name]
        if pd.isna(te):
            continue

        tx_seq = row.get("tx_sequence", "")
        utr5_size = row.get("utr5_size", 0)
        cds_size = row.get("cds_size", 0)

        if not isinstance(tx_seq, str) or not tx_seq:
            continue
        if pd.isna(utr5_size) or pd.isna(cds_size):
            continue

        codons = _extract_codons_from_cds(tx_seq, int(utr5_size), int(cds_size))
        for codon in codons:
            codon_upper = codon.upper()
            if codon_upper in codon_te:
                codon_te[codon_upper].append(te)

    # Average
    result = {}
    for codon, values in codon_te.items():
        if values:
            result[codon] = float(np.mean(values))
        else:
            result[codon] = 0.0

    return result


@lru_cache(maxsize=8)
def get_liver_codon_weights(path: str | Path | None = None) -> dict[str, float]:
    """Get per-codon liver TE weights (cached)."""
    return compute_codon_weights(LIVER_COLUMNS, path)


def get_codon_weights(cell_type: str, path: str | Path | None = None) -> dict[str, float]:
    """
    Get per-codon TE weights for any cell type column.

    Args:
        cell_type: Column name like "HepG2" or "cardiac_fibroblasts".
                   Will be prefixed with "TE_" if not already.
    """
    col = cell_type if cell_type.startswith("TE_") else f"TE_{cell_type}"
    df = load_dataset(path)
    if col not in df.columns:
        available = get_te_columns(df)
        raise ValueError(f"Cell type '{col}' not found. Available: {available}")
    return compute_codon_weights([col], path)
