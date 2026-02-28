"""Load top-TE 5'UTR seeds from the MOESM3 translation efficiency dataset.

The MOESM3_ESM.xlsx supplementary table (Supplementary Data 3 from Karollus
et al., Nature Biotechnology 2024) contains empirical ribosome-profiling TE
values for thousands of human transcripts across 77 cell types, together with
the full transcript sequence and per-region lengths.

:func:`load_top_utr5_seeds` extracts the 5'UTR prefixes of the highest-TE
transcripts (ranked by ``mean_te``) and returns them as RNA strings ready for
population seeding.  The selection is **cell-type-agnostic** — ``mean_te`` is
the arithmetic mean across all 77 cell types, so the seeds represent
universally-high-TE 5'UTRs rather than tissue-specific sequences.
"""

from __future__ import annotations

from pathlib import Path

_DEFAULT_DATA_PATH = Path(__file__).parents[2] / "data" / "MOESM3_ESM.xlsx"

# RiboNN absolute 5'UTR length limit — reject longer seeds
_MAX_UTR5_LEN = 1_381


def load_top_utr5_seeds(
    n: int = 20,
    data_path: Path = _DEFAULT_DATA_PATH,
    max_utr5_len: int = 500,
    min_utr5_len: int = 20,
) -> list[str]:
    """Return the top *n* 5'UTR sequences ranked by mean TE across all cell types.

    Sequences are returned as **RNA** strings (T→U, uppercase).  Only
    transcripts whose 5'UTR length falls in [*min_utr5_len*, *max_utr5_len*]
    are considered; very short or very long UTRs are poor seeds for a
    population seeded at ``initial_length=200``.

    Args:
        n: Maximum number of seed sequences to return.
        data_path: Path to MOESM3_ESM.xlsx.  Defaults to ``data/MOESM3_ESM.xlsx``
            relative to the repo root.
        max_utr5_len: Only include transcripts whose 5'UTR is at most this
            many nt.
        min_utr5_len: Only include transcripts whose 5'UTR is at least this
            many nt.

    Returns:
        List of up to *n* RNA 5'UTR strings, best-first.  Returns an empty
        list if the file is missing or unreadable.
    """
    try:
        import openpyxl  # soft dependency — already in the venv via pandas
    except ImportError:
        return []

    if not data_path.exists():
        return []

    try:
        wb = openpyxl.load_workbook(data_path, read_only=True, data_only=True)
        ws = wb.active

        rows_iter = ws.iter_rows(values_only=True)
        headers = list(next(rows_iter))

        mean_te_idx = headers.index("mean_te")
        tx_seq_idx = headers.index("tx_sequence")
        utr5_size_idx = headers.index("utr5_size")

        candidates: list[tuple[float, str]] = []
        for row in rows_iter:
            mean_te = row[mean_te_idx]
            tx_seq = row[tx_seq_idx]
            utr5_size = row[utr5_size_idx]

            if mean_te is None or tx_seq is None or utr5_size is None:
                continue
            try:
                mean_te = float(mean_te)
                utr5_size = int(utr5_size)
            except (TypeError, ValueError):
                continue

            if not (min_utr5_len <= utr5_size <= min(max_utr5_len, _MAX_UTR5_LEN)):
                continue

            utr5_dna = str(tx_seq)[:utr5_size].upper()
            utr5_rna = utr5_dna.replace("T", "U")

            # Skip sequences containing ambiguous bases
            if any(c not in "ACGU" for c in utr5_rna):
                continue

            candidates.append((mean_te, utr5_rna))

        wb.close()

        candidates.sort(key=lambda x: x[0], reverse=True)
        return [seq for _, seq in candidates[:n]]

    except Exception:
        return []
