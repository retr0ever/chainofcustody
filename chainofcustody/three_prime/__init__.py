"""Public API for the three_prime package.

:func:`generate_utr3` generates a 3'UTR sponge that captures miRNAs highly
expressed in a *target* cell type — suppressing expression everywhere that
cell type is not active (on-target protection strategy).

The greedy set-cover algorithm in :mod:`filtering_on_target` selects miRNAs
that are **silent** in the target but **expressed** in every other cell type,
so sponging those miRNAs de-represses translation specifically in the target
while leaving off-target expression unaffected.
"""

from chainofcustody.three_prime.filtering_on_target import (
    greedy_mirna_cover,
    load_data,
)
from chainofcustody.three_prime.generate_utr3 import generate_mrna_sponge_utr


def generate_utr3(
    target_cell_type: str,
    target_threshold: float = 10.0,
    cover_threshold: float = 1000.0,
    max_mirnas: int = 5,
    num_sites: int = 16,
) -> str:
    """Generate a 3'UTR sponge sequence tailored for a given *target* cell type.

    Uses a greedy weighted-set-cover algorithm to find miRNAs that are
    **silent** in ``target_cell_type`` (mean RPM < ``target_threshold``) but
    **expressed** in every other cell type (mean RPM ≥ ``cover_threshold``).
    Sponging these miRNAs de-represses translation in the target tissue while
    keeping off-target silencing intact.

    Parameters
    ----------
    target_cell_type : str
        Target cell type name as it appears in ``cell_type_seed_map.csv``
        (e.g. ``"Fibroblast"``).
    target_threshold : float
        Maximum mean RPM in the target for a miRNA to be a candidate.
    cover_threshold : float
        Minimum mean RPM in a non-target cell for it to count as "covered".
    max_mirnas : int
        Maximum number of miRNAs to include sponge sites for (default 5).
    num_sites : int
        Total number of sponge site repeats in the 3'UTR cassette (default 16).

    Returns
    -------
    str
        Full 3'UTR sequence (RNA, mixed case) including stop codon, sponge
        cassette, and poly-A signal.

    Raises
    ------
    ValueError
        If ``target_cell_type`` is not present in the expression database, or
        if no candidate miRNAs are found.
    """
    (
        mature_seqs,
        _seed_seqs,
        _df_seed_map,
        _df_sample_celltype_mir,
        df_mir_celltype_mean,
        _mir_to_seed,
    ) = load_data()

    result = greedy_mirna_cover(
        target_cell=target_cell_type,
        df_mir_celltype_mean=df_mir_celltype_mean,
        target_threshold=target_threshold,
        cover_threshold=cover_threshold,
        max_mirnas=max_mirnas,
    )

    selected = result["selected_mirnas"]
    if not selected:
        raise ValueError(
            f"No candidate miRNAs found for target '{target_cell_type}'. "
            "Try lowering --target-threshold or --cover-threshold."
        )

    sequences = [mature_seqs[m] for m in selected if m in mature_seqs]
    if not sequences:
        raise ValueError(
            f"Selected miRNAs for '{target_cell_type}' have no known mature sequences."
        )

    sponge = generate_mrna_sponge_utr(sequences, num_sites=num_sites)
    return sponge["full_utr"]


__all__ = ["generate_utr3"]
