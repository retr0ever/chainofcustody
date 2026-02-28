"""Public API for the three_prime package.

The single exported function ``generate_utr3`` takes an off-target cell type
name and returns an optimised 3'UTR sequence that contains sponge sites for the
miRNAs most highly and specifically expressed in that off-target cell type.
"""

from chainofcustody.three_prime.filtering import load_data, mirnas_for_off_target_cell_type
from chainofcustody.three_prime.generate_UTR3 import generate_mrna_sponge_utr


def generate_utr3(
    off_target_cell_type: str,
    threshold: float = 0.0,
    top_n: int = 5,
    num_sites: int = 16,
) -> str:
    """Generate a 3'UTR sponge sequence tailored for a given off-target cell type.

    Identifies the ``top_n`` lowest-entropy miRNAs whose mean expression in
    ``off_target_cell_type`` exceeds ``threshold``, then builds a multi-site
    sponge 3'UTR that sequesters those miRNAs to prevent off-target silencing.

    Parameters
    ----------
    off_target_cell_type : str
        Off-target cell type name as it appears in the expression database
        (e.g. ``"Hepatocyte_derived"``).
    threshold : float
        Minimum mean expression level for a miRNA to be included.
    top_n : int
        Number of miRNAs to include sponge sites for (default 5).
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
        If ``off_target_cell_type`` is not present in the expression database.
    RuntimeError
        If no miRNAs pass the expression ``threshold`` for the given
        off-target cell type.
    """
    mature_seqs, df_mirna_expr, df_grouped = load_data()

    mirnas = mirnas_for_off_target_cell_type(
        off_target_cell_type=off_target_cell_type,
        mature_seqs=mature_seqs,
        df_mirna_expr=df_mirna_expr,
        df_grouped=df_grouped,
        threshold=threshold,
        top_n=top_n,
    )

    if mirnas.empty:
        raise RuntimeError(
            f"No miRNAs found for off-target cell type '{off_target_cell_type}' "
            f"with mean expression >= {threshold}."
        )

    sequences = mirnas["mature_sequence"].tolist()
    result = generate_mrna_sponge_utr(sequences, num_sites=num_sites)
    return result["full_utr"]


__all__ = ["generate_utr3"]
