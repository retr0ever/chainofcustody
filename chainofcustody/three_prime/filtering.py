"""Print miRNAs for a fixed cell type filtered by expression threshold and entropy.

Usage
-----
    python -m chainofcustody.three_prime.analysis --cell-type Hepatocyte_derived
    python -m chainofcustody.three_prime.analysis --cell-type Hepatocyte_derived --threshold 500 --top 5
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import entropy as shannon_entropy


DB_DIR = Path(__file__).resolve().parent / "db"


# ── data loading ────────────────────────────────────────────────────────────

def _load_mature_sequences(db_dir: Path = DB_DIR) -> dict[str, str]:
    """Return a mapping MiRBase_ID → mature sequence for human miRNAs."""
    df_family = pd.read_csv(
        db_dir / "miR_Family_Info.txt",
        sep="\t",
        usecols=["Species ID", "MiRBase ID", "Mature sequence"],
    )
    df_human = df_family[df_family["Species ID"] == 9606].drop_duplicates(
        subset="MiRBase ID"
    )
    return dict(zip(df_human["MiRBase ID"], df_human["Mature sequence"]))


def load_data(db_dir: Path = DB_DIR):
    """Return mature-sequence map, sample-level miRNA-expression matrix, and
    grouped cell type × miRNA matrix with Shannon entropy column."""

    mature_seqs = _load_mature_sequences(db_dir)

    df_samples = pd.read_csv(db_dir / "expression_matrix.csv", index_col=0)
    df_metadata = pd.read_csv(db_dir / "sample_metadata.csv", index_col=0)

    # map sample columns → cell type labels
    map_sample_celltype = df_metadata["CellType"].to_dict()
    df_sample_celltype = df_samples.copy()
    df_sample_celltype.columns = df_sample_celltype.columns.map(map_sample_celltype)

    # keep only miRNAs with max count > 100
    df_sample_celltype = df_sample_celltype[df_sample_celltype.max(axis=1) > 100]

    # miRNA × sample expression matrix (index = MiRBase_ID)
    df_mirna_expr = df_sample_celltype[df_sample_celltype.sum(axis=1) > 0]

    # group duplicate cell type columns → mean per cell type
    df_grouped = df_mirna_expr.T.groupby(level=0).mean().T

    # Shannon entropy (base-2) per miRNA across cell types
    row_sums = df_grouped.sum(axis=1)
    prob = df_grouped.div(row_sums.replace(0, np.nan), axis=0)
    df_grouped["shannon_entropy"] = (
        prob.apply(lambda row: shannon_entropy(row.dropna(), base=2), axis=1)
        .fillna(0.0)
    )

    return mature_seqs, df_mirna_expr, df_grouped


# ── core logic ──────────────────────────────────────────────────────────────

def mirnas_for_cell_type(
    cell_type: str,
    mature_seqs: dict[str, str],
    df_mirna_expr: pd.DataFrame,
    df_grouped: pd.DataFrame,
    threshold: float = 0.0,
    top_n: int = 10,
) -> pd.DataFrame:
    """Return the *top_n* lowest-entropy miRNAs whose mean expression in
    *cell_type* is above *threshold*.

    Parameters
    ----------
    cell_type : str
        Cell type name (must match a column in *df_grouped*).
    mature_seqs : dict[str, str]
        Mapping MiRBase_ID → mature sequence.
    df_mirna_expr : DataFrame
        MiRNA × sample expression matrix (columns = cell type labels).
    df_grouped : DataFrame
        MiRNA × cell type matrix with a ``shannon_entropy`` column.
    threshold : float
        Minimum mean expression in *cell_type* for a miRNA to be reported.
    top_n : int
        Number of lowest-entropy miRNAs to return.

    Returns
    -------
    DataFrame
        Columns: MiRBase_ID, mature_sequence, mean_expr, shannon_entropy
        — sorted by entropy ascending.
    """
    if cell_type not in df_grouped.columns:
        available = sorted(c for c in df_grouped.columns if c != "shannon_entropy")
        raise ValueError(
            f"Cell type '{cell_type}' not found. Available cell types:\n"
            + "\n".join(f"  {t}" for t in available)
        )

    # mean expression per miRNA in the chosen cell type
    cell_type_expr = df_grouped[cell_type]

    # filter by threshold
    mask = cell_type_expr >= threshold
    candidates = df_grouped.loc[mask].copy()

    if candidates.empty:
        return pd.DataFrame(
            columns=["MiRBase_ID", "mature_sequence", "mean_expr", "shannon_entropy"]
        )

    # sort by entropy (ascending) and pick top_n
    candidates = candidates.sort_values("shannon_entropy").head(top_n)

    result = pd.DataFrame({
        "MiRBase_ID": candidates.index,
        "mature_sequence": [mature_seqs.get(mid, "") for mid in candidates.index],
        "mean_expr": cell_type_expr.loc[candidates.index].values,
        "shannon_entropy": candidates["shannon_entropy"].values,
    }).reset_index(drop=True)

    return result


# ── plotting ────────────────────────────────────────────────────────────────

def plot_mirnas_boxplot(
    mirna_ids: list[str],
    cell_type: str,
    df_mirna_expr: pd.DataFrame,
    df_grouped: pd.DataFrame,
    top_celltypes: int = 20,
) -> None:
    """Show a boxplot + strip-plot per miRNA.

    Parameters
    ----------
    mirna_ids : list[str]
        MiRBase IDs to plot.
    cell_type : str
        Cell type name (highlighted in the title).
    df_mirna_expr : DataFrame
        MiRNA × sample expression matrix (columns = cell type labels).
    df_grouped : DataFrame
        MiRNA × cell type matrix with ``shannon_entropy`` column.
    top_celltypes : int
        Max cell types to show per panel (ordered by median expression).
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    n = len(mirna_ids)
    fig, axes = plt.subplots(n, 1, figsize=(10, 1 * n), constrained_layout=True,
                             squeeze=False)
    axes = axes.ravel()

    for ax, mid in zip(axes, mirna_ids):
        mir_vals = df_mirna_expr.loc[mid]
        plot_df = pd.DataFrame({
            "Cell type": mir_vals.index,
            "Expression": mir_vals.values,
        })
        plot_df = plot_df[plot_df["Expression"] > 0]

        # keep top N cell types by median expression
        top_ct = (
            plot_df.groupby("Cell type")["Expression"]
            .median()
            .sort_values(ascending=False)
            .head(top_celltypes)
            .index
        )
        plot_df = plot_df[plot_df["Cell type"].isin(top_ct)]

        sns.boxplot(data=plot_df, x="Cell type", y="Expression",
                    ax=ax, color="steelblue")
        sns.stripplot(data=plot_df, x="Cell type", y="Expression",
                      ax=ax, color="black", size=2, alpha=0.5)

        h = df_grouped.loc[mid, "shannon_entropy"]
        ax.set_title(f"{mid}  (entropy = {h:.3f})")
        ax.set_xlabel("Cell type")
        ax.set_ylabel("Expression")
        ax.tick_params(axis="x", rotation=90)

    fig.suptitle(
        f"Lowest-entropy miRNA distributions — cell type: {cell_type}",
        fontsize=13, fontweight="bold", y=1.02,
    )
    plt.show()


# ── CLI ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Print lowest-entropy miRNAs for a given cell type."
    )
    parser.add_argument(
        "--cell-type", "-t", required=True,
        help="Cell type name (e.g. Hepatocyte_derived).",
    )
    parser.add_argument(
        "--threshold", type=float, default=0.0,
        help="Minimum mean expression in the cell type (default: 0).",
    )
    parser.add_argument(
        "--top", "-n", type=int, default=10,
        help="Number of lowest-entropy miRNAs to report (default: 10).",
    )
    parser.add_argument(
        "--list-cell-types", action="store_true",
        help="List available cell type names and exit.",
    )
    parser.add_argument(
        "--plot", action="store_true",
        help="Show debug boxplots of the selected miRNAs across cell types.",
    )
    args = parser.parse_args()

    print("Loading data …")
    mature_seqs, df_mirna_expr, df_grouped = load_data()

    if args.list_cell_types:
        cell_types = sorted(c for c in df_grouped.columns if c != "shannon_entropy")
        print(f"\n{len(cell_types)} available cell types:")
        for t in cell_types:
            print(f"  {t}")
        return

    result = mirnas_for_cell_type(
        cell_type=args.cell_type,
        mature_seqs=mature_seqs,
        df_mirna_expr=df_mirna_expr,
        df_grouped=df_grouped,
        threshold=args.threshold,
        top_n=args.top,
    )

    if result.empty:
        print(f"\nNo miRNAs found for cell type '{args.cell_type}' "
              f"with mean expression >= {args.threshold}.")
        return

    print(f"\nCell type : {args.cell_type}")
    print(f"Threshold : {args.threshold}")
    print(f"Top {args.top} lowest-entropy miRNAs:\n")
    print(f"{'MiRBase_ID':<22} {'Mature sequence':<28} {'Mean expr':>12} {'Entropy (H)':>12}")
    print("-" * 76)
    for _, row in result.iterrows():
        print(
            f"{row['MiRBase_ID']:<22} {row['mature_sequence']:<28} "
            f"{row['mean_expr']:>12.2f} {row['shannon_entropy']:>12.4f}"
        )

    if args.plot:
        plot_mirnas_boxplot(
            mirna_ids=result["MiRBase_ID"].tolist(),
            cell_type=args.cell_type,
            df_mirna_expr=df_mirna_expr,
            df_grouped=df_grouped,
        )


if __name__ == "__main__":
    main()
