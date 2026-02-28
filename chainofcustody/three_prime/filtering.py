"""Print miRNAs for a fixed off-target cell type filtered by expression threshold and entropy.

Usage
-----
    python -m chainofcustody.three_prime.filtering --off-target-cell-type Hepatocyte_derived
    python -m chainofcustody.three_prime.filtering --off-target-cell-type Hepatocyte_derived --threshold 500 --top 5
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import entropy as shannon_entropy


DB_DIR = Path(__file__).resolve().parent / "db"


# ── data loading ────────────────────────────────────────────────────────────

def _load_mature_sequences(db_dir: Path = DB_DIR) -> tuple[dict[str, str], dict[str, str]]:
    """Return mappings MiRBase_ID → mature sequence and MiRBase_ID → seed."""
    df_family = pd.read_csv(
        db_dir / "miR_Family_Info.txt",
        sep="\t",
        usecols=["Species ID", "MiRBase ID", "Mature sequence", "Seed+m8"],
    )
    df_human = df_family[df_family["Species ID"] == 9606].drop_duplicates(
        subset="MiRBase ID"
    )
    mature_seqs = dict(zip(df_human["MiRBase ID"], df_human["Mature sequence"]))
    seed_seqs = dict(zip(df_human["MiRBase ID"], df_human["Seed+m8"]))
    return mature_seqs, seed_seqs


def load_data(db_dir: Path = DB_DIR):
    """Return mature-sequence map, sample-level miRNA-expression matrix, and
    grouped off-target cell type × miRNA matrix with Shannon entropy column."""

    mature_seqs, seed_seqs = _load_mature_sequences(db_dir)

    df_samples = pd.read_csv(db_dir / "expression_matrix.csv", index_col=0)
    df_metadata = pd.read_csv(db_dir / "sample_metadata.csv", index_col=0)

    # map sample columns → off-target cell type labels
    map_sample_celltype = df_metadata["CellType"].to_dict()
    df_sample_celltype = df_samples.copy()
    df_sample_celltype.columns = df_sample_celltype.columns.map(map_sample_celltype)

    # keep only miRNAs with max count > 100
    df_sample_celltype = df_sample_celltype[df_sample_celltype.max(axis=1) > 100]

    # miRNA × sample expression matrix (index = MiRBase_ID)
    df_mirna_expr = df_sample_celltype[df_sample_celltype.sum(axis=1) > 0]

    # group duplicate off-target cell type columns → mean per off-target cell type
    df_grouped = df_mirna_expr.T.groupby(level=0).mean().T

    # Shannon entropy (base-2) per miRNA across off-target cell types
    row_sums = df_grouped.sum(axis=1)
    prob = df_grouped.div(row_sums.replace(0, np.nan), axis=0)
    df_grouped["shannon_entropy"] = (
        prob.apply(lambda row: shannon_entropy(row.dropna(), base=2), axis=1)
        .fillna(0.0)
    )

    return mature_seqs, seed_seqs, df_mirna_expr, df_grouped


# ── core logic ──────────────────────────────────────────────────────────────

def mirnas_for_off_target_cell_type(
    off_target_cell_type: str,
    mature_seqs: dict[str, str],
    seed_seqs: dict[str, str],
    df_mirna_expr: pd.DataFrame,
    df_grouped: pd.DataFrame,
    threshold: float = 0.0,
    top_n: int = 10,
) -> pd.DataFrame:
    """Return the *top_n* lowest-entropy miRNAs whose mean expression in
    *off_target_cell_type* is above *threshold*.

    Parameters
    ----------
    off_target_cell_type : str
        Off-target cell type name (must match a column in *df_grouped*).
    mature_seqs : dict[str, str]
        Mapping MiRBase_ID → mature sequence.
    seed_seqs : dict[str, str]
        Mapping MiRBase_ID → seed (nt 2-8) sequence.
    df_mirna_expr : DataFrame
        MiRNA × sample expression matrix (columns = off-target cell type labels).
    df_grouped : DataFrame
        MiRNA × off-target cell type matrix with a ``shannon_entropy`` column.
    threshold : float
        Minimum mean expression in *off_target_cell_type* for a miRNA to be reported.
    top_n : int
        Number of lowest-entropy miRNAs to return.

    Returns
    -------
    DataFrame
        Columns: MiRBase_ID, mature_sequence, seed, mean_expr, shannon_entropy
        — sorted by entropy ascending.
    """
    if off_target_cell_type not in df_grouped.columns:
        available = sorted(c for c in df_grouped.columns if c != "shannon_entropy")
        raise ValueError(
            f"Off-target cell type '{off_target_cell_type}' not found. "
            f"Available off-target cell types:\n"
            + "\n".join(f"  {t}" for t in available)
        )

    # mean expression per miRNA in the chosen off-target cell type
    off_target_cell_type_expr = df_grouped[off_target_cell_type]

    # filter by threshold
    mask = off_target_cell_type_expr >= threshold
    candidates = df_grouped.loc[mask].copy()

    if candidates.empty:
        return pd.DataFrame(
            columns=["MiRBase_ID", "mature_sequence", "seed", "mean_expr", "shannon_entropy"]
        )

    # sort by entropy (ascending) and pick top_n
    candidates = candidates.sort_values("shannon_entropy").head(top_n)

    result = pd.DataFrame({
        "MiRBase_ID": candidates.index,
        "mature_sequence": [mature_seqs.get(mid, "") for mid in candidates.index],
        "seed": [seed_seqs.get(mid, "") for mid in candidates.index],
        "mean_expr": off_target_cell_type_expr.loc[candidates.index].values,
        "shannon_entropy": candidates["shannon_entropy"].values,
    }).reset_index(drop=True)

    # Drop rows whose sequence is missing from the lookup — they cannot be used
    # downstream to build sponge sites and would cause a ValueError there.
    result = result[result["mature_sequence"] != ""].reset_index(drop=True)

    return result


# ── plotting ────────────────────────────────────────────────────────────────

def plot_mirnas_boxplot(
    mirna_ids: list[str],
    off_target_cell_type: str,
    df_mirna_expr: pd.DataFrame,
    df_grouped: pd.DataFrame,
    top_celltypes: int = 20,
) -> None:
    """Show a boxplot + strip-plot per miRNA.

    Parameters
    ----------
    mirna_ids : list[str]
        MiRBase IDs to plot.
    off_target_cell_type : str
        Off-target cell type name (highlighted in the title).
    df_mirna_expr : DataFrame
        MiRNA × sample expression matrix (columns = off-target cell type labels).
    df_grouped : DataFrame
        MiRNA × off-target cell type matrix with ``shannon_entropy`` column.
    top_celltypes : int
        Max off-target cell types to show per panel (ordered by median expression).
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
            "Off-target cell type": mir_vals.index,
            "Expression": mir_vals.values,
        })
        plot_df = plot_df[plot_df["Expression"] > 0]

        # keep top N off-target cell types by median expression
        top_ct = (
            plot_df.groupby("Off-target cell type")["Expression"]
            .median()
            .sort_values(ascending=False)
            .head(top_celltypes)
            .index
        )
        plot_df = plot_df[plot_df["Off-target cell type"].isin(top_ct)]

        sns.boxplot(data=plot_df, x="Off-target cell type", y="Expression",
                    ax=ax, color="steelblue")
        sns.stripplot(data=plot_df, x="Off-target cell type", y="Expression",
                      ax=ax, color="black", size=2, alpha=0.5)

        h = df_grouped.loc[mid, "shannon_entropy"]
        ax.set_title(f"{mid}  (entropy = {h:.3f})")
        ax.set_xlabel("Off-target cell type")
        ax.set_ylabel("Expression")
        ax.tick_params(axis="x", rotation=90)

    fig.suptitle(
        f"Lowest-entropy miRNA distributions — off-target cell type: {off_target_cell_type}",
        fontsize=13, fontweight="bold", y=1.02,
    )
    plt.show()


# ── CLI ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Print lowest-entropy miRNAs for a given off-target cell type."
    )
    parser.add_argument(
        "--off-target-cell-type", "-t", required=True,
        help="Off-target cell type name (e.g. Hepatocyte_derived).",
    )
    parser.add_argument(
        "--threshold", type=float, default=0.0,
        help="Minimum mean expression in the off-target cell type (default: 0).",
    )
    parser.add_argument(
        "--top", "-n", type=int, default=10,
        help="Number of lowest-entropy miRNAs to report (default: 10).",
    )
    parser.add_argument(
        "--list-off-target-cell-types", action="store_true",
        help="List available off-target cell type names and exit.",
    )
    parser.add_argument(
        "--plot", action="store_true",
        help="Show debug boxplots of the selected miRNAs across off-target cell types.",
    )
    args = parser.parse_args()

    print("Loading data …")
    mature_seqs, seed_seqs, df_mirna_expr, df_grouped = load_data()

    if args.list_off_target_cell_types:
        off_target_cell_types = sorted(c for c in df_grouped.columns if c != "shannon_entropy")
        print(f"\n{len(off_target_cell_types)} available off-target cell types:")
        for t in off_target_cell_types:
            print(f"  {t}")
        return

    result = mirnas_for_off_target_cell_type(
        off_target_cell_type=args.off_target_cell_type,
        mature_seqs=mature_seqs,
        seed_seqs=seed_seqs,
        df_mirna_expr=df_mirna_expr,
        df_grouped=df_grouped,
        threshold=args.threshold,
        top_n=args.top,
    )

    if result.empty:
        print(f"\nNo miRNAs found for off-target cell type '{args.off_target_cell_type}' "
              f"with mean expression >= {args.threshold}.")
        return

    print(f"\nOff-target cell type : {args.off_target_cell_type}")
    print(f"Threshold            : {args.threshold}")
    print(f"Top {args.top} lowest-entropy miRNAs:\n")
    print(f"{'MiRBase_ID':<22} {'Mature sequence':<28} {'Seed':<10} {'Mean expr':>12} {'Entropy (H)':>12}")
    print("-" * 86)
    for _, row in result.iterrows():
        print(
            f"{row['MiRBase_ID']:<22} {row['mature_sequence']:<28} "
            f"{row['seed']:<10} "
            f"{row['mean_expr']:>12.2f} {row['shannon_entropy']:>12.4f}"
        )

    if args.plot:
        plot_mirnas_boxplot(
            mirna_ids=result["MiRBase_ID"].tolist(),
            off_target_cell_type=args.off_target_cell_type,
            df_mirna_expr=df_mirna_expr,
            df_grouped=df_grouped,
        )


if __name__ == "__main__":
    main()
