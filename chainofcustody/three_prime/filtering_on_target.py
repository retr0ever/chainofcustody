"""Find a minimal set of miRNAs that silence every cell type *except* a target.

For a given target cell type the algorithm finds miRNAs (by MiRBase_ID) that are
**silent** in the target (mean RPM < ``target_threshold``) but **collectively
active** in every other cell type (mean RPM ≥ ``cover_threshold``).  This is a
greedy weighted-set-cover approach operating at the individual miRNA level
(not aggregated by seed).

Usage
-----
    python -m chainofcustody.three_prime.filtering_on_target --target Hepatocyte_derived
    python -m chainofcustody.three_prime.filtering_on_target --target Astrocyte --target-thresh 50 --cover-thresh 500 --plot
    python -m chainofcustody.three_prime.filtering_on_target --list-targets
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DB_DIR = Path(__file__).resolve().parent / "db"


# ── data loading ────────────────────────────────────────────────────────────

def _load_mature_sequences(
    db_dir: Path = DB_DIR,
) -> tuple[dict[str, str], dict[str, str]]:
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


def load_data(db_dir: Path = DB_DIR) -> tuple[
    dict[str, str],
    dict[str, str],
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, str],
]:
    """Load expression data and build miRNA-level (MiRBase_ID) matrices.

    Returns
    -------
    mature_seqs : dict
        MiRBase_ID → mature sequence.
    seed_seqs : dict
        MiRBase_ID → seed (nt 2-8).
    df_seed_map : DataFrame
        The ``cell_type_seed_map.csv`` contents.
    df_sample_celltype_mir : DataFrame
        MiRBase_ID × sample expression matrix (columns are cell-type labels).
    df_mir_celltype_mean : DataFrame
        MiRBase_ID × cell-type mean-expression matrix.
    mir_to_seed : dict
        MiRBase_ID → seed string.
    """
    mature_seqs, seed_seqs = _load_mature_sequences(db_dir)

    df_seed_map = pd.read_csv(db_dir / "cell_type_seed_map.csv")
    df_samples = pd.read_csv(db_dir / "expression_matrix.csv", index_col=0)
    df_metadata = pd.read_csv(db_dir / "sample_metadata.csv", index_col=0)

    # MiRBase_ID → seed lookup from the seed map
    mir_to_seed: dict[str, str] = dict(
        zip(df_seed_map["MiRBase_ID"], df_seed_map["seed"])
    )

    # Map sample columns → cell-type labels
    map_sample_celltype = df_metadata["CellType"].to_dict()
    df_sample_celltype = df_samples.copy()
    df_sample_celltype.columns = df_sample_celltype.columns.map(map_sample_celltype)

    # Keep only miRNAs with max RPM > 100 in at least one sample
    df_sample_celltype = df_sample_celltype[df_sample_celltype.max(axis=1) > 100]

    # Restrict to miRNAs present in the seed map
    known_mirs = set(df_seed_map["MiRBase_ID"].unique())
    keep_idx = df_sample_celltype.index.intersection(sorted(known_mirs))
    df_sample_celltype_mir = df_sample_celltype.loc[keep_idx]

    # Drop miRNAs with zero total expression
    df_sample_celltype_mir = df_sample_celltype_mir[
        df_sample_celltype_mir.sum(axis=1) > 0
    ]

    # Mean expression per MiRBase_ID per cell type
    df_mir_celltype_mean = df_sample_celltype_mir.T.groupby(level=0).mean().T

    return (
        mature_seqs,
        seed_seqs,
        df_seed_map,
        df_sample_celltype_mir,
        df_mir_celltype_mean,
        mir_to_seed,
    )


# ── core algorithm ──────────────────────────────────────────────────────────

def greedy_mirna_cover(
    target_cell: str,
    df_mir_celltype_mean: pd.DataFrame,
    target_threshold: float = 10.0,
    cover_threshold: float = 1000.0,
    max_mirnas: int = 20,
) -> dict:
    """Find a minimal set of miRNAs (MiRBase_IDs) silent in *target_cell* but
    collectively expressed (≥ *cover_threshold*) in every other cell type.

    Parameters
    ----------
    target_cell : str
        The cell type to protect (miRNAs must be silent here).
    df_mir_celltype_mean : DataFrame
        MiRBase_ID × cell-type mean expression matrix.
    target_threshold : float
        Maximum mean RPM in *target_cell* for a miRNA to be a candidate.
    cover_threshold : float
        Minimum mean RPM in a non-target cell for it to count as "covered".
    max_mirnas : int
        Give up after selecting this many miRNAs.

    Returns
    -------
    dict with keys
        success : bool
        selected_mirnas : list[str]
        covered_per_step : list[set]
        uncovered : set
        all_other_cells : set
    """
    all_cells = set(df_mir_celltype_mean.columns)
    if target_cell not in all_cells:
        available = sorted(all_cells)
        raise ValueError(
            f"Target '{target_cell}' not found.  Available targets:\n"
            + "\n".join(f"  {t}" for t in available)
        )

    other_cells = all_cells - {target_cell}

    # Candidate miRNAs: silent in target
    target_expr = df_mir_celltype_mean[target_cell]
    candidate_mirnas = set(target_expr[target_expr < target_threshold].index)

    # Pre-compute coverage sets per miRNA
    mirna_covers: dict[str, set] = {}
    for mirna in candidate_mirnas:
        mirna_covers[mirna] = {
            cell
            for cell in other_cells
            if df_mir_celltype_mean.loc[mirna, cell] >= cover_threshold
        }

    uncovered = set(other_cells)
    selected: list[str] = []
    covered_per_step: list[set] = []

    while uncovered and len(selected) < max_mirnas:
        best_mirna = None
        best_new: set[str] = set()
        for mirna in candidate_mirnas - set(selected):
            new_cover = mirna_covers[mirna] & uncovered
            if len(new_cover) > len(best_new):
                best_mirna = mirna
                best_new = new_cover

        if best_mirna is None or len(best_new) == 0:
            break

        selected.append(best_mirna)
        covered_per_step.append(best_new)
        uncovered -= best_new

    return {
        "success": len(uncovered) == 0,
        "selected_mirnas": selected,
        "covered_per_step": covered_per_step,
        "uncovered": uncovered,
        "all_other_cells": other_cells,
    }


def build_result_table(
    result: dict,
    target_cell: str,
    df_mir_celltype_mean: pd.DataFrame,
    mir_to_seed: dict[str, str],
    mature_seqs: dict[str, str],
    seed_seqs: dict[str, str],
) -> pd.DataFrame:
    """Build a summary DataFrame for the selected miRNAs.

    Columns: step, seed, target_RPM, n_covered, MiRBase_IDs,
             mature_sequences, seed_sequences.
    """
    rows: list[dict] = []
    for i, (mirna_id, covered) in enumerate(
        zip(result["selected_mirnas"], result["covered_per_step"]), 1
    ):
        seed = mir_to_seed.get(mirna_id, "?")
        rows.append(
            {
                "step": i,
                "seed": seed,
                "target_RPM": df_mir_celltype_mean.loc[mirna_id, target_cell],
                "n_covered": len(covered),
                "MiRBase_IDs": mirna_id,
                "mature_sequences": mature_seqs.get(mirna_id, "?"),
                "seed_sequences": seed_seqs.get(mirna_id, seed),
            }
        )
    return pd.DataFrame(rows)


# ── plotting ────────────────────────────────────────────────────────────────

def plot_on_target(
    target_cell: str,
    selected_mirnas: list[str],
    df_sample_celltype_mir: pd.DataFrame,
    top_n_celltypes: int = 20,
) -> None:
    """Boxplot of combined miRNA expression — bottom *top_n_celltypes* cell
    types by mean, with the target cell highlighted in red."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    combined = df_sample_celltype_mir.loc[
        df_sample_celltype_mir.index.intersection(selected_mirnas)
    ].sum(axis=0)

    plot_df = pd.DataFrame(
        {"Cell type": combined.index, "Combined seed RPM": combined.values}
    )

    # Bottom N cell types by mean expression
    ct_means = plot_df.groupby("Cell type")["Combined seed RPM"].mean()
    bottom_cts = set(ct_means.nsmallest(top_n_celltypes).index)
    bottom_cts.add(target_cell)  # always include target

    plot_df = plot_df[plot_df["Cell type"].isin(bottom_cts)]

    ct_order = (
        plot_df.groupby("Cell type")["Combined seed RPM"]
        .mean()
        .sort_values(ascending=True)
        .index.tolist()
    )

    fig, ax = plt.subplots(
        figsize=(max(12, len(ct_order) * 0.55), 5), constrained_layout=True
    )

    sns.boxplot(
        data=plot_df,
        x="Cell type",
        y="Combined seed RPM",
        order=ct_order,
        ax=ax,
        color="steelblue",
        fliersize=2,
    )
    sns.stripplot(
        data=plot_df,
        x="Cell type",
        y="Combined seed RPM",
        order=ct_order,
        ax=ax,
        color="black",
        size=2,
        alpha=0.4,
    )

    for lbl in ax.get_xticklabels():
        if lbl.get_text() == target_cell:
            lbl.set_color("red")
            lbl.set_fontweight("bold")

    ax.set_title(
        f"Target: {target_cell}  —  {len(selected_mirnas)} miRNA(s): "
        f"{', '.join(selected_mirnas)}  (bottom {top_n_celltypes} by mean)",
        fontsize=11,
        fontweight="bold",
    )
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Combined seed expression (RPM)")
    ax.tick_params(axis="x", rotation=90)
    plt.show()


# ── CLI ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Find a minimal set of miRNAs (by MiRBase_ID) that are silent in "
            "a target cell type but collectively expressed in all other cell "
            "types."
        ),
    )
    parser.add_argument(
        "--target",
        "-t",
        required=False,
        default=None,
        help="Target cell type to protect (e.g. Hepatocyte_derived).",
    )
    parser.add_argument(
        "--target-thresh",
        type=float,
        default=10.0,
        help="Max mean RPM in target for a miRNA to be a candidate (default: 10).",
    )
    parser.add_argument(
        "--cover-thresh",
        type=float,
        default=1000.0,
        help="Min mean RPM in a non-target cell to count as covered (default: 1000).",
    )
    parser.add_argument(
        "--max-seeds",
        type=int,
        default=20,
        help="Maximum number of miRNAs to select (default: 20).",
    )
    parser.add_argument(
        "--list-targets",
        action="store_true",
        help="List available target cell types and exit.",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Show a boxplot of the bottom-20 cell types by combined miRNA expression.",
    )
    parser.add_argument(
        "--plot-top-n",
        type=int,
        default=20,
        help="Number of lowest-expression cell types to show in the plot (default: 20).",
    )
    args = parser.parse_args()

    print("Loading data …")
    (
        mature_seqs,
        seed_seqs,
        df_seed_map,
        df_sample_celltype_mir,
        df_mir_celltype_mean,
        mir_to_seed,
    ) = load_data()

    if args.list_targets:
        targets = sorted(df_mir_celltype_mean.columns)
        print(f"\n{len(targets)} available target cell types:")
        for t in targets:
            print(f"  {t}")
        return

    if args.target is None:
        parser.error("--target is required (use --list-targets to see options)")

    # ── Run greedy set cover ────────────────────────────────────────────────
    result = greedy_mirna_cover(
        target_cell=args.target,
        df_mir_celltype_mean=df_mir_celltype_mean,
        target_threshold=args.target_thresh,
        cover_threshold=args.cover_thresh,
        max_mirnas=args.max_seeds,
    )

    n_other = len(result["all_other_cells"])
    n_covered = n_other - len(result["uncovered"])

    print(f"\nTarget cell     : {args.target}")
    print(f"Silent thresh   : < {args.target_thresh} RPM")
    print(f"Cover thresh    : ≥ {args.cover_thresh} RPM")
    print(f"Success         : {result['success']}")
    print(f"Seeds selected  : {len(result['selected_mirnas'])}")
    print(f"% Off Target   : {(1-n_covered/n_other)*100:.2f}")

    if result["uncovered"]:
        print(f"Uncovered cells : {', '.join(sorted(result['uncovered']))}")

    if not result["selected_mirnas"]:
        print(
            "\nNo candidate miRNAs found — try lowering --target-thresh "
            "or --cover-thresh."
        )
        return

    # ── Result table ────────────────────────────────────────────────────────
    df_result = build_result_table(
        result=result,
        target_cell=args.target,
        df_mir_celltype_mean=df_mir_celltype_mean,
        mir_to_seed=mir_to_seed,
        mature_seqs=mature_seqs,
        seed_seqs=seed_seqs,
    )

    print(
        f"\n{'Step':>4}  {'Seed':<10} {'Target RPM':>10}  "
        f"{'Covered':>7}  {'MiRBase_ID':<26} {'Mature sequence':<28} {'Seed seq':<10}"
    )
    print("-" * 100)

    for _, row in df_result.iterrows():
        mir_ids = [m.strip() for m in row["MiRBase_IDs"].split(",")]
        matures = [m.strip() for m in row["mature_sequences"].split(",")]
        seeds_display = [m.strip() for m in row["seed_sequences"].split(",")]

        for j, (mid, mat, sd) in enumerate(zip(mir_ids, matures, seeds_display)):
            if j == 0:
                print(
                    f"{row['step']:>4}  {row['seed']:<10} "
                    f"{row['target_RPM']:>10.1f}  {row['n_covered']:>7}  "
                    f"{mid:<26} {mat:<28} {sd:<10}"
                )
            else:
                print(
                    f"{'':>4}  {'':10} {'':>10}  {'':>7}  "
                    f"{mid:<26} {mat:<28} {sd:<10}"
                )

    # ── Plot ────────────────────────────────────────────────────────────────
    if args.plot:
        plot_on_target(
            target_cell=args.target,
            selected_mirnas=result["selected_mirnas"],
            df_sample_celltype_mir=df_sample_celltype_mir,
            top_n_celltypes=args.plot_top_n,
        )


if __name__ == "__main__":
    main()
