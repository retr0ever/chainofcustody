"""One-time data preparation script: download and merge miRNA databases.

This script is NOT part of the installed ``chainofcustody`` package.
Run it manually to regenerate the flat-file database used at runtime.

Requirements (not listed in pyproject.toml — install separately):
    - R (≥ 4.0) with Bioconductor and the ``microRNAome`` package
    - rpy2 (pip install rpy2)

Usage
-----
    python scripts/merge_db.py

Outputs written to chainofcustody/three_prime/db/:
    expression_matrix.csv   — miRNA × sample RPM matrix
    sample_metadata.csv     — sample → cell type mapping
    cell_type_seed_map.csv  — pre-joined cell-type × seed table

Pipeline
--------
1. Load the microRNAome Bioconductor dataset (expression matrix + sample
   metadata with ``CellType``) via rpy2.
2. Parse TargetScan's ``miR_Family_Info.txt`` to get a human miRNA → seed
   lookup (species ID 9606).
3. Parse miRBase ``mature.fa`` to get a human miRNA name → accession →
   sequence lookup (useful for cross-referencing).
4. Join expression data with the seed map so that every (cell_type, miRNA)
   pair carries its seed sequence.  The expression matrix is normalised to
   **RPM** (Reads Per Million) before any downstream use.  Only miRNAs with
   RPM > 10 in at least one sample of the cell type are retained.
"""

from pathlib import Path

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parent.parent
DB_DIR = _REPO_ROOT / "chainofcustody" / "three_prime" / "db"
FAMILY_INFO_PATH = DB_DIR / "miR_Family_Info.txt"
MATURE_FA_PATH = DB_DIR / "mature.fa"
HUMAN_SPECIES_ID = 9606

# Classes that are not true human cells (acellular / non-somatic)
EXCLUDED_CLASSES = {"Plasma", "Sperm"}


# ---------------------------------------------------------------------------
# 1. microRNAome expression data & metadata (via R / Bioconductor)
# ---------------------------------------------------------------------------
def load_microRNAome() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Download (once) and return the microRNAome dataset, RPM-normalised.

    Returns
    -------
    expr_df : pd.DataFrame
        RPM-normalised counts – rows are miRNAs, columns are SRR sample IDs.
        RPM = (raw_count / library_size) × 1 000 000.
    meta_df : pd.DataFrame
        Per-sample metadata including ``CellType``.
    """
    print("Setting up R environment and downloading Bioconductor package...")

    ro.r('''
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos="http://cran.us.r-project.org")
    }
    BiocManager::install("microRNAome", update=FALSE, ask=FALSE)
    ''')

    print("Loading data in R...")

    ro.r('''
    library(microRNAome)
    data("microRNAome")
    expr_matrix <- as.data.frame(assay(microRNAome))
    metadata    <- as.data.frame(colData(microRNAome))
    ''')

    print("Converting to Python Pandas DataFrames...")
    with localconverter(ro.default_converter + pandas2ri.converter):
        expr_df = ro.conversion.rpy2py(ro.globalenv['expr_matrix'])
        meta_df = ro.conversion.rpy2py(ro.globalenv['metadata'])

    # ---- RPM normalisation ------------------------------------------------
    library_sizes = expr_df.sum(axis=0)           # total counts per sample
    expr_df = expr_df.div(library_sizes, axis=1) * 1e6
    print("Expression matrix normalised to RPM.")

    return expr_df, meta_df


# ---------------------------------------------------------------------------
# 2. TargetScan miR_Family_Info.txt → human miRNA-to-seed lookup
# ---------------------------------------------------------------------------
def load_human_seed_map(
    family_info_path: Path = FAMILY_INFO_PATH,
) -> pd.DataFrame:
    """Return a DataFrame of human miRNAs and their seed (nt 2-8) sequences.

    Columns: ``miR_family``, ``seed``, ``MiRBase_ID``, ``MiRBase_Accession``.
    """
    df = pd.read_csv(family_info_path, sep="\t", dtype=str)

    # Normalise column names: strip whitespace, replace spaces/symbols
    df.columns = [c.strip().replace(" ", "_").replace("+", "_") for c in df.columns]

    human = df[df["Species_ID"] == str(HUMAN_SPECIES_ID)].copy()
    human = human[["miR_family", "Seed_m8", "MiRBase_ID", "MiRBase_Accession"]]
    human = human.rename(columns={"Seed_m8": "seed"}).reset_index(drop=True)
    return human


# ---------------------------------------------------------------------------
# 3. miRBase mature.fa → name / accession / sequence (human only)
# ---------------------------------------------------------------------------
def parse_mature_fa(fasta_path: Path = MATURE_FA_PATH) -> pd.DataFrame:
    """Parse miRBase ``mature.fa`` and return human entries.

    Columns: ``mirna_name``, ``accession``, ``description``, ``sequence``.
    """
    records: list[tuple[str, str | None, str | None, str]] = []
    name = accession = desc = None
    seq_parts: list[str] = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, accession, desc, "".join(seq_parts)))
                parts = line[1:].split(maxsplit=2)
                name = parts[0]
                accession = parts[1] if len(parts) > 1 else None
                desc = parts[2] if len(parts) > 2 else None
                seq_parts = []
            else:
                seq_parts.append(line)
        if name is not None:
            records.append((name, accession, desc, "".join(seq_parts)))

    fa_df = pd.DataFrame(
        records, columns=["mirna_name", "accession", "description", "sequence"]
    )
    return fa_df[fa_df["mirna_name"].str.startswith("hsa-")].reset_index(drop=True)


# ---------------------------------------------------------------------------
# 4. Build the cell-type ↔ seed mapping
# ---------------------------------------------------------------------------
def build_cell_type_seed_map(
    expr_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    seed_map: pd.DataFrame,
) -> pd.DataFrame:
    """For each cell type, find expressed miRNAs and annotate with their seed.

    A miRNA is considered *expressed* in a cell type when its RPM is > 10
    in **at least one** sample belonging to that cell type.

    For every miRNA the distribution of RPM values across the experiments
    (samples) belonging to that cell type is summarised with: mean, median,
    standard deviation, min, max, the number of samples that express it
    (RPM > 0), and the total number of samples in the cell type.

    Returns
    -------
    pd.DataFrame
        Columns: ``cell_type``, ``MiRBase_ID``, ``miR_family``, ``seed``,
        ``mean_count``, ``median_count``, ``std_count``, ``min_count``,
        ``max_count``, ``n_expressing``, ``n_samples``.
    """
    # Quick lookup: MiRBase ID → (family, seed)
    id_to_seed = (
        seed_map.set_index("MiRBase_ID")[["miR_family", "seed"]]
        .loc[lambda d: ~d.index.duplicated(keep="first")]
    )

    rows: list[dict] = []
    for cell_type in meta_df["CellType"].unique():
        sample_ids = meta_df.loc[meta_df["CellType"] == cell_type].index
        subset = expr_df[expr_df.columns.intersection(sample_ids)]
        n_samples = subset.shape[1]
        expressed_mirnas = subset.index[subset.max(axis=1) > 10]

        for mirna in expressed_mirnas:
            if mirna in id_to_seed.index:
                family, seed = id_to_seed.loc[mirna]
                counts = subset.loc[mirna]
                rows.append(
                    {
                        "cell_type": cell_type,
                        "MiRBase_ID": mirna,
                        "miR_family": family,
                        "seed": seed,
                        "mean_count": counts.mean(),
                        "median_count": counts.median(),
                        "std_count": counts.std(),
                        "min_count": counts.min(),
                        "max_count": counts.max(),
                        "n_expressing": int((counts > 0).sum()),
                        "n_samples": n_samples,
                    }
                )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main – download, match, and save
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # ---- 1. microRNAome expression + metadata ----------------------------
    expr_df, meta_df = load_microRNAome()
    print(
        f"\nExpression matrix : {expr_df.shape[0]} miRNAs × {expr_df.shape[1]} samples"
    )

    # Filter to human cells only (drop acellular / non-somatic classes)
    cell_mask = ~meta_df["Class"].isin(EXCLUDED_CLASSES)
    meta_df = meta_df[cell_mask]
    expr_df = expr_df[expr_df.columns.intersection(meta_df.index)]
    print(
        f"After filtering   : {expr_df.shape[1]} samples, "
        f"{meta_df['CellType'].nunique()} cell types "
        f"(excluded {', '.join(sorted(EXCLUDED_CLASSES))})"
    )
    print(meta_df[["CellType"]].head(10))

    # ---- 2. Human seed map from TargetScan -------------------------------
    seed_map = load_human_seed_map()
    print(
        f"\nHuman seed map    : {len(seed_map)} miRNA entries, "
        f"{seed_map['seed'].nunique()} unique seeds"
    )
    print(seed_map.head(10))

    # ---- 3. mature.fa (optional enrichment) ------------------------------
    if MATURE_FA_PATH.exists():
        mature_df = parse_mature_fa()
        print(f"\nMature miRNAs (hsa): {len(mature_df)} entries")
        print(mature_df.head(10))

    # ---- 4. Cell-type ↔ seed join ----------------------------------------
    ct_seed = build_cell_type_seed_map(expr_df, meta_df, seed_map)
    print(f"\nCell-type × seed associations: {len(ct_seed)} rows")
    print(ct_seed.head(20))

    # ---- Save all outputs ------------------------------------------------
    output_path = DB_DIR / "cell_type_seed_map.csv"
    ct_seed.to_csv(output_path, index=False)
    print(f"\nSaved cell-type ↔ seed map to {output_path}")

    expr_path = DB_DIR / "expression_matrix.csv"
    expr_df.to_csv(expr_path)
    print(f"Saved expression matrix to {expr_path}")

    meta_path = DB_DIR / "sample_metadata.csv"
    meta_df.to_csv(meta_path)
    print(f"Saved sample metadata to {meta_path}")
