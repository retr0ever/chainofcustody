"""Export the pre-computed miRNA expression data for the coc-web frontend.

Loads the per-cell-type mean expression matrix via the existing
filtering_on_target.load_data() pipeline, then serialises it as
coc-web/public/mirna_data.json.

Run once (or whenever the underlying CSVs change):
    uv run python scripts/export_data.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

# Allow running from the repo root without installing the package
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from chainofcustody.three_prime.filtering_on_target import load_data  # noqa: E402

OUTPUT_PATH = REPO_ROOT.parent / "coc-web" / "public" / "mirna_data.json"


def main() -> None:
    print("Loading expression data …")
    mature_seqs, seed_seqs, _df_seed_map, _df_sample, df_mir_celltype_mean, mir_to_seed = (
        load_data()
    )

    cell_types = sorted(df_mir_celltype_mean.columns.tolist())
    mirnas = df_mir_celltype_mean.index.tolist()

    print(f"  {len(mirnas)} miRNAs × {len(cell_types)} cell types")

    # Build the mean matrix as a nested dict: mirna_id → {cell_type → mean_rpm}
    # Round to 2 decimal places to keep the file small.
    mean_matrix: dict[str, dict[str, float]] = {}
    for mirna in mirnas:
        row = df_mir_celltype_mean.loc[mirna]
        mean_matrix[mirna] = {
            ct: round(float(row[ct]), 2) for ct in cell_types if float(row[ct]) != 0.0
        }

    payload = {
        "cell_types": cell_types,
        "mirnas": mirnas,
        "mean_matrix": mean_matrix,
        "mir_to_seed": {k: v for k, v in mir_to_seed.items() if k in set(mirnas)},
        "mature_seqs": {k: v for k, v in mature_seqs.items() if k in set(mirnas)},
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w") as fh:
        json.dump(payload, fh, separators=(",", ":"))

    size_mb = OUTPUT_PATH.stat().st_size / 1_000_000
    print(f"Written → {OUTPUT_PATH}  ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
