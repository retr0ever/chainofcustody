"""Tests for chainofcustody.three_prime.filtering — pure-logic functions only.

All tests use in-memory DataFrames so no database files or network access are
needed.
"""

import numpy as np
import pandas as pd
import pytest

from chainofcustody.three_prime.filtering import mirnas_for_off_target_cell_type


# ── fixtures ─────────────────────────────────────────────────────────────────

def _make_data(
    mirna_ids: list[str],
    off_target_cell_types: list[str],
    expr_values: list[list[float]],
    entropies: list[float],
) -> tuple[dict, dict, pd.DataFrame, pd.DataFrame]:
    """Build minimal (mature_seqs, seed_seqs, df_mirna_expr, df_grouped) fixtures."""
    mature_seqs = {mid: f"ACGU{'GC' * i}" for i, mid in enumerate(mirna_ids)}
    seed_seqs = {mid: f"ACGU{i}" for i, mid in enumerate(mirna_ids)}

    # df_mirna_expr: miRNA × sample (one sample per off-target cell type for simplicity)
    df_mirna_expr = pd.DataFrame(
        expr_values,
        index=mirna_ids,
        columns=off_target_cell_types,
    )

    # df_grouped: miRNA × off-target cell type + shannon_entropy column
    df_grouped = pd.DataFrame(
        expr_values,
        index=mirna_ids,
        columns=off_target_cell_types,
    )
    df_grouped["shannon_entropy"] = entropies

    return mature_seqs, seed_seqs, df_mirna_expr, df_grouped


@pytest.fixture()
def three_mirna_data():
    mirna_ids              = ["hsa-miR-1", "hsa-miR-2", "hsa-miR-3"]
    off_target_cell_types  = ["Liver", "Brain"]
    expr                   = [
        [900.0, 10.0],   # hsa-miR-1: high in Liver
        [200.0, 800.0],  # hsa-miR-2: expressed in both
        [50.0,   5.0],   # hsa-miR-3: low in both
    ]
    entropies = [0.2, 0.9, 0.5]
    return _make_data(mirna_ids, off_target_cell_types, expr, entropies)


# ── mirnas_for_off_target_cell_type ───────────────────────────────────────────

class TestMirnasForOffTargetCellType:

    def test_returns_dataframe(self, three_mirna_data):
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        result = mirnas_for_off_target_cell_type("Liver", mature_seqs, seed_seqs, df_expr, df_grouped)
        assert isinstance(result, pd.DataFrame)

    def test_expected_columns(self, three_mirna_data):
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        result = mirnas_for_off_target_cell_type("Liver", mature_seqs, seed_seqs, df_expr, df_grouped)
        assert set(result.columns) == {"MiRBase_ID", "mature_sequence", "seed", "mean_expr", "shannon_entropy"}

    def test_returns_all_passing_mirnas(self, three_mirna_data):
        """With threshold=0 all 3 miRNAs should be returned (up to top_n)."""
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        result = mirnas_for_off_target_cell_type("Liver", mature_seqs, seed_seqs, df_expr, df_grouped, threshold=0.0, top_n=10)
        assert len(result) == 3

    def test_threshold_filters_low_expression(self, three_mirna_data):
        """miR-3 has mean_expr=50 in Liver, should be excluded at threshold=100."""
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        result = mirnas_for_off_target_cell_type("Liver", mature_seqs, seed_seqs, df_expr, df_grouped, threshold=100.0)
        assert "hsa-miR-3" not in result["MiRBase_ID"].values

    def test_threshold_exact_boundary_included(self, three_mirna_data):
        """Threshold uses >=, so a miRNA at exactly the threshold is included."""
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        result = mirnas_for_off_target_cell_type("Liver", mature_seqs, seed_seqs, df_expr, df_grouped, threshold=50.0)
        assert "hsa-miR-3" in result["MiRBase_ID"].values

    def test_sorted_by_entropy_ascending(self, three_mirna_data):
        """Result must be sorted by shannon_entropy ascending."""
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        result = mirnas_for_off_target_cell_type("Liver", mature_seqs, seed_seqs, df_expr, df_grouped, threshold=0.0, top_n=10)
        assert result["shannon_entropy"].is_monotonic_increasing

    def test_top_n_limits_result_size(self, three_mirna_data):
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        result = mirnas_for_off_target_cell_type("Liver", mature_seqs, seed_seqs, df_expr, df_grouped, threshold=0.0, top_n=2)
        assert len(result) == 2

    def test_top_n_picks_lowest_entropy(self, three_mirna_data):
        """top_n=1 should return the miRNA with the lowest entropy."""
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        result = mirnas_for_off_target_cell_type("Liver", mature_seqs, seed_seqs, df_expr, df_grouped, threshold=0.0, top_n=1)
        assert result.iloc[0]["MiRBase_ID"] == "hsa-miR-1"  # entropy 0.2 is lowest

    def test_mature_sequence_populated(self, three_mirna_data):
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        result = mirnas_for_off_target_cell_type("Liver", mature_seqs, seed_seqs, df_expr, df_grouped)
        for _, row in result.iterrows():
            assert row["mature_sequence"] == mature_seqs[row["MiRBase_ID"]]

    def test_seed_populated(self, three_mirna_data):
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        result = mirnas_for_off_target_cell_type("Liver", mature_seqs, seed_seqs, df_expr, df_grouped)
        for _, row in result.iterrows():
            assert row["seed"] == seed_seqs[row["MiRBase_ID"]]

    def test_mirnas_missing_from_sequence_lookup_are_excluded(self, three_mirna_data):
        """miRNAs absent from mature_seqs must be dropped, not returned with empty sequences."""
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        sparse_seqs = {}  # nothing in lookup
        result = mirnas_for_off_target_cell_type("Liver", sparse_seqs, seed_seqs, df_expr, df_grouped)
        assert result.empty

    def test_partial_sequence_lookup_excludes_missing_only(self, three_mirna_data):
        """Only miRNAs that have a known sequence are returned."""
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        partial_seqs = {"hsa-miR-1": "UGGAGUGUGACAAUGGUGUUUG"}  # only miR-1 known
        result = mirnas_for_off_target_cell_type("Liver", partial_seqs, seed_seqs, df_expr, df_grouped)
        assert list(result["MiRBase_ID"]) == ["hsa-miR-1"]
        assert (result["mature_sequence"] != "").all()

    def test_invalid_off_target_cell_type_raises_value_error(self, three_mirna_data):
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        with pytest.raises(ValueError, match="Kidney"):
            mirnas_for_off_target_cell_type("Kidney", mature_seqs, seed_seqs, df_expr, df_grouped)

    def test_error_message_lists_available_off_target_cell_types(self, three_mirna_data):
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        with pytest.raises(ValueError) as exc_info:
            mirnas_for_off_target_cell_type("Kidney", mature_seqs, seed_seqs, df_expr, df_grouped)
        msg = str(exc_info.value)
        assert "Liver" in msg
        assert "Brain" in msg

    def test_empty_result_when_nothing_passes_threshold(self, three_mirna_data):
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        result = mirnas_for_off_target_cell_type("Liver", mature_seqs, seed_seqs, df_expr, df_grouped, threshold=9999.0)
        assert result.empty
        assert list(result.columns) == ["MiRBase_ID", "mature_sequence", "seed", "mean_expr", "shannon_entropy"]

    def test_mean_expr_values_match_off_target_cell_type_column(self, three_mirna_data):
        """mean_expr in result should equal the df_grouped value for the off-target cell type."""
        mature_seqs, seed_seqs, df_expr, df_grouped = three_mirna_data
        result = mirnas_for_off_target_cell_type("Brain", mature_seqs, seed_seqs, df_expr, df_grouped, threshold=0.0, top_n=10)
        for _, row in result.iterrows():
            expected = df_grouped.loc[row["MiRBase_ID"], "Brain"]
            assert row["mean_expr"] == pytest.approx(expected)

    def test_single_mirna_dataset(self):
        """Edge case: only one miRNA in the data."""
        mature_seqs = {"hsa-miR-X": "AUGCAUGC"}
        seed_seqs = {"hsa-miR-X": "AUGC"}
        df_grouped = pd.DataFrame(
            {"Liver": [300.0], "shannon_entropy": [0.1]},
            index=["hsa-miR-X"],
        )
        df_expr = pd.DataFrame({"Liver": [300.0]}, index=["hsa-miR-X"])
        result = mirnas_for_off_target_cell_type("Liver", mature_seqs, seed_seqs, df_expr, df_grouped)
        assert len(result) == 1
        assert result.iloc[0]["MiRBase_ID"] == "hsa-miR-X"
