"""Tests for manufacturing checks, focusing on uORF detection."""
import pytest

from chainofcustody.sequence import mRNASequence
from chainofcustody.evaluation.manufacturing import (
    check_uorfs,
    score_manufacturing,
)


# ── check_uorfs ──────────────────────────────────────────────────────────────

def test_uorf_no_aug():
    result = check_uorfs("AAACCCUUUGGG")
    assert result["pass"] is True
    assert result["count"] == 0
    assert result["positions"] == []


def test_uorf_single_aug():
    result = check_uorfs("AAAUGCCC")
    assert result["pass"] is False
    assert result["count"] == 1
    assert result["positions"] == [2]


def test_uorf_multiple_augs():
    result = check_uorfs("AUGCCCAUGAAA")
    assert result["count"] == 2
    assert 0 in result["positions"]
    assert 6 in result["positions"]


def test_uorf_dna_input_normalised():
    """T in input should be treated as U (DNA → RNA normalisation)."""
    result = check_uorfs("AATGCCC")   # ATG in DNA = AUG in RNA
    assert result["count"] == 1


def test_uorf_empty_sequence():
    result = check_uorfs("")
    assert result["pass"] is True
    assert result["count"] == 0


def test_uorf_violations_list_has_position_and_sequence():
    result = check_uorfs("CCAUGGG")
    assert len(result["violations"]) == 1
    v = result["violations"][0]
    assert v["position"] == 2
    assert v["sequence"] == "AUG"


# ── score_manufacturing integration ──────────────────────────────────────────

def test_score_manufacturing_includes_uorfs_key():
    seq = mRNASequence(utr5="AAACCCUUU", cds="AUGAAGUAA", utr3="")
    result = score_manufacturing(seq)
    assert "uorfs" in result


def test_score_manufacturing_uorf_counted_in_utr5_violations():
    """An AUG in the 5'UTR must increase utr5_violations."""
    seq_clean = mRNASequence(utr5="AAACCCUUU", cds="AUGAAGUAA", utr3="")
    seq_uorf  = mRNASequence(utr5="AAAUGCCC",  cds="AUGAAGUAA", utr3="")

    clean = score_manufacturing(seq_clean)
    uorf  = score_manufacturing(seq_uorf)

    assert clean["utr5_violations"] < uorf["utr5_violations"]
    assert uorf["uorfs"]["count"] == 1


def test_score_manufacturing_uorf_not_counted_in_total_violations():
    """uORF count must NOT inflate total_violations (which covers the full sequence)."""
    seq = mRNASequence(utr5="AAAUGCCC", cds="AUGAAGUAA", utr3="")
    result = score_manufacturing(seq)
    # total_violations covers gc_windows + homopolymers + restriction_sites on the
    # full sequence; uORFs are 5'UTR-only and tracked separately.
    total_from_standard = (
        len(result["gc_windows"]["violations"])
        + len(result["homopolymers"]["violations"])
        + len(result["restriction_sites"]["violations"])
    )
    assert result["total_violations"] == total_from_standard


def test_score_manufacturing_clean_utr5_zero_uorf_violations():
    seq = mRNASequence(utr5="CCCCCCUUUUUU", cds="AUGAAGUAA", utr3="")
    result = score_manufacturing(seq)
    assert result["uorfs"]["count"] == 0
    assert result["uorfs"]["pass"] is True
