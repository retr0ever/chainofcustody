"""Shared sequence utilities for the evaluation package."""


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of an RNA sequence."""
    comp = {"A": "U", "U": "A", "G": "C", "C": "G"}
    return "".join(comp[nt] for nt in reversed(seq))
