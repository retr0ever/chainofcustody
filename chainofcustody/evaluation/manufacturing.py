"""Metric 5: Manufacturability checks for DNA synthesis."""

import re

import numpy as np

from chainofcustody.sequence import mRNASequence
from .utils import reverse_complement

# Common restriction enzyme recognition sites to avoid
RESTRICTION_SITES = {
    "BsaI": "GGUCUC",
    "BsmBI": "CGUCUC",
    "EcoRI": "GAAUUC",
    "BamHI": "GGAUCC",
    "HindIII": "AAGCUU",
    "NotI": "GCGGCCGC",
}


def check_gc_windows(seq: str, window: int = 50, min_gc: float = 0.30, max_gc: float = 0.70) -> dict:
    """
    Check GC content in sliding windows across the sequence.
    Returns violations where GC falls outside the acceptable range.
    """
    if len(seq) < window:
        gc = sum(1 for nt in seq if nt in "GC") / len(seq)
        return {
            "pass": min_gc <= gc <= max_gc,
            "violations": [],
            "windows_checked": 1,
        }

    # Vectorised sliding-window GC count via cumulative sum — O(n), no Python loop.
    is_gc = np.frombuffer(seq.encode(), dtype=np.uint8)
    is_gc = ((is_gc == ord("G")) | (is_gc == ord("C"))).astype(np.int32)
    cumsum = np.empty(len(is_gc) + 1, dtype=np.int32)
    cumsum[0] = 0
    np.cumsum(is_gc, out=cumsum[1:])
    window_gc_counts = cumsum[window:] - cumsum[:-window]
    gc_fractions = window_gc_counts / window

    n_windows = len(gc_fractions)
    violation_mask = (gc_fractions < min_gc) | (gc_fractions > max_gc)
    violation_indices = np.where(violation_mask)[0]

    violations = []
    last_pos = -window
    for i in violation_indices.tolist():
        if i - last_pos >= window:
            gc_val = float(gc_fractions[i])
            violations.append({
                "position": i,
                "gc_content": round(gc_val * 100, 1),
                "issue": "too_low" if gc_val < min_gc else "too_high",
            })
            last_pos = i

    return {
        "pass": len(violations) == 0,
        "violations": violations,
        "windows_checked": n_windows,
    }


def check_homopolymers(seq: str, max_run: int = 8) -> dict:
    """Check for homopolymer runs (e.g. AAAAAAAAA) exceeding max_run length."""
    violations = []
    for nt in "AUGC":
        pattern = f"{nt}{{{max_run + 1},}}"
        for m in re.finditer(pattern, seq):
            violations.append({
                "position": m.start(),
                "nucleotide": nt,
                "length": m.end() - m.start(),
                "sequence": m.group(),
            })

    return {
        "pass": len(violations) == 0,
        "violations": violations,
        "max_allowed": max_run,
    }


def check_restriction_sites(seq: str, sites: dict[str, str] | None = None) -> dict:
    """Scan for restriction enzyme recognition sites."""
    sites = sites or RESTRICTION_SITES
    violations = []

    for enzyme, recognition_seq in sites.items():
        # Check forward strand
        for m in re.finditer(re.escape(recognition_seq), seq):
            violations.append({
                "position": m.start(),
                "enzyme": enzyme,
                "sequence": recognition_seq,
                "strand": "forward",
            })
        # Check reverse complement
        rc = reverse_complement(recognition_seq)
        if rc != recognition_seq:  # skip palindromes (already found)
            for m in re.finditer(re.escape(rc), seq):
                violations.append({
                    "position": m.start(),
                    "enzyme": enzyme,
                    "sequence": rc,
                    "strand": "reverse",
                })

    return {
        "pass": len(violations) == 0,
        "violations": violations,
        "enzymes_checked": list(sites.keys()),
    }


def score_manufacturing(parsed: mRNASequence) -> dict:
    """Run all manufacturability checks.

    Violations are computed on the full assembled mRNA for reporting purposes.
    ``utr5_violations`` counts only violations within the 5'UTR — the region
    the optimiser actually controls — and is used by the fitness normaliser.
    """
    seq = str(parsed)

    gc = check_gc_windows(seq)
    homopolymers = check_homopolymers(seq)
    restriction = check_restriction_sites(seq)
    total_violations = (
        len(gc["violations"]) + len(homopolymers["violations"]) + len(restriction["violations"])
    )

    # 5'UTR-only violations (the evolved region the optimizer controls)
    utr5_seq = parsed.utr5
    if utr5_seq:
        gc_u = check_gc_windows(utr5_seq)
        hp_u = check_homopolymers(utr5_seq)
        rs_u = check_restriction_sites(utr5_seq)
        utr5_violations = len(gc_u["violations"]) + len(hp_u["violations"]) + len(rs_u["violations"])
    else:
        utr5_violations = 0

    return {
        "gc_windows": gc,
        "homopolymers": homopolymers,
        "restriction_sites": restriction,
        "total_violations": total_violations,
        "utr5_violations": utr5_violations,
        "overall_pass": total_violations == 0,
    }
