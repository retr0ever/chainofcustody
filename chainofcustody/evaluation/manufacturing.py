"""Metric 5: Manufacturability checks for DNA synthesis."""

import re

from .parser import ParsedSequence

# Common restriction enzyme recognition sites to avoid
RESTRICTION_SITES = {
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
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

    violations = []
    for i in range(len(seq) - window + 1):
        w = seq[i:i + window]
        gc = sum(1 for nt in w if nt in "GC") / window
        if gc < min_gc or gc > max_gc:
            violations.append({
                "position": i,
                "gc_content": round(gc * 100, 1),
                "issue": "too_low" if gc < min_gc else "too_high",
            })

    # Deduplicate: only report the worst violation per 50bp region
    deduped = []
    last_pos = -window
    for v in violations:
        if v["position"] - last_pos >= window:
            deduped.append(v)
            last_pos = v["position"]

    return {
        "pass": len(deduped) == 0,
        "violations": deduped,
        "windows_checked": len(seq) - window + 1,
    }


def check_homopolymers(seq: str, max_run: int = 8) -> dict:
    """Check for homopolymer runs (e.g. AAAAAAAAA) exceeding max_run length."""
    violations = []
    for nt in "ATGC":
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
        comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
        rc = "".join(comp[nt] for nt in reversed(recognition_seq))
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


def score_manufacturing(parsed: ParsedSequence) -> dict:
    """Run all manufacturability checks."""
    seq = parsed.raw

    gc = check_gc_windows(seq)
    homopolymers = check_homopolymers(seq)
    restriction = check_restriction_sites(seq)

    total_violations = len(gc["violations"]) + len(homopolymers["violations"]) + len(restriction["violations"])

    return {
        "gc_windows": gc,
        "homopolymers": homopolymers,
        "restriction_sites": restriction,
        "total_violations": total_violations,
        "overall_pass": total_violations == 0,
    }
