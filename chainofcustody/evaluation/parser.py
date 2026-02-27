"""Parse mRNA sequences into 5'UTR / CDS / 3'UTR regions."""

from dataclasses import dataclass
from pathlib import Path
import re

STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass
class ParsedSequence:
    """An mRNA sequence split into its functional regions."""
    raw: str
    utr5: str
    cds: str
    utr3: str

    @property
    def codons(self) -> list[str]:
        """Extract codons from the CDS."""
        return [self.cds[i:i+3] for i in range(0, len(self.cds), 3)]

    @property
    def cds_start(self) -> int:
        return len(self.utr5)

    @property
    def cds_end(self) -> int:
        return len(self.utr5) + len(self.cds)


def load_fasta(path: str | Path) -> str:
    """Load a sequence from a FASTA file, stripping headers and whitespace."""
    lines = Path(path).read_text().strip().splitlines()
    seq_lines = [line.strip() for line in lines if not line.startswith(">")]
    return "".join(seq_lines).upper()


def clean_sequence(seq: str) -> str:
    """Normalise a sequence: uppercase, strip whitespace, convert U→T."""
    seq = re.sub(r"\s+", "", seq).upper()
    seq = seq.replace("U", "T")
    if not re.match(r"^[ATGC]+$", seq):
        invalid = set(seq) - set("ATGC")
        raise ValueError(f"Sequence contains invalid characters: {invalid}")
    return seq


def find_cds(seq: str) -> tuple[int, int]:
    """Find CDS boundaries: first ATG to first in-frame stop codon."""
    atg_pos = seq.find("ATG")
    if atg_pos == -1:
        raise ValueError("No ATG start codon found in sequence")

    # Walk in-frame from the ATG looking for a stop codon
    for i in range(atg_pos, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in STOP_CODONS:
            # CDS includes the stop codon
            return atg_pos, i + 3

    raise ValueError("No in-frame stop codon found after ATG")


def parse_sequence(
    seq: str,
    utr5_end: int | None = None,
    cds_end: int | None = None,
) -> ParsedSequence:
    """
    Parse an mRNA sequence into regions.

    Args:
        seq: Raw DNA/RNA sequence string or path to FASTA file.
        utr5_end: Manual 5'UTR end position (0-indexed). If None, auto-detect from first ATG.
        cds_end: Manual CDS end position (0-indexed, exclusive). If None, auto-detect from stop codon.
    """
    # Handle file input — only try path check for short strings that look like filenames
    # (no whitespace, no DNA-only content, reasonable length for a path)
    is_path = (
        len(seq) < 300
        and "\n" not in seq
        and not re.match(r"^[ATGCU\s]+$", seq, re.IGNORECASE)
    )
    if is_path:
        try:
            p = Path(seq)
            if p.exists():
                seq = load_fasta(seq)
            else:
                seq = clean_sequence(seq)
        except OSError:
            seq = clean_sequence(seq)
    else:
        seq = clean_sequence(seq)

    if utr5_end is not None and cds_end is not None:
        # Manual boundaries
        cds_start = utr5_end
        cds_stop = cds_end
    else:
        # Auto-detect from ATG + stop codon
        cds_start, cds_stop = find_cds(seq)

    utr5 = seq[:cds_start]
    cds = seq[cds_start:cds_stop]
    utr3 = seq[cds_stop:]

    # Validate CDS
    if len(cds) % 3 != 0:
        raise ValueError(f"CDS length ({len(cds)}) is not divisible by 3")
    if not cds.startswith("ATG"):
        raise ValueError(f"CDS does not start with ATG: starts with {cds[:3]}")

    return ParsedSequence(raw=seq, utr5=utr5, cds=cds, utr3=utr3)
