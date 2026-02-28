"""mRNA sequence dataclass and Kozak consensus constant."""

from dataclasses import dataclass

# Kozak consensus sequence inserted between the 5'UTR and the CDS start codon.
KOZAK = "GCCACC"

# 5' cap analogue (CleanCap-AG / ARCA): contributes a GGG trinucleotide at the
# transcript 5' end in sequence notation.  The actual structure is a
# m7GpppG(2'-O-methyl) cap; GGG is the canonical sequence representation used
# in mRNA therapeutic design.
CAP5 = "GGG"

# Poly-A tail length standard for human mRNA therapeutics (120 nt).
POLY_A_LENGTH = 120
POLY_A = "A" * POLY_A_LENGTH


@dataclass
class mRNASequence:
    """An mRNA sequence split into its three functional regions.

    ``__str__`` / ``__len__`` return only the core transcript
    (5'UTR + CDS + 3'UTR), which is the region used for structure
    prediction and all scoring metrics.

    Use :attr:`full_sequence` to obtain the complete molecule
    including the 5' cap and poly-A tail.
    """
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

    @property
    def full_sequence(self) -> str:
        """Complete mRNA molecule: 5' cap + core transcript + poly-A tail."""
        return CAP5 + self.utr5 + self.cds + self.utr3 + POLY_A

    @property
    def full_length(self) -> int:
        """Total length including 5' cap and poly-A tail."""
        return len(CAP5) + len(self) + POLY_A_LENGTH

    def __len__(self) -> int:
        return len(self.utr5) + len(self.cds) + len(self.utr3)

    def __repr__(self) -> str:
        return (
            f"mRNASequence("
            f"cap={len(CAP5)}nt + "
            f"utr5={len(self.utr5)}nt + "
            f"cds={len(self.cds)}nt + "
            f"utr3={len(self.utr3)}nt + "
            f"polyA={POLY_A_LENGTH}nt, "
            f"total={self.full_length}nt)"
        )

    def __str__(self) -> str:
        return self.utr5 + self.cds + self.utr3
