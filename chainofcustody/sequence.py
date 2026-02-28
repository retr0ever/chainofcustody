"""mRNA sequence dataclass and Kozak consensus constant."""

from dataclasses import dataclass

# Kozak consensus sequence inserted between the 5'UTR and the CDS start codon.
KOZAK = "GCCACC"


@dataclass
class mRNASequence:
    """An mRNA sequence split into its three functional regions."""
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

    def __len__(self) -> int:
        return len(self.utr5) + len(self.cds) + len(self.utr3)

    def __repr__(self) -> str:
        return (
            f"mRNASequence("
            f"utr5={len(self.utr5)}nt, "
            f"cds={len(self.cds)}nt, "
            f"utr3={len(self.utr3)}nt, "
            f"total={len(self)}nt)"
        )

    def __str__(self) -> str:
        return self.utr5 + self.cds + self.utr3
