# Chain of Custody

A pipeline for mRNA sequence construction and evaluation. Design tissue-specific mRNA sequences and score them across 5 metrics.

## Installation

This project uses [uv](https://docs.astral.sh/uv/) for dependency management.

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
uv sync
```

## Usage

### Fetch a gene's CDS

```bash
uv run chainofcustody fetch POU5F1
```

### Score a single sequence

```bash
uv run chainofcustody evaluate sequence.fasta
uv run chainofcustody evaluate --gene POU5F1
uv run chainofcustody fetch POU5F1 | uv run chainofcustody evaluate
```

### Score and rank multiple candidates

```bash
uv run chainofcustody evaluate a.fasta b.fasta c.fasta
uv run chainofcustody evaluate candidates/
```

### Output formats

```bash
uv run chainofcustody evaluate seq.fasta                  # rich terminal (default)
uv run chainofcustody evaluate seq.fasta --output json     # JSON
uv run chainofcustody evaluate seq.fasta --output markdown # markdown
```

### Options

```
--gene TEXT        Gene symbol — fetches CDS and scores it directly
--target TEXT      Target cell type column from dataset
--offtarget TEXT   Off-target cell type (default: HepG2 / liver)
--utr5-end INT    Manual 5'UTR end position
--cds-end INT     Manual CDS end position
--output FORMAT   rich | json | markdown
```

## Pipeline Overview

```
Gene Input
    |
    v
+---------+
| initial |  Fetch canonical CDS from Ensembl
+---------+
    |
    +---------------------+
    v                     v
+-----------+       +------------+
| five_prime|       | three_prime|
|  (5' UTR) |       |  (3' UTR)  |
+-----------+       +------------+
    |                     |
    +---------+-----------+
              v
    +---------------------+
    |  Concatenation      |
    |  5'UTR + CDS + 3'UTR|
    +---------------------+
              |
              v
       +------------+
       | evaluation |  Score across 5 metrics
       +------------+
```

## Evaluation Metrics

| # | Metric | What it measures |
|---|--------|-----------------|
| 1 | Codon quality | CAI, GC content, rare codon clusters |
| 2 | Cell-type selectivity | Liver detargeting score, target selectivity ratio |
| 3 | miRNA detargeting | miR-122 sites in 3'UTR, accidental off-target sites |
| 4 | Structure | 5'UTR accessibility (MFE), global folding energy |
| 5 | Manufacturability | GC windows, homopolymers, restriction enzyme sites |

Each metric is normalised to 0-1 and weighted for an overall fitness score. Candidates are ranked by fitness when evaluating in batch.

## Python API

```python
from chainofcustody.evaluation import evaluate_candidate, evaluate_batch

# Single sequence
result = evaluate_candidate("ATGCCC...TAA", label="v3")
result["fitness"]["scores"]    # per-metric normalised scores
result["suggestions"]          # actionable improvement hints

# Batch — returns ranked list (best first)
results = evaluate_batch([
    {"seq": seq1, "label": "v1"},
    {"seq": seq2, "label": "v2"},
])
```

## Modules

### `initial`
Fetches the canonical CDS for a human gene symbol from Ensembl.

### `five_prime`
Generates a 5'UTR sequence to prepend to the CDS.

### `three_prime`
Generates a 3'UTR sequence to append to the CDS.

### `evaluation`
Scores the concatenated construct across 5 metrics and returns a report with traffic-light status, fitness scores, and improvement suggestions.

## Development

```bash
uv sync
uv run pytest
```
