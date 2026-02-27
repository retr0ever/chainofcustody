# Chain of Custody

A pipeline for RNA sequence construction and evaluation.

## Installation

This project uses [uv](https://docs.astral.sh/uv/) for dependency management.

**Install uv** (macOS/Linux):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Or via Homebrew:

```bash
brew install uv
```

**Install dependencies and run tests:**

```bash
uv sync
uv run pytest
```

## Pipeline Overview

```
Gene Input
    │
    ▼
┌─────────┐
│ initial │  Transcribes a gene into an RNA sequence
└─────────┘
    │
    ├──────────────────────┐
    ▼                      ▼
┌───────────┐        ┌────────────┐
│ five_prime│        │ three_prime│
│  (5' UTR) │        │  (3' UTR)  │
└───────────┘        └────────────┘
    │                      │
    └──────────┬───────────┘
               ▼
    ┌─────────────────────┐
    │  Concatenation      │
    │  5' + RNA + 3'      │
    └─────────────────────┘
               │
               ▼
        ┌────────────┐
        │ evaluation │  Evaluates the final construct
        └────────────┘
```

## Modules

### `initial`
Takes a gene as input and produces an RNA sequence (transcription).

### `five_prime`
Generates a 5' UTR sequence to be prepended to the RNA.

### `three_prime`
Generates a 3' UTR sequence to be appended to the RNA.

### `evaluation`
Receives the concatenated construct (`5' UTR + RNA + 3' UTR`) and evaluates it.
