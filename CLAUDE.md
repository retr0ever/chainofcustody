# Chain of Custody

mRNA sequence design and evaluation pipeline for the Serova x Berlin Biohack challenge.

## Project structure

```
chainofcustody/
  cli.py              # Single CLI entry point (chainofcustody command, optimize only)
  sequence.py         # mRNASequence dataclass + KOZAK constant
  cds/                # Gene fetching from Ensembl (get_canonical_cds)
  evaluation/          # 4-metric scoring pipeline
    structure.py       # Metric 1: ViennaRNA folding (5'UTR accessibility, global MFE)
    manufacturing.py   # Metric 2: GC windows, homopolymers, restriction sites
    stability.py       # Metric 3: GC3, MFE/nt, AU-rich elements
    ribonn.py          # Metric 4: RiboNN translation efficiency (optional, subprocess)
    scoring.py         # Pipeline orchestrator: runs all metrics, builds report + summary
    fitness.py         # Normalisation (all metrics to 0-1), weighted scoring, suggestion engine
    report.py          # Rich terminal output + markdown/JSON formatting
    utils.py           # Shared utilities (reverse_complement)
  optimization/        # pymoo NSGA-III multi-objective optimizer
    problem.py         # SequenceProblem: 4 objectives = 1 - normalised metric score each
    operators.py       # NucleotideSampling, NucleotideMutation (random, not protein-preserving)
    algorithm.py       # NSGA-III setup and run()
  three_prime/         # 3'UTR generation (miRNA database, cell-type seed maps)
data/
  MOESM3_ESM.xlsx      # Translation efficiency dataset (human, multiple cell types)
scripts/               # One-time data preparation utilities (not part of the package)
  merge_db.py          # Builds three_prime/db/ CSVs from Bioconductor microRNAome (requires R + rpy2)
vendor/
  RiboNN/              # git submodule — Sanofi RiboNN deep-CNN (translation efficiency)
tests/                 # pytest suite (skip integration tests by default)
```

## Critical conventions

### RNA throughout

All internal sequence handling uses **RNA nucleotides (A/U/G/C)**, not DNA. Callers must normalise DNA to RNA (replace T→U) before constructing an `mRNASequence`. Everything downstream operates in RNA.

- Start codon: `AUG` (not ATG)
- Stop codons: `UAA`, `UAG`, `UGA`
- ViennaRNA receives RNA directly (no conversion needed)
- The optimizer encodes nucleotides as integers: 0=A, 1=C, 2=G, 3=U
- The Ensembl API and FASTA files return DNA — convert at the call site

Do **not** add T->U or U->T conversion hacks in individual modules. If something needs DNA, convert at the boundary, not inside the pipeline.

### Metric display

The summary table shows **normalised 0-1 scores** (higher = better) for all metrics except 5'UTR accessibility (which shows raw MFE in kcal/mol). Raw values appear in the hint column. This keeps output comparable at a glance.

### Fitness scoring

4 metrics, weighted sum, all normalised to 0-1:

| Metric | Weight | GREEN threshold | What it measures |
|---|---|---|---|
| utr5_accessibility | 22% | MFE < -30 kcal/mol | 5'UTR structure stability |
| manufacturability | 30% | 0 violations | DNA synthesis feasibility |
| stability | 35% | Combined >= 0.7 | mRNA half-life (GC3, MFE/nt, AREs) |
| translation_efficiency | 13% | mean TE >= 1.5 | RiboNN predicted translation efficiency (optional) |

### Tests

```bash
uv run pytest              # unit tests only (default: skips integration)
uv run pytest -m integration  # integration tests (requires network)
uv run pytest -x -v         # stop on first failure, verbose
```

**Running tests:** All tests must pass before pushing. Run `uv run pytest` from the repo root.

**Mocking rules:**
- Mock `score_parsed` (not `compute_fitness`) in CLI tests — `compute_fitness` is pure arithmetic and should run on mock report data to catch display/formatting bugs.
- Patch imports at their **usage** location (e.g. `chainofcustody.cli.run`, not `chainofcustody.optimization.run`) since the CLI uses top-level imports.
- The mock report dict must include all top-level keys: `sequence_info`, `structure_scores`, `manufacturing_scores`, `stability_scores`, `ribonn_scores`, and `summary`. `ribonn_scores` must contain `mean_te` (float), `status`, and `message`.


**Test structure:**
- `tests/test_cli.py` — CLI integration tests (invoke via CliRunner, mock scoring pipeline)
- `tests/optimization/test_optimization.py` — problem dimensions, operator shapes, end-to-end NSGA-III
- `tests/cds/test_integration.py` — live Ensembl fetch (marked `integration`, skipped by default)
- `tests/*/test_dummy.py` — placeholder tests for modules under development

**Writing new tests:**
- Evaluation metric tests need ViennaRNA installed (structure.py, stability.py). Mock `RNA.fold` if testing without it.
- Sequences in tests should use RNA (`AUG`, not `ATG`).
- The optimizer uses `N_OBJECTIVES = len(METRIC_NAMES)` (currently 4). Never hardcode objective counts — import `N_OBJECTIVES` or `METRIC_NAMES` from `optimization.problem`.
- Use `pytest-mock`'s `mocker` fixture for patching (already a dev dependency).

## Python API

```python
from chainofcustody.optimization import mRNASequence, score_parsed, assemble_mrna
from chainofcustody.evaluation import compute_fitness

# Assemble and score a sequence (Kozak is inserted automatically by assemble_mrna)
full_seq = assemble_mrna(utr5="AAAAAAA", cds="AUGCCCAAAUGGG...", utr3="CCCGGG")

# Construct an mRNASequence and score it
seq = mRNASequence(utr5="AAAAAAA", cds="AUGCCCAAAUGGG...", utr3="CCCGGG")
report = score_parsed(seq)
fitness = compute_fitness(report)
```

## Incomplete modules

- `five_prime/` — deleted placeholder, no 5'UTR generator yet
- `optimization/operators.py` — random nucleotide mutation, does **not** preserve protein sequence

## Dependencies requiring system install

- **ViennaRNA** (`viennarna`): RNA secondary structure folding. Required for metrics 1 (structure) and 3 (stability). Install via conda or system package manager if pip install fails.
- **RiboNN**: Translation efficiency prediction. Bundled as a git submodule at `vendor/RiboNN`. Initialise with `git submodule update --init vendor/RiboNN`. Requires PyTorch (already in the venv).
