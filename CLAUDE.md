# Chain of Custody

mRNA sequence design and evaluation pipeline for the Serova x Berlin Biohack challenge.

## Project structure

```
chainofcustody/
  cli.py              # Single CLI entry point (chainofcustody command, optimize only)
  cds/                # Gene fetching from Ensembl (get_canonical_cds)
  evaluation/          # 6-metric scoring pipeline
    parser.py          # Sequence parsing: FASTA loading, T->U normalisation, CDS boundary detection
    codons.py          # Metric 1: CAI, GC content, rare codon clusters, liver/target selectivity
    mirna.py           # Metric 3: miR-122 seed site scanning + accidental match warnings
    structure.py       # Metric 4: ViennaRNA folding (5'UTR accessibility, global MFE)
    manufacturing.py   # Metric 5: GC windows, homopolymers, restriction sites
    stability.py       # Metric 6: GC3, MFE/nt, AU-rich elements
    fitness.py         # Normalisation (all metrics to 0-1), weighted scoring, suggestion engine
    report.py          # Rich terminal output + markdown/JSON formatting
    data.py            # MOESM3_ESM.xlsx loader, per-codon TE weight computation
    mutations.py       # Protein-preserving operators: synonymous swap, crossover, miR-122 insertion
    evolve.py          # Evolutionary loop (selection, crossover, mutation, early stopping)
  optimization/        # pymoo NSGA-III multi-objective optimizer
    problem.py         # SequenceProblem: 6 objectives = 1 - normalised metric score each
    operators.py       # NucleotideSampling, NucleotideMutation (random, not protein-preserving)
    algorithm.py       # NSGA-III setup and run()
  three_prime/         # 3'UTR generation (miRNA database, cell-type seed maps)
data/
  MOESM3_ESM.xlsx      # Translation efficiency dataset (human, multiple cell types)
scripts/               # One-time data preparation utilities (not part of the package)
  merge_db.py          # Builds three_prime/db/ CSVs from Bioconductor microRNAome (requires R + rpy2)
tests/                 # pytest suite (skip integration tests by default)
```

## Critical conventions

### RNA throughout

All internal sequence handling uses **RNA nucleotides (A/U/G/C)**, not DNA. The only T->U conversion point is `parser.py:clean_sequence()`, which normalises any input (DNA or RNA) to RNA on entry. Everything downstream operates in RNA.

- Start codon: `AUG` (not ATG)
- Stop codons: `UAA`, `UAG`, `UGA`
- ViennaRNA receives RNA directly (no conversion needed)
- The optimizer encodes nucleotides as integers: 0=A, 1=C, 2=G, 3=U
- The Ensembl API and FASTA files return DNA — `clean_sequence()` handles conversion

Do **not** add T->U or U->T conversion hacks in individual modules. If something needs DNA, convert at the boundary, not inside the pipeline.

### Metric display

The summary table shows **normalised 0-1 scores** (higher = better) for all metrics except 5'UTR accessibility (which shows raw MFE in kcal/mol). Raw values appear in the hint column. This keeps output comparable at a glance.

### Fitness scoring

6 metrics, weighted sum, all normalised to 0-1:

| Metric | Weight | GREEN threshold | What it measures |
|---|---|---|---|
| codon_quality | 15% | CAI >= 0.8 | Human codon adaptation |
| gc_content | 10% | CDS 40-60% | Nucleotide composition balance |
| mir122_detargeting | 25% | 3+ sites in 3'UTR | Liver-specific silencing |
| utr5_accessibility | 10% | MFE < -30 kcal/mol | 5'UTR structure stability |
| manufacturability | 15% | 0 violations | DNA synthesis feasibility |
| stability | 25% | Combined >= 0.7 | mRNA half-life (GC3, MFE/nt, AREs) |

### Tests

```bash
uv run pytest              # unit tests only (default: skips integration)
uv run pytest -m integration  # integration tests (requires network)
uv run pytest -x -v         # stop on first failure, verbose
```

**Running tests:** All tests must pass before pushing. Run `uv run pytest` from the repo root.

**Mocking rules:**
- Mock `score_sequence` (not `compute_fitness`) in CLI tests — `compute_fitness` is pure arithmetic and should run on mock report data to catch display/formatting bugs.
- Patch imports at their **usage** location (e.g. `chainofcustody.cli.run`, not `chainofcustody.optimization.run`) since the CLI uses top-level imports.
- The mock report dict must include all 6 top-level keys: `sequence_info`, `codon_scores`, `mirna_scores`, `structure_scores`, `manufacturing_scores`, `stability_scores`, and `summary`. Missing keys will cause KeyError in `compute_fitness` or `print_report`.

**Test structure:**
- `tests/test_cli.py` — CLI integration tests (invoke via CliRunner, mock scoring pipeline)
- `tests/optimization/test_optimization.py` — problem dimensions, operator shapes, end-to-end NSGA-III
- `tests/cds/test_integration.py` — live Ensembl fetch (marked `integration`, skipped by default)
- `tests/*/test_dummy.py` — placeholder tests for modules under development

**Writing new tests:**
- Evaluation metric tests need ViennaRNA installed (structure.py, stability.py). Mock `RNA.fold` if testing without it.
- Sequences in tests should use RNA (`AUG`, not `ATG`) unless testing the T->U conversion in `clean_sequence` itself.
- The optimizer uses `N_OBJECTIVES = len(METRIC_NAMES)` (currently 6). Never hardcode objective counts — import `N_OBJECTIVES` or `METRIC_NAMES` from `optimization.problem`.
- Use `pytest-mock`'s `mocker` fixture for patching (already a dev dependency).

## Python API

```python
from chainofcustody.evaluation import score_sequence, compute_fitness, evaluate_candidate, evaluate_batch, evolve

# Score a single sequence (accepts DNA or RNA)
report = score_sequence("AUGCCCAAAGGG...")
fitness = compute_fitness(report)

# Or use the convenience wrapper
result = evaluate_candidate("AUGCCCAAAGGG...", label="my_seq")

# Batch scoring with ranking
results = evaluate_batch(["seq1", "seq2", "seq3"])

# Evolutionary optimisation (protein-preserving mutations)
result = evolve(population=["AUGCCC..."], generations=20, mutation_rate=0.05)
```

## Incomplete modules

- `five_prime/` — deleted placeholder, no 5'UTR generator yet
- `optimization/operators.py` — random nucleotide mutation, does **not** preserve protein sequence (unlike `evaluation/mutations.py` which does)
- No concatenation step (5'UTR + CDS + 3'UTR assembly)

## Dependencies requiring system install

- **ViennaRNA** (`viennarna`): RNA secondary structure folding. Required for metrics 4 (structure) and 6 (stability). Install via conda or system package manager if pip install fails.
