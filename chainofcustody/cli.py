import csv
import json
from pathlib import Path

import numpy as np
import rich_click as click
from rich.console import Console

from chainofcustody.evaluation.fitness import compute_fitness
from chainofcustody.evaluation.report import print_batch_report, print_report, score_sequence
from chainofcustody.optimization import METRIC_NAMES, SequenceProblem, run

console = Console()

_CSV_COLUMNS = ["rank", "label", "sequence", *METRIC_NAMES, "overall"]


def _write_csv(path: Path, results: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_CSV_COLUMNS)
        writer.writeheader()
        for rank, r in enumerate(results, start=1):
            scores = r["fitness"]["scores"]
            writer.writerow({
                "rank": rank,
                "label": r["label"],
                "sequence": r["sequence"],
                **{m: round(scores[m]["value"], 4) for m in METRIC_NAMES},
                "overall": r["fitness"]["overall"],
            })


@click.command()
@click.option("--seq-len", type=int, default=100, show_default=True, help="Length of each candidate sequence.")
@click.option("--pop-size", type=int, default=128, show_default=True, help="Population size (must be ≥ number of reference directions; default covers Das-Dennis n_partitions=4 → 126 directions).")
@click.option("--n-gen", type=int, default=50, show_default=True, help="Number of generations.")
@click.option("--mutation-rate", type=float, default=0.01, show_default=True, help="Per-position mutation probability.")
@click.option("--seed", type=int, default=None, help="Random seed for reproducibility.")
@click.option("--workers", type=int, default=None, help="Parallel worker processes for fitness evaluation (default: all CPU cores). Pass 1 to disable parallelism.")
@click.option("--output", "output_fmt", type=click.Choice(["summary", "json"]), default="summary", show_default=True, help="Output format.")
@click.option("--csv", "csv_path", type=click.Path(dir_okay=False, writable=True, path_type=Path), default=None, help="Write Pareto-front results to a CSV file.")
def main(seq_len: int, pop_size: int, n_gen: int, mutation_rate: float, seed: int | None, workers: int | None, output_fmt: str, csv_path: Path | None) -> None:
    """Run NSGA3 to evolve an optimal nucleotide sequence population."""
    console.print(
        f"\nRunning NSGA3 — seq_len=[bold]{seq_len}[/bold]  "
        f"pop_size=[bold]{pop_size}[/bold]  "
        f"n_gen=[bold]{n_gen}[/bold]  "
        f"mutation_rate=[bold]{mutation_rate}[/bold]  "
        f"workers=[bold]{workers if workers is not None else 'auto'}[/bold]\n"
    )

    X, F = run(seq_len=seq_len, pop_size=pop_size, n_gen=n_gen, mutation_rate=mutation_rate, seed=seed, n_workers=workers)

    problem = SequenceProblem(seq_len=seq_len)
    sequences = problem.decode(X)

    results = []
    for i, seq in enumerate(sequences):
        try:
            report = score_sequence(seq)
            fitness = compute_fitness(report)
            results.append({"label": f"pareto_{i + 1}", "sequence": seq, "report": report, "fitness": fitness})
        except Exception:
            continue

    if not results:
        console.print("[bold red]Error:[/bold red] No sequences could be scored.")
        raise SystemExit(1)

    results.sort(key=lambda r: r["fitness"]["overall"], reverse=True)

    if csv_path:
        _write_csv(csv_path, results)
        console.print(f"Results written to [bold]{csv_path}[/bold]\n")

    if output_fmt == "json":
        out = []
        for r in results:
            data = dict(r["report"])
            data["fitness"] = r["fitness"]
            data["label"] = r["label"]
            out.append(data)
        console.print_json(json.dumps(out, indent=2, default=str))
    else:
        console.print(f"Pareto front: [bold]{len(results)}[/bold] scored sequences\n")
        if len(results) == 1:
            print_report(console, results[0]["report"], label=results[0]["label"])
        else:
            print_batch_report(console, results)
