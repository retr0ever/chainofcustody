import rich_click as click
from rich.console import Console

console = Console()


@click.command()
@click.option("--seq-len", type=int, default=100, show_default=True, help="Length of each candidate sequence.")
@click.option("--pop-size", type=int, default=100, show_default=True, help="Population size.")
@click.option("--n-gen", type=int, default=50, show_default=True, help="Number of generations.")
@click.option("--mutation-rate", type=float, default=0.01, show_default=True, help="Per-position mutation probability.")
@click.option("--seed", type=int, default=None, help="Random seed for reproducibility.")
@click.option("--workers", type=int, default=None, help="Parallel worker processes for fitness evaluation (default: all CPU cores). Pass 1 to disable parallelism.")
@click.option("--output", "output_fmt", type=click.Choice(["summary", "json"]), default="summary", show_default=True, help="Output format.")
def main(seq_len: int, pop_size: int, n_gen: int, mutation_rate: float, seed: int | None, workers: int | None, output_fmt: str) -> None:
    """Run NSGA3 to evolve an optimal nucleotide sequence population."""
    import json
    import numpy as np
    from chainofcustody.optimization import run, METRIC_NAMES, SequenceProblem
    from chainofcustody.evaluation.report import score_sequence, print_report, print_batch_report
    from chainofcustody.evaluation.fitness import compute_fitness

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
            results.append({"label": f"pareto_{i + 1}", "report": report, "fitness": fitness})
        except Exception:
            continue

    if not results:
        console.print("[bold red]Error:[/bold red] No sequences could be scored.")
        raise SystemExit(1)

    results.sort(key=lambda r: r["fitness"]["overall"], reverse=True)

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
