import csv
import json
from pathlib import Path

import rich_click as click
from rich.console import Console
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn

from chainofcustody.evaluation.fitness import compute_fitness
from chainofcustody.evaluation.report import print_batch_report, print_report
from chainofcustody.cds import GeneNotFoundError, get_canonical_cds
from chainofcustody.optimization import KOZAK, METRIC_NAMES, mRNASequence, SequenceProblem, assemble_mrna, run, score_parsed
from chainofcustody.three_prime import generate_utr3

console = Console()

_CSV_COLUMNS = ["generation", "sequence", *METRIC_NAMES, "overall"]


def _write_csv(path: Path, history: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(history)


def _to_rna(seq: str) -> str:
    """Convert a DNA sequence to RNA (T → U, uppercase)."""
    return seq.upper().replace("T", "U")


@click.command()
@click.option("--gene", default="POU5F1", show_default=True, help="HGNC gene symbol whose canonical CDS to optimise (e.g. TP53).")
@click.option("--off-target-cell-type", default="Dendritic_cell", show_default=True, help="Off-target cell type used to generate the 3'UTR sponge (must match a name in the expression database).")
@click.option("--utr5-min", type=int, default=10, show_default=True, help="Minimum 5'UTR length the GA can produce.")
@click.option("--utr5-max", type=int, default=100, show_default=True, help="Maximum 5'UTR length the GA can produce.")
@click.option("--pop-size", type=int, default=128, show_default=True, help="Population size (must be ≥ number of reference directions; default covers Das-Dennis n_partitions=4 → 126 directions).")
@click.option("--n-gen", type=int, default=50, show_default=True, help="Number of generations.")
@click.option("--mutation-rate", type=float, default=0.01, show_default=True, help="Per-position mutation probability.")
@click.option("--seed", type=int, default=None, help="Random seed for reproducibility.")
@click.option("--workers", type=int, default=None, help="Parallel worker processes for fitness evaluation (default: all CPU cores). Pass 1 to disable parallelism.")
@click.option("--output", "output_fmt", type=click.Choice(["summary", "json"]), default="summary", show_default=True, help="Output format.")
@click.option("--csv", "csv_path", type=click.Path(dir_okay=False, writable=True, path_type=Path), default=None, help="Write Pareto-front results to a CSV file.")
def main(gene: str, off_target_cell_type: str, utr5_min: int, utr5_max: int, pop_size: int, n_gen: int, mutation_rate: float, seed: int | None, workers: int | None, output_fmt: str, csv_path: Path | None) -> None:
    """Run NSGA3 to evolve an optimal 5'UTR for a given gene."""
    if utr5_min > utr5_max:
        console.print(f"[bold red]Error:[/bold red] --utr5-min ({utr5_min}) must be ≤ --utr5-max ({utr5_max}).")
        raise SystemExit(1)

    console.print(f"\nResolving canonical CDS for gene [bold]{gene}[/bold]…")
    try:
        cds = _to_rna(get_canonical_cds(gene))
    except GeneNotFoundError as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        raise SystemExit(1)

    console.print(f"Generating 3'UTR sponge for off-target cell type [bold]{off_target_cell_type}[/bold]…")
    try:
        utr3 = generate_utr3(off_target_cell_type)
    except (ValueError, RuntimeError) as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        raise SystemExit(1)

    console.print(
        f"CDS resolved: [bold]{len(cds)} nt[/bold]  "
        f"3'UTR: [bold]{len(utr3)} nt[/bold]  "
        f"5'UTR length: [bold]{utr5_min}–{utr5_max} nt[/bold] (evolved)\n"
    )
    console.print(
        f"Running NSGA3 — "
        f"pop_size=[bold]{pop_size}[/bold]  "
        f"n_gen=[bold]{n_gen}[/bold]  "
        f"mutation_rate=[bold]{mutation_rate}[/bold]  "
        f"workers=[bold]{workers if workers is not None else 'auto'}[/bold]\n"
    )

    with Progress(
        SpinnerColumn(),
        TextColumn("[bold blue]Evolving[/bold blue] {task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("5'UTR", total=n_gen)
        X, F, history = run(
            utr5_min=utr5_min, utr5_max=utr5_max, cds=cds, utr3=utr3,
            pop_size=pop_size, n_gen=n_gen,
            mutation_rate=mutation_rate, seed=seed, n_workers=workers,
            progress=progress, progress_task=task,
        )

    problem = SequenceProblem(utr5_min=utr5_min, utr5_max=utr5_max, cds=cds, utr3=utr3)
    sequences = problem.decode(X)

    results = []
    for i, seq in enumerate(sequences):
        utr5_len = int(X[i][0])
        utr5 = seq[: utr5_len + len(KOZAK)]
        parsed = mRNASequence(utr5=utr5, cds=cds, utr3=utr3)
        try:
            report = score_parsed(parsed)
            fitness = compute_fitness(report)
            results.append({"label": f"pareto_{i + 1}", "sequence": seq, "report": report, "fitness": fitness})
        except Exception:
            continue

    if not results:
        console.print("[bold red]Error:[/bold red] No sequences could be scored.")
        raise SystemExit(1)

    results.sort(key=lambda r: r["fitness"]["overall"], reverse=True)

    if csv_path:
        _write_csv(csv_path, history)
        console.print(f"History written to [bold]{csv_path}[/bold] ({len(history)} rows)\n")

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
