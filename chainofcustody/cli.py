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
from chainofcustody.progress import set_status_callback, set_best_score_callback
from chainofcustody.evaluation.ribonn import get_valid_tissue_names

console = Console()

_CSV_COLUMNS = ["generation", "sequence", *METRIC_NAMES, "overall"]


def _write_csv(path: Path, history: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(history)


def _write_ribonn_csv(path: Path, results: list[dict]) -> None:
    """Write per-tissue RiboNN predictions for each Pareto-front candidate."""
    rows = []
    tissue_cols: list[str] = []
    for r in results:
        per_tissue: dict | None = r["report"]["ribonn_scores"].get("per_tissue")
        if per_tissue and not tissue_cols:
            tissue_cols = list(per_tissue.keys())
        row: dict = {
            "label": r["label"],
            "mean_te": r["report"]["ribonn_scores"]["mean_te"],
            "status": r["report"]["ribonn_scores"]["status"],
        }
        if per_tissue:
            row.update(per_tissue)
        rows.append(row)

    fieldnames = ["label", "mean_te", "status", *tissue_cols]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _to_rna(seq: str) -> str:
    """Convert a DNA sequence to RNA (T → U, uppercase)."""
    return seq.upper().replace("T", "U")


@click.command()
@click.option("--gene", default="POU5F1", show_default=True, help="HGNC gene symbol whose canonical CDS to optimise (e.g. TP53).")
@click.option("--off-target-cell-type", default="Hepatocyte", show_default=True, help="Off-target cell type used to generate the 3'UTR sponge (must match a name in the expression database).")
@click.option("--utr5-min", type=int, default=20, show_default=True, help="Minimum 5'UTR length the GA can produce.")
@click.option("--utr5-max", type=int, default=1000, show_default=True, help="Maximum 5'UTR length the GA can produce.")
@click.option("--utr5-init", type=int, default=200, show_default=True, help="Initial 5'UTR length seed for the population (Gaussian spread ±10%).")
@click.option("--pop-size", type=int, default=128, show_default=True, help="Population size (must be ≥ number of reference directions; default covers Das-Dennis n_partitions=4 → 126 directions).")
@click.option("--n-gen", type=int, default=50, show_default=True, help="Number of generations.")
@click.option("--mutation-rate", type=float, default=0.05, show_default=True, help="Per-position mutation probability.")
@click.option("--max-length-delta", type=int, default=50, show_default=True, help="Maximum change in 5'UTR length per mutation event.")
@click.option("--seed", type=int, default=None, help="Random seed for reproducibility.")
@click.option("--workers", type=int, default=None, help="Parallel worker processes for fitness evaluation (default: all CPU cores). Pass 1 to disable parallelism.")
@click.option("--output", "output_fmt", type=click.Choice(["summary", "json"]), default="summary", show_default=True, help="Output format.")
@click.option("--csv", "csv_path", type=click.Path(dir_okay=False, writable=True, path_type=Path), default=None, help="Write Pareto-front results to a CSV file.")
@click.option("--ribonn-output", "ribonn_path", type=click.Path(dir_okay=False, writable=True, path_type=Path), default=None, help="Write per-tissue RiboNN predictions for Pareto-front candidates to a CSV file.")
@click.option("--target-cell-type", default="megakaryocytes", show_default=True, help="Cell type in which high translation efficiency is desired (must exist in RiboNN tissue columns). All others are treated as off-target.")
def main(gene: str, off_target_cell_type: str, utr5_min: int, utr5_max: int, utr5_init: int, pop_size: int, n_gen: int, mutation_rate: float, max_length_delta: int, seed: int | None, workers: int | None, output_fmt: str, csv_path: Path | None, ribonn_path: Path | None, target_cell_type: str) -> None:
    """Run NSGA3 to evolve an optimal 5'UTR for a given gene."""
    if utr5_min > utr5_max:
        console.print(f"[bold red]Error:[/bold red] --utr5-min ({utr5_min}) must be ≤ --utr5-max ({utr5_max}).")
        raise SystemExit(1)

    valid_tissues = get_valid_tissue_names()
    if target_cell_type not in valid_tissues:
        console.print(
            f"[bold red]Error:[/bold red] Unknown --target-cell-type {target_cell_type!r}.\n"
            f"Valid options: {', '.join(sorted(valid_tissues))}"
        )
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
        f"Target cell type: [bold]{target_cell_type}[/bold]\n"
    )
    console.print(
        f"Running NSGA3 — "
        f"pop_size=[bold]{pop_size}[/bold]  "
        f"n_gen=[bold]{n_gen}[/bold]  "
        f"mutation_rate=[bold]{mutation_rate}[/bold]  "
        f"max_length_delta=[bold]{max_length_delta}[/bold]  "
        f"workers=[bold]{workers if workers is not None else 'auto'}[/bold]\n"
    )

    with Progress(
        SpinnerColumn(),
        TextColumn("[bold blue]{task.description}[/bold blue]"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        TextColumn("[bold green]best {task.fields[best_score]}[/bold green]"),
        TextColumn("[dim]{task.fields[status]}[/dim]"),
        console=console,
        transient=True,
    ) as progress:
        gen_task = progress.add_task(
            f"Evolving 5'UTR  (gen 0/{n_gen})", total=n_gen,
            status="starting…", best_score="—",
        )

        def _on_status(msg: str) -> None:
            progress.update(gen_task, status=msg)

        def _on_best_score(score: float) -> None:
            progress.update(gen_task, best_score=f"{score:.3f}")

        set_status_callback(_on_status)
        set_best_score_callback(_on_best_score)
        try:
            X, F, history = run(
                utr5_min=utr5_min, utr5_max=utr5_max, cds=cds, utr3=utr3,
                pop_size=pop_size, n_gen=n_gen,
                mutation_rate=mutation_rate, seed=seed, n_workers=workers,
                progress=progress, progress_task=gen_task,
                target_cell_type=target_cell_type,
                initial_length=utr5_init,
                max_length_delta=max_length_delta,
            )
        finally:
            set_status_callback(None)
            set_best_score_callback(None)

    problem = SequenceProblem(utr5_min=utr5_min, utr5_max=utr5_max, cds=cds, utr3=utr3)
    sequences = problem.decode(X)

    results = []
    for i, seq in enumerate(sequences):
        utr5_len = int(X[i][0])
        utr5 = seq[: utr5_len + len(KOZAK)]
        parsed = mRNASequence(utr5=utr5, cds=cds, utr3=utr3)
        try:
            report = score_parsed(parsed, target_cell_type=target_cell_type)
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

    if ribonn_path:
        _write_ribonn_csv(ribonn_path, results)
        console.print(f"RiboNN predictions written to [bold]{ribonn_path}[/bold] ({len(results)} candidates)\n")

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
