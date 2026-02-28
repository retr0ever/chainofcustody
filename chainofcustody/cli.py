import csv
import json
from pathlib import Path

import rich_click as click
from rich.console import Console
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn

from chainofcustody.evaluation.fitness import compute_fitness
from chainofcustody.evaluation.report import print_batch_report, print_report
from chainofcustody.cds import GeneNotFoundError, get_canonical_cds
from chainofcustody.optimization import KOZAK, METRIC_NAMES, mRNASequence, SequenceProblem, run, run_rl, score_parsed
from chainofcustody.three_prime import generate_utr3
from chainofcustody.three_prime.cell_type_map import SEED_MAP_TO_RIBONN, seed_map_to_ribonn
from chainofcustody.progress import set_status_callback, set_best_score_callback

console = Console()

_CSV_COLUMNS = ["generation", "sequence", *METRIC_NAMES, "overall"]
_DEFAULT_TARGET = "Fibroblast"


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
    """Convert a DNA sequence to RNA (T -> U, uppercase)."""
    return seq.upper().replace("T", "U")


def _run_rl_pipeline(
    cds: str,
    utr3: str,
    ribonn_cell_type: str,
    utr5_min: int,
    utr5_max: int,
    rl_episodes: int,
    rl_batch_size: int,
    rl_lr: float,
    seed: int | None,
    output_fmt: str,
    csv_path: Path | None,
    ribonn_path: Path | None,
) -> None:
    """Run the RL (PPO) optimisation pipeline and display results."""
    n_batches = max(1, rl_episodes // rl_batch_size)
    console.print(
        f"Running RL (PPO) - "
        f"episodes=[bold]{rl_episodes}[/bold]  "
        f"batch_size=[bold]{rl_batch_size}[/bold]  "
        f"lr=[bold]{rl_lr}[/bold]  "
        f"batches=[bold]{n_batches}[/bold]\n"
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
        rl_task = progress.add_task(
            f"RL PPO  (batch 0/{n_batches})", total=n_batches,
            status="starting...", best_score="--",
        )

        def _on_status(msg: str) -> None:
            progress.update(rl_task, status=msg)

        def _on_best_score(score: float) -> None:
            progress.update(rl_task, best_score=f"{score:.3f}")

        set_status_callback(_on_status)
        set_best_score_callback(_on_best_score)
        try:
            best_sequences, best_fitnesses, history = run_rl(
                cds=cds,
                utr3=utr3,
                target_cell_type=ribonn_cell_type,
                utr5_min=utr5_min,
                utr5_max=utr5_max,
                n_episodes=rl_episodes,
                batch_size=rl_batch_size,
                lr=rl_lr,
                seed=seed,
                progress=progress,
                progress_task=rl_task,
            )
        finally:
            set_status_callback(None)
            set_best_score_callback(None)

    if not best_sequences:
        console.print("[bold red]Error:[/bold red] No sequences could be scored.")
        raise SystemExit(1)

    results = []
    for i, (utr5, fit_data) in enumerate(zip(best_sequences, best_fitnesses)):
        report = fit_data["report"]
        fitness = fit_data["fitness"]
        seq = utr5 + KOZAK + cds + utr3
        results.append({
            "label": f"rl_{i + 1}",
            "sequence": seq,
            "report": report,
            "fitness": fitness,
        })

    if csv_path:
        # Write RL training history (batch-level stats)
        _rl_csv_columns = ["batch", "episodes", "mean_reward", "max_reward", "best_reward"]
        with csv_path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=_rl_csv_columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(history)
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
        console.print(f"RL candidates: [bold]{len(results)}[/bold] scored sequences\n")
        if len(results) == 1:
            print_report(console, results[0]["report"], label=results[0]["label"])
        else:
            print_batch_report(console, results)


@click.command()
@click.option("--gene", default="POU5F1", show_default=True, help="HGNC gene symbol whose canonical CDS to optimise (e.g. TP53).")
@click.option(
    "--target",
    default=_DEFAULT_TARGET,
    show_default=True,
    help=(
        "Target cell type for both 3'UTR sponge design and translation efficiency scoring. "
        "Must appear in cell_type_seed_map.csv with a RiboNN mapping. "
        f"Available: {', '.join(sorted(SEED_MAP_TO_RIBONN))}."
    ),
)
@click.option("--method", type=click.Choice(["nsga3", "rl"]), default="nsga3", show_default=True, help="Optimisation method: nsga3 (evolutionary) or rl (PPO reinforcement learning).")
@click.option("--utr5-min", type=int, default=20, show_default=True, help="Minimum 5'UTR length the optimizer can produce.")
@click.option("--utr5-max", type=int, default=1000, show_default=True, help="Maximum 5'UTR length the optimizer can produce.")
@click.option("--utr5-init", type=int, default=200, show_default=True, help="Initial 5'UTR length seed for the population (Gaussian spread +-10%).")
@click.option("--pop-size", type=int, default=128, show_default=True, help="Population size (must be >= number of reference directions; default covers Das-Dennis n_partitions=4 -> 126 directions).")
@click.option("--n-gen", type=int, default=50, show_default=True, help="Number of generations.")
@click.option("--mutation-rate", type=float, default=0.05, show_default=True, help="Per-position mutation probability.")
@click.option("--max-length-delta", type=int, default=50, show_default=True, help="Maximum change in 5'UTR length per mutation event.")
@click.option("--seed", type=int, default=None, help="Random seed for reproducibility.")
@click.option("--workers", type=int, default=None, help="Parallel worker processes for fitness evaluation (default: all CPU cores). Pass 1 to disable parallelism.")
@click.option("--seed-from-data/--no-seed-from-data", default=True, show_default=True, help="Warm-start a portion of the initial population with top-TE 5'UTR sequences from the MOESM3 dataset.")
@click.option("--gradient-seed-steps", type=int, default=0, show_default=True, help="Run this many gradient-ascent steps through RiboNN to design warm-start 5'UTR seeds before NSGA-III (0 = disabled).")
@click.option("--rl-episodes", type=int, default=2000, show_default=True, help="[RL] Total number of 5'UTR episodes to generate during PPO training.")
@click.option("--rl-batch-size", type=int, default=64, show_default=True, help="[RL] Episodes per rollout batch (scored together for GPU efficiency).")
@click.option("--rl-lr", type=float, default=3e-4, show_default=True, help="[RL] Adam learning rate for PPO.")
@click.option("--output", "output_fmt", type=click.Choice(["summary", "json"]), default="summary", show_default=True, help="Output format.")
@click.option("--csv", "csv_path", type=click.Path(dir_okay=False, writable=True, path_type=Path), default=None, help="Write Pareto-front results to a CSV file.")
@click.option("--ribonn-output", "ribonn_path", type=click.Path(dir_okay=False, writable=True, path_type=Path), default=None, help="Write per-tissue RiboNN predictions for Pareto-front candidates to a CSV file.")
def main(gene: str, target: str, method: str, utr5_min: int, utr5_max: int, utr5_init: int, pop_size: int, n_gen: int, mutation_rate: float, max_length_delta: int, seed: int | None, workers: int | None, seed_from_data: bool, gradient_seed_steps: int, rl_episodes: int, rl_batch_size: int, rl_lr: float, output_fmt: str, csv_path: Path | None, ribonn_path: Path | None) -> None:
    """Run optimisation to evolve an optimal 5'UTR for a given gene."""
    if utr5_min > utr5_max:
        console.print(f"[bold red]Error:[/bold red] --utr5-min ({utr5_min}) must be <= --utr5-max ({utr5_max}).")
        raise SystemExit(1)

    if target not in SEED_MAP_TO_RIBONN:
        available = ", ".join(sorted(SEED_MAP_TO_RIBONN))
        console.print(
            f"[bold red]Error:[/bold red] Unknown --target {target!r}.\n"
            f"Valid options: {available}"
        )
        raise SystemExit(1)

    ribonn_cell_type = seed_map_to_ribonn(target)

    console.print(f"\nResolving canonical CDS for gene [bold]{gene}[/bold]...")
    try:
        cds = _to_rna(get_canonical_cds(gene))
    except GeneNotFoundError as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        raise SystemExit(1)

    console.print(f"Generating 3'UTR sponge for target cell type [bold]{target}[/bold]...")
    try:
        utr3 = generate_utr3(target)
    except (ValueError, RuntimeError) as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        raise SystemExit(1)

    console.print(
        f"CDS resolved: [bold]{len(cds)} nt[/bold]  "
        f"3'UTR: [bold]{len(utr3)} nt[/bold]  "
        f"5'UTR length: [bold]{utr5_min}-{utr5_max} nt[/bold] (evolved)\n"
        f"Target cell type: [bold]{target}[/bold] (RiboNN: [bold]{ribonn_cell_type}[/bold])\n"
    )

    if method == "rl":
        _run_rl_pipeline(
            cds=cds, utr3=utr3, ribonn_cell_type=ribonn_cell_type,
            utr5_min=utr5_min, utr5_max=utr5_max,
            rl_episodes=rl_episodes, rl_batch_size=rl_batch_size, rl_lr=rl_lr,
            seed=seed, output_fmt=output_fmt, csv_path=csv_path, ribonn_path=ribonn_path,
        )
        return

    console.print(
        f"Running NSGA3 - "
        f"pop_size=[bold]{pop_size}[/bold]  "
        f"n_gen=[bold]{n_gen}[/bold]  "
        f"mutation_rate=[bold]{mutation_rate}[/bold]  "
        f"max_length_delta=[bold]{max_length_delta}[/bold]  "
        f"workers=[bold]{workers if workers is not None else 'auto'}[/bold]  "
        f"seed_from_data=[bold]{seed_from_data}[/bold]  "
        f"gradient_seed_steps=[bold]{gradient_seed_steps}[/bold]\n"
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
            status="starting...", best_score="--",
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
                target_cell_type=ribonn_cell_type,
                initial_length=utr5_init,
                max_length_delta=max_length_delta,
                seed_from_data=seed_from_data,
                gradient_seed_steps=gradient_seed_steps,
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
            report = score_parsed(parsed, target_cell_type=ribonn_cell_type)
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
