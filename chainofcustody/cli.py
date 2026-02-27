import sys
from pathlib import Path

import rich_click as click
from rich.console import Console

console = Console()


@click.group()
def main() -> None:
    """Chain of Custody — mRNA sequence design and evaluation."""


@main.command()
@click.argument("gene_symbol")
def fetch(gene_symbol: str) -> None:
    """Fetch the canonical CDS for a given GENE_SYMBOL."""
    from chainofcustody.initial import GeneNotFoundError, get_canonical_cds

    try:
        cds = get_canonical_cds(gene_symbol)
        console.print(cds)
    except GeneNotFoundError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise SystemExit(1)


@main.command()
@click.argument("inputs", nargs=-1)
@click.option("--gene", default=None, help="Gene symbol — fetches the CDS and scores it directly.")
@click.option("--offtarget", default="HepG2", help="Off-target cell type (default: HepG2 / liver).")
@click.option("--target", default=None, help="Target cell type column name from the dataset.")
@click.option("--utr5-end", type=int, default=None, help="Manual 5'UTR end position (0-indexed).")
@click.option("--cds-end", type=int, default=None, help="Manual CDS end position (0-indexed, exclusive).")
@click.option("--output", "output_fmt", type=click.Choice(["markdown", "json", "rich"]), default="rich", help="Output format (default: rich).")
def evaluate(inputs: tuple[str, ...], gene: str | None, offtarget: str, target: str | None, utr5_end: int | None, cds_end: int | None, output_fmt: str) -> None:
    """Score mRNA sequences across 5 metrics.

    \b
    Accepts input in several ways:
      chainofcustody evaluate sequence.fasta                     (one file)
      chainofcustody evaluate a.fasta b.fasta c.fasta           (multiple files)
      chainofcustody evaluate candidates/                        (directory of FASTA files)
      chainofcustody evaluate --gene POU5F1                     (fetch + score)
      chainofcustody fetch POU5F1 | chainofcustody evaluate     (pipe)
    """
    from chainofcustody.evaluation.report import (
        score_sequence, format_report, report_to_json, print_report, print_batch_report,
    )
    from chainofcustody.evaluation.fitness import compute_fitness

    sequences = _resolve_inputs(inputs, gene)

    if not sequences:
        console.print("[bold red]Error:[/bold red] No input. Provide a sequence, file path, directory, or --gene flag.")
        console.print()
        console.print("  chainofcustody evaluate sequence.fasta")
        console.print("  chainofcustody evaluate candidates/")
        console.print("  chainofcustody evaluate --gene POU5F1")
        raise SystemExit(1)

    # Score all sequences
    results = []
    for label, seq in sequences:
        try:
            report = score_sequence(seq=seq, offtarget=offtarget, target=target, utr5_end=utr5_end, cds_end=cds_end)
            fitness = compute_fitness(report)
            results.append({"label": label, "report": report, "fitness": fitness})
        except Exception as e:
            console.print(f"[bold red]Error scoring {label}:[/bold red] {e}")

    if not results:
        raise SystemExit(1)

    # Sort by fitness
    results.sort(key=lambda r: r["fitness"]["overall"], reverse=True)

    # Output
    if len(results) == 1:
        r = results[0]
        if output_fmt == "json":
            data = dict(r["report"])
            data["fitness"] = r["fitness"]
            console.print_json(report_to_json(r["report"]))
        elif output_fmt == "markdown":
            console.print(format_report(r["report"]))
        else:
            print_report(console, r["report"], label=r["label"])
    else:
        if output_fmt == "json":
            import json
            out = []
            for r in results:
                data = dict(r["report"])
                data["fitness"] = r["fitness"]
                data["label"] = r["label"]
                out.append(data)
            console.print_json(json.dumps(out, indent=2, default=str))
        elif output_fmt == "markdown":
            for r in results:
                console.print(f"\n## {r['label']}\n")
                console.print(format_report(r["report"]))
        else:
            print_batch_report(console, results)


@main.command()
@click.option("--seq-len", type=int, default=100, show_default=True, help="Length of each candidate sequence.")
@click.option("--pop-size", type=int, default=100, show_default=True, help="Population size.")
@click.option("--n-gen", type=int, default=50, show_default=True, help="Number of generations.")
@click.option("--mutation-rate", type=float, default=0.01, show_default=True, help="Per-position mutation probability.")
@click.option("--seed", type=int, default=None, help="Random seed for reproducibility.")
@click.option("--output", "output_fmt", type=click.Choice(["summary", "json"]), default="summary", show_default=True, help="Output format.")
def optimize(seq_len: int, pop_size: int, n_gen: int, mutation_rate: float, seed: int | None, output_fmt: str) -> None:
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
        f"mutation_rate=[bold]{mutation_rate}[/bold]\n"
    )

    X, F = run(seq_len=seq_len, pop_size=pop_size, n_gen=n_gen, mutation_rate=mutation_rate, seed=seed)

    problem = SequenceProblem(seq_len=seq_len)
    sequences = problem.decode(X)

    # Score each Pareto-front sequence with the full evaluation pipeline
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


def _resolve_inputs(inputs: tuple[str, ...], gene: str | None) -> list[tuple[str, str]]:
    """Resolve CLI inputs into a list of (label, sequence) pairs."""
    sequences = []

    if gene:
        if inputs:
            console.print("[bold red]Error:[/bold red] Provide either files or --gene, not both.")
            raise SystemExit(1)
        from chainofcustody.initial import GeneNotFoundError, get_canonical_cds
        console.print(f"Fetching CDS for [bold]{gene}[/bold]...")
        try:
            seq = get_canonical_cds(gene)
            return [(gene, seq)]
        except GeneNotFoundError as e:
            console.print(f"[bold red]Error:[/bold red] {e}")
            raise SystemExit(1)
        except Exception as e:
            console.print(f"[bold red]Error fetching gene:[/bold red] {e}")
            raise SystemExit(1)

    if inputs:
        for inp in inputs:
            p = Path(inp)
            if p.is_dir():
                # Collect all FASTA-like files in directory
                for f in sorted(p.iterdir()):
                    if f.suffix.lower() in (".fasta", ".fa", ".fna", ".txt", ".seq"):
                        sequences.append((f.stem, str(f)))
            elif p.is_file():
                sequences.append((p.stem, str(p)))
            else:
                # Treat as raw sequence
                sequences.append(("input", inp))
        return sequences

    # No arguments — try stdin
    if not sys.stdin.isatty():
        data = sys.stdin.read().strip()
        if data:
            return [("stdin", data)]

    return []
