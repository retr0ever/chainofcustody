import rich_click as click
from rich.console import Console

from chainofcustody.initial import GeneNotFoundError, get_canonical_cds

console = Console()


@click.command()
@click.argument("gene_symbol")
def main(gene_symbol: str) -> None:
    """Fetch the canonical CDS for a given GENE_SYMBOL."""
    try:
        cds = get_canonical_cds(gene_symbol)
        console.print(cds)
    except GeneNotFoundError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
