import mygene
import requests

_ENSEMBL_API = "https://rest.ensembl.org"
_HEADERS = {"Content-Type": "application/json"}
_mg = mygene.MyGeneInfo()


class GeneNotFoundError(Exception):
    pass


def _resolve_ensembl_gene_id(gene_symbol: str) -> str:
    result = _mg.query(gene_symbol, species="human", fields="ensembl.gene", size=1)
    hits = result.get("hits", [])
    if not hits:
        raise GeneNotFoundError(f"Gene '{gene_symbol}' not found")

    ensembl = hits[0].get("ensembl")
    if not ensembl:
        raise GeneNotFoundError(f"No Ensembl entry for gene '{gene_symbol}'")

    # ensembl field can be a single dict or a list when multiple Ensembl entries exist
    if isinstance(ensembl, list):
        return ensembl[0]["gene"]
    return ensembl["gene"]


def _lookup_canonical_transcript(ensembl_gene_id: str, gene_symbol: str) -> str:
    url = f"{_ENSEMBL_API}/lookup/id/{ensembl_gene_id}"
    response = requests.get(url, headers=_HEADERS, params={"expand": 1})
    response.raise_for_status()

    canonical_id = response.json().get("canonical_transcript")
    if not canonical_id:
        raise GeneNotFoundError(f"No canonical transcript found for gene '{gene_symbol}'")

    return canonical_id


def _fetch_cds(transcript_id: str) -> str:
    url = f"{_ENSEMBL_API}/sequence/id/{transcript_id}"
    response = requests.get(url, headers=_HEADERS, params={"type": "cds"})
    response.raise_for_status()
    return response.json()["seq"]


def get_canonical_cds(gene_symbol: str) -> str:
    """Return the canonical CDS sequence for a human gene symbol.

    Resolves the gene symbol via MyGene.info, fetches the canonical transcript
    from Ensembl, and returns its CDS (start codon through stop codon, no UTRs).

    Args:
        gene_symbol: HGNC gene symbol, e.g. "BRCA1".

    Returns:
        The CDS nucleotide sequence as a string.

    Raises:
        GeneNotFoundError: If the gene or its canonical transcript cannot be found.
        requests.HTTPError: On unexpected Ensembl API errors.
    """
    ensembl_gene_id = _resolve_ensembl_gene_id(gene_symbol)
    transcript_id = _lookup_canonical_transcript(ensembl_gene_id, gene_symbol)
    return _fetch_cds(transcript_id)
