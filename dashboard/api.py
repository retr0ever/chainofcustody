"""
FastAPI server for mRNA structure visualisation and optimisation.

Wraps the plotting functions from plot_secondary_structure.py to serve
1-D linear maps and 2-D secondary-structure plots as SVG over HTTP,
and exposes the GA optimisation pipeline.

Run:
    uv run uvicorn dashboard.api:app --port 8000
"""

from __future__ import annotations

import io
import re
import sys
import tempfile
import traceback
from pathlib import Path
from typing import Sequence, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plot_secondary_structure import (
    _clean_rna,
    _fold_mfe,
    _plot_structure_naview,
    PlotStyle,
    plot_mrna_construct,
    predict_and_plot_full_and_utrs,
)

# Import CDS lookup and Optimisation logic
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
try:
    from chainofcustody.cds import get_canonical_cds, GeneNotFoundError
    from chainofcustody.dashboard_api.api import optimize_and_plot
except ImportError:
    get_canonical_cds = None
    GeneNotFoundError = Exception
    optimize_and_plot = None

app = FastAPI(title="Chain of Custody – Structure & Optimisation API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Distinguishable colours for sponge sites, one per unique miRNA.
SITE_PALETTE = [
    "#6ee7b7", "#93c5fd", "#fca5a5", "#fde68a",
    "#c4b5fd", "#fdba74", "#67e88c", "#f9a8d4",
    "#a5f3fc", "#d9f99d", "#e9d5ff", "#fecdd3",
]
SPACER_COLOR = "#475569"


def _site_color_map(mirna_names: Sequence[str] | None) -> dict[str, str]:
    if not mirna_names:
        return {}
    unique = list(dict.fromkeys(mirna_names))
    return {n: SITE_PALETTE[i % len(SITE_PALETTE)] for i, n in enumerate(unique)}


def _detect_sites(seq: str) -> list[tuple[int, int, int]]:
    """Same regex as teammate's plot_mrna_construct: uppercase runs 15-30 nt."""
    return [(m.start(), m.end(), i) for i, m in enumerate(re.finditer(r"[A-Z]{15,30}", seq))]


class FoldRequest(BaseModel):
    utr5: str = ""
    cds: str = ""
    utr3: str
    mirna_names: list[str] | None = None


class FoldResponse(BaseModel):
    plot_1d: str | None = None
    plot_2d: str | None = None
    dot_bracket: str | None = None
    mfe: float | None = None


class GeneCdsResponse(BaseModel):
    ok: bool
    gene: str
    cds: str = ""
    error: str = ""


class OptimizeRequest(BaseModel):
    gene: str
    target_cell_type: str = "Dendritic_cell"
    utr5_min: int = 20
    utr5_max: int = 150
    n_gen: int = 10


def _fig_to_svg(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", bbox_inches="tight", transparent=True)
    plt.close(fig)
    buf.seek(0)
    return buf.read().decode("utf-8")


# ── 2-D: uses teammate's _plot_structure_naview directly ──────────────

def _plot_2d_custom(
    seq_5utr: str,
    seq_cds: str,
    seq_3utr: str,
    mirna_names: Sequence[str] | None = None,
) -> tuple[plt.Figure, str, float]:
    """
    Renders 2D secondary structure using predict_and_plot_full_and_utrs logic.
    """
    full_seq = _clean_rna(seq_5utr + seq_cds + seq_3utr)
    n5, ncds = len(seq_5utr), len(seq_cds)
    n_full = len(full_seq)
    
    # Fold
    structure, mfe = _fold_mfe(full_seq, temperature_c=37.0)
    
    # Build per-nucleotide colour array
    nt_colors = [SPACER_COLOR] * n_full
    for i in range(n5): nt_colors[i] = "#4A8DAD"
    for i in range(n5, n5 + ncds): nt_colors[i] = "#3D6880"
    for i in range(n5 + ncds, n_full): nt_colors[i] = "#D4635A"

    cmap = _site_color_map(mirna_names)
    sites = _detect_sites(seq_3utr)
    for start, end, idx in sites:
        name = mirna_names[idx % len(mirna_names)] if mirna_names else None
        color = cmap.get(name, "gold") if name else "gold"
        for pos in range(n5 + ncds + start, min(n5 + ncds + end, n_full)):
            nt_colors[pos] = color

    import RNA
    from plot_secondary_structure import _extract_xy, _dist, PlotStyle

    style = PlotStyle(
        max_pair_span=250,
        backbone_max_dist=15.0,
        node_size=6.0,
        backbone_width=0.7,
        pair_width=0.6,
        backbone_alpha=0.3,
        pair_alpha=0.6,
        backbone_color="#64748b",
        pair_color="#94a3b8",
        figsize=(12, 12),
    )

    pt = RNA.ptable(structure)
    coords = RNA.naview_xy_coordinates(structure)
    x, y = _extract_xy(coords, n_full)

    fig, ax = plt.subplots(figsize=style.figsize)
    fig.patch.set_alpha(0.0)
    ax.set_facecolor("none")
    ax.set_aspect("equal")
    ax.set_axis_off()

    # Backbone
    for i in range(1, n_full):
        d = _dist(x[i], y[i], x[i + 1], y[i + 1])
        if style.backbone_max_dist and d > style.backbone_max_dist:
            continue
        ax.plot(
            [x[i], x[i + 1]], [y[i], y[i + 1]],
            linewidth=style.backbone_width, alpha=style.backbone_alpha,
            color=style.backbone_color, zorder=1,
        )

    # Base pairs
    for i in range(1, n_full + 1):
        j = int(pt[i])
        if j > i:
            if style.max_pair_span and abs(j - i) > style.max_pair_span:
                continue
            ax.plot(
                [x[i], x[j]], [y[i], y[j]],
                linewidth=style.pair_width, alpha=style.pair_alpha,
                color=style.pair_color, zorder=2,
            )

    # Nucleotide dots
    ax.scatter(
        [x[i] for i in range(1, n_full + 1)],
        [y[i] for i in range(1, n_full + 1)],
        s=style.node_size,
        c=[nt_colors[i - 1] for i in range(1, n_full + 1)],
        zorder=3, edgecolors="none",
    )

    # Legend
    handles = [
        patches.Patch(facecolor="#4A8DAD", label="5' UTR"),
        patches.Patch(facecolor="#3D6880", label="CDS"),
        patches.Patch(facecolor="#D4635A", label="3' UTR"),
    ]
    if mirna_names:
        unique = list(dict.fromkeys(mirna_names))
        for n in unique:
            handles.append(patches.Patch(facecolor=cmap[n], label=n))
            
    leg = ax.legend(
        handles=handles, loc="lower right", fontsize=8,
        frameon=True, facecolor="#1e293b", edgecolor="#334155",
        labelcolor="white", ncol=2 if len(handles) > 6 else 1
    )
    leg.get_frame().set_alpha(0.85)

    ax.set_title(
        f"mRNA Secondary Structure  |  MFE = {mfe:.2f} kcal/mol  |  {n_full} nt",
        fontsize=14, color="white", pad=12,
    )
    plt.tight_layout()
    return fig, structure, mfe


@app.post("/api/fold")
def fold(req: FoldRequest) -> FoldResponse:
    # 1-D linear map (Full Gene)
    fig_1d = plot_mrna_construct(
        seq_5utr=req.utr5,
        seq_cds=req.cds,
        seq_3utr=req.utr3,
        mirna_names=req.mirna_names,
        show=False
    )
    # Fix 1D colors for dark mode dashboard
    for text in fig_1d.findobj(matplotlib.text.Text):
        text.set_color("white")
    for patch in fig_1d.findobj(patches.Rectangle):
        patch.set_edgecolor("white")
    fig_1d.patch.set_alpha(0.0)
    fig_1d.gca().set_facecolor("none")
    
    svg_1d = _fig_to_svg(fig_1d)

    # 2-D secondary structure
    fig_2d, structure, mfe = _plot_2d_custom(
        seq_5utr=req.utr5,
        seq_cds=req.cds,
        seq_3utr=req.utr3,
        mirna_names=req.mirna_names
    )
    svg_2d = _fig_to_svg(fig_2d)

    return FoldResponse(
        plot_1d=svg_1d,
        plot_2d=svg_2d,
        dot_bracket=structure,
        mfe=mfe,
    )


@app.get("/api/gene-cds/{gene_symbol}")
def get_cds(gene_symbol: str) -> GeneCdsResponse:
    if get_canonical_cds is None:
        return GeneCdsResponse(ok=False, gene=gene_symbol, error="CDS lookup module not available")
    try:
        cds = get_canonical_cds(gene_symbol)
        cds_rna = cds.replace("T", "U")
        return GeneCdsResponse(ok=True, gene=gene_symbol, cds=cds_rna)
    except GeneNotFoundError as e:
        return GeneCdsResponse(ok=False, gene=gene_symbol, error=str(e))
    except Exception as e:
        return GeneCdsResponse(ok=False, gene=gene_symbol, error=f"Internal error: {str(e)}")


@app.post("/api/optimize")
def optimize(req: OptimizeRequest):
    if optimize_and_plot is None:
        return {"ok": False, "error": "Optimisation module not available"}
    try:
        result = optimize_and_plot(
            gene=req.gene,
            target_cell_type=req.target_cell_type,
            utr5_min=req.utr5_min,
            utr5_max=req.utr5_max,
            n_gen=req.n_gen
        )
        return result
    except Exception as e:
        return {"ok": False, "error": str(e), "traceback": traceback.format_exc()}
