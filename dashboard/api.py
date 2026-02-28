"""
FastAPI server for mRNA structure visualisation.

Wraps the plotting functions from plot_secondary_structure.py to serve
1-D linear maps and 2-D secondary-structure plots as SVG over HTTP.

Run:
    uv run uvicorn dashboard.api:app --port 8000
"""

from __future__ import annotations

import io
import re
import sys
import tempfile
from pathlib import Path
from typing import Sequence

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
)

app = FastAPI(title="Chain of Custody – Structure API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Distinguishable colours for sponge sites, one per unique miRNA.
SITE_PALETTE = [
    "#6ee7b7", "#93c5fd", "#fca5a5", "#fde68a",
    "#c4b5fd", "#fdba74", "#67e8f9", "#f9a8d4",
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
    sequence: str
    mirna_names: list[str] | None = None


class FoldResponse(BaseModel):
    plot_1d: str | None = None
    plot_2d: str | None = None
    dot_bracket: str | None = None
    mfe: float | None = None


def _fig_to_svg(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", bbox_inches="tight", transparent=True)
    plt.close(fig)
    buf.seek(0)
    return buf.read().decode("utf-8")


# ── 1-D: based on teammate's plot_mrna_construct ──────────────────────
# Follows the same structure but adapted for cassette-only view
# (no 5'UTR / CDS) and with per-miRNA site colours + legend.

def _plot_1d(
    seq_3utr: str,
    mirna_names: Sequence[str] | None = None,
) -> plt.Figure:
    total_len = len(seq_3utr)
    cmap = _site_color_map(mirna_names)
    sites = _detect_sites(seq_3utr)

    fig, ax = plt.subplots(figsize=(18, 5))
    fig.patch.set_alpha(0.0)
    ax.set_facecolor("none")
    ax.set_xlim(-total_len * 0.05, total_len * 1.05)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Backbone (same as teammate)
    ax.plot([0, total_len], [0.4, 0.4], color="white", linewidth=2)

    # 3'UTR domain block (teammate draws all three; we only have this one)
    rect = patches.Rectangle(
        (0, 0.3), total_len, 0.2,
        facecolor="lightcoral", edgecolor="white", alpha=0.35, zorder=2,
    )
    ax.add_patch(rect)
    ax.text(
        total_len / 2, 0.55, "3\u2032 UTR",
        ha="center", va="bottom", fontsize=12, fontweight="bold", color="white",
    )

    # Detect and label sponge sites — same regex as teammate's code
    site_counter = 0
    for match in re.finditer(r"[A-Z]{15,30}", seq_3utr):
        site_start = match.start()
        site_width = len(match.group())

        # Per-miRNA colour instead of uniform gold
        mirna_name = mirna_names[site_counter % len(mirna_names)] if mirna_names else None
        color = cmap.get(mirna_name, "gold") if mirna_name else "gold"

        site_rect = patches.Rectangle(
            (site_start, 0.3), site_width, 0.2,
            facecolor=color, edgecolor="white", linewidth=0.6, zorder=3,
        )
        ax.add_patch(site_rect)

        # Label — same staggered approach as teammate
        if mirna_names:
            mirna_label = mirna_names[site_counter % len(mirna_names)]
            label = f"{mirna_label} (#{site_counter + 1})"
        else:
            label = f"Site {site_counter + 1}"

        y_text = 0.65 if site_counter % 2 == 0 else 0.85
        ax.text(
            site_start + site_width / 2, y_text, label,
            ha="center", va="bottom", fontsize=8, rotation=45, color="white",
        )
        ax.plot(
            [site_start + site_width / 2, site_start + site_width / 2],
            [0.5, y_text - 0.02],
            color="gray", lw=0.5, zorder=1,
        )
        site_counter += 1

    # Legend — coloured patches for each miRNA
    if mirna_names:
        unique = list(dict.fromkeys(mirna_names))
        handles = [
            patches.Patch(facecolor=cmap[n], edgecolor="white", linewidth=0.5, label=n)
            for n in unique
        ]
        handles.append(patches.Patch(facecolor="lightcoral", alpha=0.35, edgecolor="white", linewidth=0.5, label="Spacer"))
        leg = ax.legend(
            handles=handles, loc="lower center",
            ncol=min(len(handles), 5), fontsize=8,
            frameon=True, facecolor="#1e293b", edgecolor="#334155",
            labelcolor="white", bbox_to_anchor=(0.5, -0.04),
        )
        leg.get_frame().set_alpha(0.85)

    ax.set_title(
        f"mRNA Sponge 3\u2032UTR  |  {total_len} nt  |  {site_counter} binding sites",
        fontsize=14, pad=20, color="white",
    )
    plt.tight_layout()
    return fig


# ── 2-D: uses teammate's _plot_structure_naview directly ──────────────
# Only change: colour function maps nucleotides by sponge site identity
# instead of by base (A/U/G/C), plus transparent background and legend.

def _plot_2d(
    seq: str,
    structure: str,
    mfe: float,
    original_seq: str,
    mirna_names: Sequence[str] | None = None,
) -> plt.Figure:
    n = len(seq)
    cmap = _site_color_map(mirna_names)
    sites = _detect_sites(original_seq)

    # Build per-nucleotide colour array: site colour or spacer grey
    nt_colors = [SPACER_COLOR] * n
    for start, end, idx in sites:
        name = mirna_names[idx % len(mirna_names)] if mirna_names else None
        color = cmap.get(name, "gold") if name else "gold"
        for pos in range(start, min(end, n)):
            nt_colors[pos] = color

    # Use teammate's _plot_structure_naview with custom colour function
    style = PlotStyle(
        node_size=8.0,
        backbone_width=0.9,
        pair_width=0.8,
        backbone_alpha=0.3,
        pair_alpha=0.5,
        backbone_color="#64748b",
        pair_color="#94a3b8",
        figsize=(10.0, 10.0),
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        color_func = lambda i: nt_colors[i - 1]

        files = _plot_structure_naview(
            seq=seq,
            structure=structure,
            title=f"mRNA Sponge 3\u2032UTR Fold  |  MFE = {mfe:.2f} kcal/mol  |  {n} nt",
            out_prefix="structure_2d",
            out_dir=Path(tmpdir),
            color_func=color_func,
            style=style,
            save_png=False,
            save_svg=True,
            show=False,
        )

        # Read the SVG, then patch in transparent background and legend
        with open(files["svg"]) as f:
            svg_raw = f.read()

    # The teammate's function saves to file with white bg. We need to
    # re-render with transparent bg and a legend, so we re-plot using
    # the same approach but return the figure directly.
    import RNA
    from plot_secondary_structure import _extract_xy, _dist

    pt = RNA.ptable(structure)
    coords = RNA.naview_xy_coordinates(structure)
    x, y = _extract_xy(coords, n)

    fig, ax = plt.subplots(figsize=style.figsize)
    fig.patch.set_alpha(0.0)
    ax.set_facecolor("none")
    ax.set_aspect("equal")
    ax.set_axis_off()

    # Backbone
    for i in range(1, n):
        d = _dist(x[i], y[i], x[i + 1], y[i + 1])
        if d > 15.0:
            continue
        ax.plot(
            [x[i], x[i + 1]], [y[i], y[i + 1]],
            linewidth=style.backbone_width, alpha=style.backbone_alpha,
            color=style.backbone_color, zorder=1,
        )

    # Base pairs
    for i in range(1, n + 1):
        j = int(pt[i])
        if j > i:
            ax.plot(
                [x[i], x[j]], [y[i], y[j]],
                linewidth=style.pair_width, alpha=style.pair_alpha,
                color=style.pair_color, zorder=2,
            )

    # Nucleotide dots
    ax.scatter(
        [x[i] for i in range(1, n + 1)],
        [y[i] for i in range(1, n + 1)],
        s=style.node_size,
        c=[nt_colors[i - 1] for i in range(1, n + 1)],
        zorder=3, edgecolors="none",
    )

    # Legend
    if mirna_names:
        unique = list(dict.fromkeys(mirna_names))
        handles = [
            patches.Patch(facecolor=cmap[n], edgecolor="none", label=n)
            for n in unique
        ]
        handles.append(patches.Patch(facecolor=SPACER_COLOR, edgecolor="none", label="Spacer"))
        leg = ax.legend(
            handles=handles, loc="lower right", fontsize=8,
            frameon=True, facecolor="#1e293b", edgecolor="#334155",
            labelcolor="white",
        )
        leg.get_frame().set_alpha(0.85)

    ax.set_title(
        f"mRNA Sponge 3\u2032UTR Fold  |  MFE = {mfe:.2f} kcal/mol  |  {n} nt",
        fontsize=13, color="white", pad=12,
    )
    plt.tight_layout()
    return fig


@app.post("/api/fold")
def fold(req: FoldRequest) -> FoldResponse:
    seq = _clean_rna(req.sequence)

    # 1-D linear map
    fig_1d = _plot_1d(seq_3utr=req.sequence, mirna_names=req.mirna_names)
    svg_1d = _fig_to_svg(fig_1d)

    # 2-D secondary structure
    structure, mfe = _fold_mfe(seq, temperature_c=37.0)
    fig_2d = _plot_2d(
        seq=seq, structure=structure, mfe=mfe,
        original_seq=req.sequence, mirna_names=req.mirna_names,
    )
    svg_2d = _fig_to_svg(fig_2d)

    return FoldResponse(
        plot_1d=svg_1d,
        plot_2d=svg_2d,
        dot_bracket=structure,
        mfe=mfe,
    )
