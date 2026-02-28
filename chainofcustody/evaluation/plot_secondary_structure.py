from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple, Any, List

import math
import matplotlib.pyplot as plt
import RNA


# ----------------------------
# Utilities
# ----------------------------

def _clean_rna(seq: str) -> str:
    seq = (seq or "").strip().upper().replace("T", "U")
    bad = set(seq) - set("AUGCN")
    if bad:
        raise ValueError(f"Invalid characters {sorted(bad)}. Allowed: A,U,G,C,N (T->U).")
    return seq


def _normalize_coords(coords: Any, n: int) -> List[Any]:
    """
    ViennaRNA NAView coordinate vectors sometimes have length n+1 with a dummy entry.
    Normalize to a list of length n.
    """
    m = len(coords)
    if m == n:
        return list(coords)
    if m == n + 1:
        return list(coords[1:])  # drop dummy at index 0
    raise ValueError(f"Unexpected coordinate vector length: got {m}, expected {n} (or {n+1}).")


def _extract_xy(coords: Any, n: int) -> Tuple[List[float], List[float]]:
    """
    Extract x/y arrays (1-based, index 0 unused) from ViennaRNA coordinate objects.
    Handles common SWIG wrapper variants.
    """
    coords_n = _normalize_coords(coords, n)
    x = [0.0] * (n + 1)
    y = [0.0] * (n + 1)

    for i in range(1, n + 1):
        c = coords_n[i - 1]

        if hasattr(c, "X") and hasattr(c, "Y"):
            x[i] = float(c.X)
            y[i] = float(c.Y)
            continue

        if hasattr(c, "x") and hasattr(c, "y"):
            x[i] = float(c.x)
            y[i] = float(c.y)
            continue

        try:
            x[i] = float(c[0])
            y[i] = float(c[1])
            continue
        except Exception:
            pass

        attrs = [a for a in dir(c) if not a.startswith("_")]
        raise TypeError(
            f"Cannot extract numeric x/y from coord element type={type(c)}. "
            f"Public attrs sample: {attrs[:40]}"
        )

    return x, y


def _fold_mfe(seq: str, temperature_c: float) -> Tuple[str, float]:
    md = RNA.md()
    md.temperature = float(temperature_c)
    fc = RNA.fold_compound(seq, md)
    struct, mfe = fc.mfe()
    return struct, float(mfe)


def _dist(x1: float, y1: float, x2: float, y2: float) -> float:
    return math.hypot(x1 - x2, y1 - y2)


# ----------------------------
# Plotting
# ----------------------------

@dataclass(frozen=True)
class PlotStyle:
    # Pair filtering (readability). 0 = draw all.
    max_pair_span: int = 0
    # Backbone filtering (avoid teleports). 0 disables filtering.
    backbone_max_dist: float = 0.0

    # Appearance
    node_size: float = 8.0
    backbone_width: float = 0.9
    pair_width: float = 0.8
    backbone_alpha: float = 0.35
    pair_alpha: float = 0.9

    # IMPORTANT: explicit colors to avoid matplotlib rainbow cycling
    backbone_color: str = "0.75"   # light gray
    pair_color: str = "0.25"       # dark gray

    figsize: Tuple[float, float] = (10.0, 10.0)
    dpi: int = 300


def _plot_structure_naview(
    seq: str,
    structure: str,
    title: str,
    out_prefix: str,
    out_dir: Path,
    color_func,
    style: PlotStyle,
    draw_backbone: bool = True,
    draw_pairs: bool = True,
    save_png: bool = True,
    save_svg: bool = True,
    show: bool = False,
) -> Dict[str, str]:
    """
    Plot RNA structure using ViennaRNA NAView coordinates.

    Note: “diagonal” base pair lines are normal in NAView. If you want RNAplot-style ladders,
    you need RNAplot or a different layout engine.
    """
    n = len(seq)
    pt = RNA.ptable(structure)  # 1-based
    coords = RNA.naview_xy_coordinates(structure)
    x, y = _extract_xy(coords, n)

    fig, ax = plt.subplots(figsize=style.figsize)
    ax.set_aspect("equal")
    ax.set_axis_off()

    # Backbone: edge-by-edge; optionally filter long jumps
    if draw_backbone:
        for i in range(1, n):
            d = _dist(x[i], y[i], x[i + 1], y[i + 1])
            if style.backbone_max_dist and d > style.backbone_max_dist:
                continue
            ax.plot(
                [x[i], x[i + 1]],
                [y[i], y[i + 1]],
                linewidth=style.backbone_width,
                alpha=style.backbone_alpha,
                color=style.backbone_color,
                zorder=1,
            )

    # Base pairs: single consistent color
    if draw_pairs:
        for i in range(1, n + 1):
            j = int(pt[i])
            if j > i:
                if style.max_pair_span and abs(j - i) > style.max_pair_span:
                    continue
                ax.plot(
                    [x[i], x[j]],
                    [y[i], y[j]],
                    linewidth=style.pair_width,
                    alpha=style.pair_alpha,
                    color=style.pair_color,
                    zorder=2,
                )

    # Nucleotides on top
    ax.scatter(
        [x[i] for i in range(1, n + 1)],
        [y[i] for i in range(1, n + 1)],
        s=style.node_size,
        c=[color_func(i) for i in range(1, n + 1)],
        zorder=3,
    )

    ax.set_title(title)

    out_dir.mkdir(parents=True, exist_ok=True)
    files: Dict[str, str] = {}

    if save_png:
        png = out_dir / f"{out_prefix}.png"
        fig.savefig(png, dpi=style.dpi, bbox_inches="tight")
        files["png"] = str(png)
    if save_svg:
        svg = out_dir / f"{out_prefix}.svg"
        fig.savefig(svg, bbox_inches="tight")
        files["svg"] = str(svg)

    if show:
        plt.show()
    else:
        plt.close(fig)

    return files


# ----------------------------
# Public API
# ----------------------------

def predict_and_plot_full_and_utrs(
    utr5: str,
    cds: str,
    utr3: str,
    out_dir: str = "plots",
    base_prefix: str = "transcript",
    temperature_c: float = 37.0,
    segment_colors: Optional[Dict[str, str]] = None,
    utr5_color: str = "dodgerblue",
    utr3_color: str = "tomato",
    show: bool = False,
) -> Dict[str, object]:
    """
    Produces:
      - Full transcript fold+plot (colored by segment)
      - 5'UTR-only fold+plot
      - 3'UTR-only fold+plot
    """
    utr5 = _clean_rna(utr5)
    cds = _clean_rna(cds)
    utr3 = _clean_rna(utr3)
    if not cds:
        raise ValueError("CDS is empty.")

    if segment_colors is None:
        segment_colors = {"5UTR": "dodgerblue", "CDS": "black", "3UTR": "tomato"}

    out_path = Path(out_dir)

    # Full-length: readability hacks ON (suppress long-range, skip backbone teleports)
    style_full = PlotStyle(
        max_pair_span=250,
        backbone_max_dist=15.0,
        node_size=3.0,
        backbone_width=0.6,
        pair_width=0.5,
        backbone_alpha=0.25,
        pair_alpha=0.7,
        backbone_color="0.80",
        pair_color="0.20",
    )

    # UTR-only: no filtering (clean geometry). Backbone in light gray.
    style_utr = PlotStyle(
        max_pair_span=0,
        backbone_max_dist=0.0,
        node_size=8.0,
        backbone_width=0.9,
        pair_width=0.9,
        backbone_alpha=0.30,
        pair_alpha=0.90,
        backbone_color="0.82",
        pair_color="0.25",
    )

    # ---- Full transcript ----
    full_seq = utr5 + cds + utr3
    n5, ncds = len(utr5), len(cds)
    n_full = len(full_seq)
    full_structure, full_mfe = _fold_mfe(full_seq, temperature_c)

    seg_ranges = {
        "5UTR": (1, n5),
        "CDS": (n5 + 1, n5 + ncds),
        "3UTR": (n5 + ncds + 1, n_full),
    }

    def full_color(i: int) -> str:
        if seg_ranges["5UTR"][0] <= i <= seg_ranges["5UTR"][1]:
            return segment_colors["5UTR"]
        if seg_ranges["CDS"][0] <= i <= seg_ranges["CDS"][1]:
            return segment_colors["CDS"]
        return segment_colors["3UTR"]

    full_title = f"{base_prefix}_full | MFE = {full_mfe:.2f} kcal/mol | n={n_full}"
    full_files = _plot_structure_naview(
        seq=full_seq,
        structure=full_structure,
        title=full_title,
        out_prefix=f"{base_prefix}_full",
        out_dir=out_path,
        color_func=full_color,
        style=style_full,
        draw_backbone=True,
        draw_pairs=True,
        show=show,
    )

    # ---- 5'UTR only ----
    utr5_structure, utr5_mfe = _fold_mfe(utr5, temperature_c)
    utr5_title = f"{base_prefix}_5UTR | MFE = {utr5_mfe:.2f} kcal/mol | n={len(utr5)}"
    utr5_files = _plot_structure_naview(
        seq=utr5,
        structure=utr5_structure,
        title=utr5_title,
        out_prefix=f"{base_prefix}_5UTR",
        out_dir=out_path,
        color_func=lambda _: utr5_color,
        style=style_utr,
        # If you still hate the look: set draw_backbone=False for UTR plots
        draw_backbone=True,
        draw_pairs=True,
        show=show,
    )

    # ---- 3'UTR only ----
    utr3_structure, utr3_mfe = _fold_mfe(utr3, temperature_c)
    utr3_title = f"{base_prefix}_3UTR | MFE = {utr3_mfe:.2f} kcal/mol | n={len(utr3)}"
    utr3_files = _plot_structure_naview(
        seq=utr3,
        structure=utr3_structure,
        title=utr3_title,
        out_prefix=f"{base_prefix}_3UTR",
        out_dir=out_path,
        color_func=lambda _: utr3_color,
        style=style_utr,
        draw_backbone=True,
        draw_pairs=True,
        show=show,
    )

    return {
        "full": {"sequence": full_seq, "structure": full_structure, "mfe": full_mfe, "files": full_files},
        "utr5": {"sequence": utr5, "structure": utr5_structure, "mfe": utr5_mfe, "files": utr5_files},
        "utr3": {"sequence": utr3, "structure": utr3_structure, "mfe": utr3_mfe, "files": utr3_files},
    }


# ----------------------------
# Example usage
# ----------------------------
if __name__ == "__main__":
    res = predict_and_plot_full_and_utrs(
        utr5="GAGTAGTCCCTTCGCAAGCCCTCATTTCACCAGGCCCCCGGCTTGGGGCGCCTTCCTTCCCC",
        cds="ATGGCGGGACACCTGGCTTCGGATTTCGCCTTCTCGCCCCCTCCAGGTGGTGGAGGTGATGGGCCAGGGGGGCCGGAGCCGGGCTGGGTTGATCCTCGGACCTGGCTAAGCTTCCAAGGCCCTCCTGGAGGGCCAGGAATCGGGCCGGGGGTTGGGCCAGGCTCTGAGGTGTGGGGGATTCCCCCATGCCCCCCGCCGTATGAGTTCTGTGGGGGGATGGCGTACTGTGGGCCCCAGGTTGGAGTGGGGCTAGTGCCCCAAGGCGGCTTGGAGACCTCTCAGCCTGAGGGCGAAGCAGGAGTCGGGGTGGAGAGCAACTCCGATGGGGCCTCCCCGGAGCCCTGCACCGTCACCCCTGGTGCCGTGAAGCTGGAGAAGGAGAAGCTGGAGCAAAACCCGGAGGAGTCCCAGGACATCAAAGCTCTGCAGAAAGAACTCGAGCAATTTGCCAAGCTCCTGAAGCAGAAGAGGATCACCCTGGGATATACACAGGCCGATGTGGGGCTCACCCTGGGGGTTCTATTTGGGAAGGTATTCAGCCAAACGACCATCTGCCGCTTTGAGGCTCTGCAGCTTAGCTTCAAGAACATGTGTAAGCTGCGGCCCTTGCTGCAGAAGTGGGTGGAGGAAGCTGACAACAATGAAAATCTTCAGGAGATATGCAAAGCAGAAACCCTCGTGCAGGCCCGAAAGAGAAAGCGAACCAGTATCGAGAACCGAGTGAGAGGCAACCTGGAGAATTTGTTCCTGCAGTGCCCGAAACCCACACTGCAGCAGATCAGCCACATCGCCCAGCAGCTTGGGCTCGAGAAGGATGTGGTCCGAGTGTGGTTCTGTAACCGGCGCCAGAAGGGCAAGCGATCAAGCAGCGACTATGCACAACGAGAGGATTTTGAGGCTGCTGGGTCTCCTTTCTCAGGGGGACCAGTGTCCTTTCCTCTGGCCCCAGGGCCCCATTTTGGTACCCCAGGCTATGGGAGCCCTCACTTCACTGCACTGTACTCCTCGGTCCCTTTCCCTGAGGGGGAAGCCTTTCCCCCTGTCTCCGTCACCACTCTGGGCTCTCCCATGCATTCAAACTGA",
        utr3="GGTGCCTGCCCTTCTAGGAATGGGGGACAGGGGGAGGGGAGGAGCTAGGGAAAGAAAACCTGGAGTTTGTGCCAGGGTTTTTGGGATTAAGTTCTTCATTCACTAAGGAAGGAATTGGGAACACAAAGGGTGGGGGCAGGGGAGTTTGGGGCAACTGGTTGGAGGGAAGGTGAAGTTCAATGATGCTCTTGATTTTAATCCCACATCATGTATCACTTTTTTCTTAAATAAAGAAGCCTGGGACACAGTAGATAGACACACTTA",
        out_dir="plots",
        base_prefix="Oct4_canonical",
        show=False,
    )
    print("3'UTR:", res["utr3"]["mfe"], res["utr3"]["files"])