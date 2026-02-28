"""Combine all metrics into a report with traffic-light summary."""

import json

from rich.console import Console
from rich.panel import Panel
from rich.rule import Rule
from rich.table import Table
from rich.text import Text

from .scoring import score_parsed
from .fitness import DEFAULT_WEIGHTS, compute_fitness


def format_report(report: dict) -> str:
    """Format a report dict as a readable markdown string."""
    lines = []
    info = report["sequence_info"]
    summary = report["summary"]
    structure = report["structure_scores"]
    mfg = report["manufacturing_scores"]
    stab = report.get("stability_scores", {})

    lines.append("# mRNA Sequence Scoring Report")
    lines.append("")
    lines.append(f"**Total length:** {info['total_length']}bp | "
                 f"**5'UTR:** {info['utr5_length']}bp | "
                 f"**CDS:** {info['cds_length']}bp ({info['num_codons']} codons) | "
                 f"**3'UTR:** {info['utr3_length']}bp")
    lines.append("")

    # Traffic light summary
    lines.append("## Summary")
    lines.append("")
    lines.append("| Metric | Status |")
    lines.append("|--------|--------|")
    for metric, status in summary.items():
        label = metric.replace("_", " ").title()
        lines.append(f"| {label} | {status} |")
    lines.append("")

    # Structure details
    lines.append("## 1. Structure")
    utr5 = structure.get("utr5_accessibility", {})
    if utr5.get("mfe") is not None:
        lines.append(f"- **5'UTR MFE:** {utr5['mfe']} kcal/mol ({utr5['status']})")
    gmfe = structure.get("global_mfe", {})
    lines.append(f"- **Global MFE:** {gmfe.get('mfe', 'N/A')} kcal/mol "
                 f"({gmfe.get('mfe_per_nt', 'N/A')} per nt)")
    lines.append("")

    # Manufacturing details
    lines.append("## 2. Manufacturability")
    lines.append(f"- **Overall:** {'PASS' if mfg['overall_pass'] else 'FAIL'} "
                 f"({mfg['total_violations']} violations)")
    gc_v = mfg["gc_windows"]
    gc_label = "PASS" if gc_v["pass"] else f"FAIL ({len(gc_v['violations'])} violations)"
    lines.append(f"- **GC windows:** {gc_label}")
    hp = mfg["homopolymers"]
    hp_label = "PASS" if hp["pass"] else f"FAIL ({len(hp['violations'])} runs >8nt)"
    lines.append(f"- **Homopolymers:** {hp_label}")
    rs = mfg["restriction_sites"]
    rs_label = "PASS" if rs["pass"] else f"FAIL ({len(rs['violations'])} sites found)"
    lines.append(f"- **Restriction sites:** {rs_label}")
    lines.append("")

    # Stability details
    if stab:
        lines.append("## 3. Stability")
        lines.append(f"- **GC3 (wobble position):** {stab.get('gc3', 0):.1%}")
        lines.append(f"- **MFE per nt:** {stab.get('mfe_per_nt', 0):.4f} kcal/mol/nt")
        lines.append(f"- **AU-rich elements:** {stab.get('au_rich_elements', 0)} in 3'UTR")
        lines.append(f"- **Stability score:** {stab.get('stability_score', 0):.2f}")

    return "\n".join(lines)


def report_to_json(report: dict, include_fitness: bool = True) -> str:
    """Serialise report to JSON string, optionally including fitness scores."""
    data = dict(report)
    if include_fitness:
        data["fitness"] = compute_fitness(report)
    return json.dumps(data, indent=2, default=str)


# ── Rich terminal output ─────────────────────────────────────────────────

STATUS_COLOURS = {"GREEN": "green", "AMBER": "yellow", "RED": "red", "GREY": "dim"}


def _styled_status(status: str) -> Text:
    colour = STATUS_COLOURS.get(status, "dim")
    return Text(status, style=f"bold {colour}")


def _metric_label(metric: str) -> str:
    labels = {
        "utr5_accessibility": "5'UTR accessibility",
        "manufacturability": "Manufacturability",
        "stability": "Stability",
    }
    return labels.get(metric, metric)


def _metric_value(metric: str, report: dict, fitness: dict) -> str:
    if metric == "utr5_accessibility":
        mfe = report["structure_scores"]["utr5_accessibility"].get("mfe")
        return f"{mfe} kcal" if mfe is not None else "N/A"
    return f"{fitness['scores'][metric]['value']:.2f}"


def _metric_hint(metric: str, report: dict) -> str:
    if metric == "utr5_accessibility":
        return "< -30 for green"
    if metric == "manufacturability":
        v = report["manufacturing_scores"]["total_violations"]
        return f"{v} violation{'s' if v != 1 else ''}"
    if metric == "stability":
        stab = report.get("stability_scores", {})
        return f"GC3 {stab.get('gc3', 0):.0%}, {stab.get('au_rich_elements', 0)} AREs"
    return ""


def print_report(console: Console, report: dict, label: str | None = None) -> None:
    """Print a single sequence report using Rich formatting."""
    info = report["sequence_info"]
    summary = report["summary"]
    structure = report["structure_scores"]
    mfg = report["manufacturing_scores"]
    fitness = compute_fitness(report)

    console.print()

    # ── Header
    title = f"Sequence Report{f' — {label}' if label else ''}"
    header_text = (
        f"{info['total_length']} bp    "
        f"5'UTR {info['utr5_length']}  |  "
        f"CDS {info['cds_length']} ({info['num_codons']} codons)  |  "
        f"3'UTR {info['utr3_length']}"
    )
    console.print(Panel(header_text, title=title, border_style="blue", expand=False))

    # ── Summary table
    table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
    table.add_column("Metric", style="bold")
    table.add_column("Value", justify="right")
    table.add_column("Status")
    table.add_column("", style="dim")

    for metric, status in summary.items():
        table.add_row(
            _metric_label(metric),
            _metric_value(metric, report, fitness),
            _styled_status(status),
            _metric_hint(metric, report),
        )

    console.print()
    console.print(table)
    console.print()

    # ── Detail sections (only for non-GREEN metrics)
    has_details = False

    # Structure detail
    if summary.get("utr5_accessibility") != "GREEN":
        has_details = True
        console.print(Rule("Structure", style="dim"))
        utr5 = structure.get("utr5_accessibility", {})
        if utr5.get("mfe") is not None:
            console.print(f"  5'UTR MFE  {utr5['mfe']} kcal/mol", style="bold")
            console.print(f"  {utr5.get('message', '')}", style="dim")
        gmfe = structure.get("global_mfe", {})
        console.print(f"  Global MFE  {gmfe.get('mfe', 'N/A')} kcal/mol ({gmfe.get('mfe_per_nt', 'N/A')} per nt)")
        console.print()

    # Manufacturing detail
    if summary.get("manufacturability") != "GREEN":
        has_details = True
        console.print(Rule("Manufacturability", style="dim"))
        gc_v = mfg["gc_windows"]
        if not gc_v["pass"]:
            console.print(f"  GC windows  {len(gc_v['violations'])} violations", style="yellow")
        else:
            console.print("  GC windows  pass", style="green")
        hp = mfg["homopolymers"]
        if not hp["pass"]:
            console.print(f"  Homopolymers  {len(hp['violations'])} runs >8nt", style="yellow")
        else:
            console.print("  Homopolymers  pass", style="green")
        rs = mfg["restriction_sites"]
        if not rs["pass"]:
            positions = [f"{v['enzyme']}@{v['position']}" for v in rs["violations"]]
            console.print(f"  Restriction sites  {len(rs['violations'])} found ({', '.join(positions)})", style="yellow")
        else:
            console.print("  Restriction sites  pass", style="green")
        console.print()

    # Stability detail
    stab = report.get("stability_scores", {})
    if summary.get("stability") != "GREEN" and stab:
        has_details = True
        console.print(Rule("Stability", style="dim"))
        console.print(f"  GC3 (wobble position)  {stab.get('gc3', 0):.1%}", style="bold")
        console.print(f"  MFE per nt  {stab.get('mfe_per_nt', 0):.4f} kcal/mol/nt")
        are = stab.get("au_rich_elements", 0)
        if are > 0:
            console.print(f"  AU-rich elements  {are} in 3'UTR", style="yellow")
        else:
            console.print("  AU-rich elements  none", style="green")
        console.print(f"  Combined score  {stab.get('stability_score', 0):.2f}", style="dim")
        console.print()

    # ── Suggestions
    if fitness["suggestions"]:
        console.print(Rule("What to improve", style="dim"))
        for i, s in enumerate(fitness["suggestions"], 1):
            priority_style = "bold red" if s["priority"] == "high" else "yellow" if s["priority"] == "medium" else "dim"
            line = Text(f"  {i}. ")
            line.append(s["action"], style=priority_style)
            console.print(line)
        console.print()

    # ── Overall score
    score = fitness["overall"]
    score_style = "bold green" if score >= 0.7 else "bold yellow" if score >= 0.4 else "bold red"
    console.print(Text(f"  Overall fitness  {score:.2f}", style=score_style))
    console.print()


def _score_cell(value: float, status: str) -> Text:
    """Format a normalised 0-1 score with status colour."""
    colour = STATUS_COLOURS.get(status, "dim")
    return Text(f"{value:.2f}", style=colour)


def print_batch_report(console: Console, results: list[dict]) -> None:
    """Print a ranked comparison table for multiple candidates."""
    console.print()

    table = Table(title="Candidate Ranking", show_header=True, header_style="bold", padding=(0, 1))
    table.add_column("Rank", justify="right", style="bold")
    table.add_column("Candidate")
    table.add_column("5'UTR", justify="center")
    table.add_column("Mfg", justify="right")
    table.add_column("Stab", justify="right")
    table.add_column("Overall", justify="right", style="bold")

    for i, r in enumerate(results, 1):
        report = r["report"]
        fitness = r["fitness"]
        summary = report["summary"]
        scores = fitness["scores"]

        overall = fitness["overall"]
        overall_style = "green" if overall >= 0.7 else "yellow" if overall >= 0.4 else "red"

        mfe = report["structure_scores"]["utr5_accessibility"].get("mfe")
        utr5_status = summary.get("utr5_accessibility", "GREY")
        utr5_colour = STATUS_COLOURS.get(utr5_status, "dim")
        utr5_text = Text(f"{mfe}" if mfe is not None else "N/A", style=utr5_colour)

        table.add_row(
            str(i),
            r["label"],
            utr5_text,
            _score_cell(scores["manufacturability"]["value"], summary["manufacturability"]),
            _score_cell(scores["stability"]["value"], summary.get("stability", "GREY")),
            Text(f"{overall:.2f}", style=overall_style),
        )

    console.print(table)
    _print_score_legend(console)


def _print_score_legend(console: Console) -> None:
    """Print a compact legend explaining each scored metric."""
    legend = Table(show_header=True, header_style="bold dim", padding=(0, 1), box=None, show_edge=False)
    legend.add_column("Metric", style="dim")
    legend.add_column("Weight", justify="right", style="dim")
    legend.add_column("Target / Interpretation", style="dim")

    rows = [
        ("5'UTR", "utr5_accessibility", "MFE kcal/mol; > −20 means accessible cap for translation"),
        ("Mfg",   "manufacturability",  "synthesis violations (GC windows, homopolymers, restriction sites); 0 ideal"),
        ("Stab",  "stability",          "mRNA stability 0→1 (GC3 wobble, AU-rich elements, MFE/nt)"),
        ("Score", None,                 "weighted sum of all metrics above"),
    ]

    for label, key, description in rows:
        weight = f"{DEFAULT_WEIGHTS[key]:.0%}" if key else "—"
        legend.add_row(label, weight, description)

    console.print(legend)
    console.print()
