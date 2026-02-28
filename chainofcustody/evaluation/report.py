"""Combine all metrics into a report with traffic-light summary."""

import json

from rich.console import Console
from rich.panel import Panel
from rich.rule import Rule
from rich.table import Table
from rich.text import Text

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
                 f"**Full (cap+polyA):** {info.get('full_length', info['total_length'])}bp | "
                 f"**5'UTR:** {info['utr5_length']}bp | "
                 f"**CDS:** {info['cds_length']}bp ({info['num_codons']} codons) | "
                 f"**3'UTR:** {info['utr3_length']}bp | "
                 f"**Cap:** {info.get('cap5_length', 0)}nt | "
                 f"**Poly-A:** {info.get('poly_a_length', 0)}nt")
    lines.append("")

    lines.append("## Summary")
    lines.append("")
    lines.append("| Metric | Status |")
    lines.append("|--------|--------|")
    for metric, status in summary.items():
        lines.append(f"| {metric.replace('_', ' ').title()} | {status} |")
    lines.append("")

    lines.append("## 1. Structure")
    utr5 = structure.get("utr5_accessibility", {})
    if utr5.get("mfe_per_nt") is not None:
        lines.append(f"- **5'UTR MFE/nt:** {utr5['mfe_per_nt']} kcal/mol/nt  (raw: {utr5['mfe']} kcal/mol, {utr5.get('utr5_length', '?')} nt) ({utr5['status']})")
    gmfe = structure.get("global_mfe", {})
    lines.append(f"- **Global MFE:** {gmfe.get('mfe', 'N/A')} kcal/mol "
                 f"({gmfe.get('mfe_per_nt', 'N/A')} per nt)")
    lines.append("")

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

    if stab:
        lines.append("## 3. Stability")
        lines.append(f"- **GC3 (wobble position):** {stab.get('gc3', 0):.1%}")
        lines.append(f"- **MFE per nt:** {stab.get('mfe_per_nt', 0):.4f} kcal/mol/nt")
        lines.append(f"- **Stability score:** {stab.get('stability_score', 0):.2f}")
    lines.append("")

    ribonn = report["ribonn_scores"]
    lines.append("## 4. Translation Efficiency")
    lines.append(f"- **Mean TE (all tissues):** {ribonn['mean_te']:.4f}")
    if ribonn.get("target_cell_type"):
        lines.append(f"- **Target cell type:** {ribonn['target_cell_type']}")
        lines.append(f"- **Target TE:** {ribonn.get('target_te', 0.0):.4f}")
        lines.append(f"- **Mean off-target TE:** {ribonn.get('mean_off_target_te', 0.0):.4f}")
    lines.append(f"- **Status:** {ribonn['status']}")

    return "\n".join(lines)


def report_to_json(report: dict, include_fitness: bool = True) -> str:
    """Serialise report to JSON string, optionally including fitness scores."""
    data = dict(report)
    if include_fitness:
        data["fitness"] = compute_fitness(report)
    return json.dumps(data, indent=2, default=str)


# ── Rich terminal output ─────────────────────────────────────────────────

STATUS_COLOURS = {"GREEN": "green", "AMBER": "yellow", "RED": "red", "GREY": "dim"}

_METRIC_LABELS = {
    "utr5_accessibility": "5'UTR accessibility",
    "manufacturability": "Manufacturability",
    "stability": "Stability",
    "specificity": "Translation eff.",
}


def _metric_value(metric: str, report: dict, fitness: dict) -> str:
    if metric == "utr5_accessibility":
        d = report["structure_scores"]["utr5_accessibility"]
        val = d.get("mfe_per_nt")
        return f"{val:.3f} /nt" if val is not None else "N/A"
    return f"{fitness['scores'][metric]['value']:.2f}"


def _metric_hint(metric: str, report: dict) -> str:
    if metric == "utr5_accessibility":
        return ">= -0.1/nt green, < -0.3/nt red"
    if metric == "manufacturability":
        v = report["manufacturing_scores"]["utr5_violations"]
        return f"{v} UTR violation{'s' if v != 1 else ''} (5'UTR only)"
    if metric == "stability":
        stab = report.get("stability_scores", {})
        return f"GC3 {stab.get('gc3', 0):.0%}, MFE/nt {stab.get('mfe_per_nt', 0):.3f}"
    if metric == "specificity":
        ribonn = report["ribonn_scores"]
        target = ribonn.get("target_cell_type", "target")
        target_te = ribonn.get("target_te", ribonn.get("mean_te", 0.0))
        mean_off = ribonn.get("mean_off_target_te", 0.0)
        return f"{target}: {target_te:.2f}, off-target: {mean_off:.2f}"
    return ""


def print_report(console: Console, report: dict, label: str | None = None) -> None:
    """Print a single sequence report using Rich formatting."""
    info = report["sequence_info"]
    summary = report["summary"]
    structure = report["structure_scores"]
    mfg = report["manufacturing_scores"]
    fitness = compute_fitness(report)

    console.print()

    title = f"Sequence Report{f' — {label}' if label else ''}"
    full_len = info.get("full_length", info["total_length"])
    cap_len = info.get("cap5_length", 0)
    polya_len = info.get("poly_a_length", 0)
    header = (
        f"[dim]cap {cap_len}nt[/dim] + "
        f"5'UTR {info['utr5_length']}  |  "
        f"CDS {info['cds_length']} ({info['num_codons']} codons)  |  "
        f"3'UTR {info['utr3_length']}  "
        f"+ [dim]polyA {polya_len}nt[/dim]  "
        f"= [bold]{full_len} nt[/bold]"
    )
    console.print(Panel.fit(header, title=title, border_style="blue"))

    table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
    table.add_column("Metric", style="bold")
    table.add_column("Value", justify="right")
    table.add_column("Status")
    table.add_column("", style="dim")

    for metric, status in summary.items():
        colour = STATUS_COLOURS.get(status, "dim")
        table.add_row(
            _METRIC_LABELS.get(metric, metric),
            _metric_value(metric, report, fitness),
            f"[bold {colour}]{status}[/]",
            _metric_hint(metric, report),
        )

    console.print()
    console.print(table)
    console.print()

    # Detail sections (only for non-GREEN metrics)
    if summary.get("utr5_accessibility") != "GREEN":
        console.print(Rule("Structure", style="dim"))
        utr5 = structure.get("utr5_accessibility", {})
        if utr5.get("mfe_per_nt") is not None:
            console.print(f"  5'UTR MFE/nt  [bold]{utr5['mfe_per_nt']:.4f} kcal/mol/nt[/]  (raw {utr5['mfe']} kcal/mol, {utr5.get('utr5_length','?')} nt)")
            console.print(f"  [dim]{utr5.get('message', '')}[/]")
        gmfe = structure.get("global_mfe", {})
        console.print(f"  Global MFE  {gmfe.get('mfe', 'N/A')} kcal/mol ({gmfe.get('mfe_per_nt', 'N/A')} per nt)")
        console.print()

    if summary.get("manufacturability") != "GREEN":
        console.print(Rule("Manufacturability", style="dim"))
        for label_, check in [
            ("GC windows", mfg["gc_windows"]),
            ("Homopolymers", mfg["homopolymers"]),
            ("Restriction sites", mfg["restriction_sites"]),
        ]:
            if check["pass"]:
                console.print(f"  [green]{label_}  pass[/]")
            else:
                n = len(check["violations"])
                if label_ == "Restriction sites":
                    sites = ", ".join(f"{v['enzyme']}@{v['position']}" for v in check["violations"])
                    suffix = f" ({sites})"
                elif label_ == "Homopolymers":
                    suffix = " runs >8nt"
                else:
                    suffix = " violations"
                console.print(f"  [yellow]{label_}  {n}{suffix}[/]")
        console.print()

    stab = report.get("stability_scores", {})
    if summary.get("stability") != "GREEN" and stab:
        console.print(Rule("Stability", style="dim"))
        console.print(f"  [bold]GC3 (wobble position)  {stab.get('gc3', 0):.1%}[/]")
        console.print(f"  MFE per nt  {stab.get('mfe_per_nt', 0):.4f} kcal/mol/nt")
        console.print(f"  [dim]Combined score  {stab.get('stability_score', 0):.2f}[/]")
        console.print()

    ribonn = report["ribonn_scores"]
    if summary.get("specificity") != "GREEN":
        console.print(Rule("Translation Efficiency", style="dim"))
        target = ribonn.get("target_cell_type", "target")
        target_te = ribonn.get("target_te", ribonn.get("mean_te", 0.0))
        mean_off = ribonn.get("mean_off_target_te", 0.0)
        console.print(f"  [bold]{target}  TE = {target_te:.4f}[/]")
        console.print(f"  mean off-target  {mean_off:.4f}")
        console.print(f"  [dim]{ribonn.get('message', '')}[/]")
        console.print()

    if fitness["suggestions"]:
        console.print(Rule("What to improve", style="dim"))
        for i, s in enumerate(fitness["suggestions"], 1):
            style = "bold red" if s["priority"] == "high" else "yellow" if s["priority"] == "medium" else "dim"
            console.print(f"  {i}. [{style}]{s['action']}[/]")
        console.print()

    score = fitness["overall"]
    score_style = "bold green" if score >= 0.7 else "bold yellow" if score >= 0.4 else "bold red"
    console.print(f"  [{score_style}]Overall fitness  {score:.2f}[/]")
    console.print()


def print_batch_report(console: Console, results: list[dict]) -> None:
    """Print a ranked comparison table for multiple candidates."""
    console.print()

    table = Table(title="Candidate Ranking", show_header=True, header_style="bold",
                  padding=(0, 1))
    table.add_column("Rank", justify="right", style="bold")
    table.add_column("Candidate")
    table.add_column("5'UTR", justify="center")
    table.add_column("Mfg", justify="right")
    table.add_column("Stab", justify="right")
    table.add_column("TE", justify="right")
    table.add_column("Overall", justify="right", style="bold")

    for i, r in enumerate(results, 1):
        report = r["report"]
        fitness = r["fitness"]
        summary = report["summary"]
        scores = fitness["scores"]
        overall = fitness["overall"]

        def _cell(metric: str) -> str:
            c = STATUS_COLOURS.get(summary.get(metric, "GREY"), "dim")
            return f"[{c}]{scores[metric]['value']:.2f}[/]"

        mfe = report["structure_scores"]["utr5_accessibility"].get("mfe")
        utr5_colour = STATUS_COLOURS.get(summary.get("utr5_accessibility", "GREY"), "dim")
        overall_colour = "green" if overall >= 0.7 else "yellow" if overall >= 0.4 else "red"

        table.add_row(
            str(i),
            r["label"],
            f"[{utr5_colour}]{mfe if mfe is not None else 'N/A'}[/]",
            _cell("manufacturability"),
            _cell("stability"),
            _cell("specificity"),
            f"[{overall_colour}]{overall:.2f}[/]",
        )

    console.print(table)
    _print_score_legend(console)


def _print_score_legend(console: Console) -> None:
    """Print a compact legend explaining each scored metric."""
    legend = Table(show_header=True, header_style="bold dim", padding=(0, 1),
                   box=None, show_edge=False)
    legend.add_column("Metric", style="dim")
    legend.add_column("Weight", justify="right", style="dim")
    legend.add_column("Target / Interpretation", style="dim")

    rows = [
        ("5'UTR", "utr5_accessibility",    "MFE/nt (5'UTR only); >= -0.1 accessible (green), < -0.3 over-structured (red)"),
        ("Mfg",   "manufacturability",     "5'UTR synthesis violations (GC windows, homopolymers, restriction sites); 0 ideal"),
        ("Stab",  "stability",             "mRNA stability 0→1 (GC3 wobble, MFE/nt)"),
        ("TE",    "specificity", "RiboNN selectivity: target TE − mean off-target TE; differential >= 1.5 ideal"),
        ("Score", None,                    "weighted sum of all metrics above"),
    ]

    for label, key, description in rows:
        legend.add_row(label, f"{DEFAULT_WEIGHTS[key]:.0%}" if key else "—", description)

    console.print(legend)
    console.print()
