"""Combine all metrics into a report with traffic-light summary."""

import json

from rich.console import Console
from rich.panel import Panel
from rich.rule import Rule
from rich.table import Table
from rich.text import Text

from .parser import ParsedSequence, parse_sequence
from .codons import score_codons
from .mirna import score_mirna
from .structure import score_structure
from .manufacturing import score_manufacturing
from .fitness import compute_fitness


def _traffic_light(value: float | None, green_range: tuple, amber_range: tuple) -> str:
    """Assign GREEN/AMBER/RED based on thresholds."""
    if value is None:
        return "GREY"
    g_lo, g_hi = green_range
    a_lo, a_hi = amber_range
    if g_lo <= value <= g_hi:
        return "GREEN"
    elif a_lo <= value <= a_hi:
        return "AMBER"
    return "RED"


def score_sequence(
    seq: str,
    offtarget: str = "HepG2",
    target: str | None = None,
    utr5_end: int | None = None,
    cds_end: int | None = None,
) -> dict:
    """
    Score an mRNA sequence across all 5 metrics.

    Args:
        seq: Raw sequence string or path to FASTA file.
        offtarget: Off-target cell type (default: HepG2 / liver).
        target: Target cell type column name (optional, placeholder).
        utr5_end: Manual 5'UTR boundary.
        cds_end: Manual CDS boundary.

    Returns: Full report dict.
    """
    parsed = parse_sequence(seq, utr5_end=utr5_end, cds_end=cds_end)

    # Metric 1 + 2: Codon quality and selectivity
    codon_scores = score_codons(parsed, target_cell_type=target)

    # Metric 3: miRNA detargeting
    mirna_scores = score_mirna(parsed)

    # Get miR-122 site positions for structure check
    mir122_positions = []
    if "miR-122-5p" in mirna_scores.get("detargeting", {}):
        mir122_positions = mirna_scores["detargeting"]["miR-122-5p"].get("positions", [])

    # Metric 4: Structure
    structure_scores = score_structure(parsed, mirna_site_positions=mir122_positions)

    # Metric 5: Manufacturability
    manufacturing_scores = score_manufacturing(parsed)

    # Traffic light summary
    cai = codon_scores.get("cai", 0)
    gc = codon_scores.get("gc_content", {}).get("cds", 0)
    mir122_count = mirna_scores.get("detargeting", {}).get("miR-122-5p", {}).get("utr3_sites", 0)
    utr5_mfe = structure_scores.get("utr5_accessibility", {}).get("mfe")
    mfg_violations = manufacturing_scores.get("total_violations", 0)

    summary = {
        "codon_quality": _traffic_light(cai, (0.8, 1.0), (0.6, 1.0)),
        "gc_content": _traffic_light(gc, (40, 60), (30, 70)),
        "mir122_detargeting": _traffic_light(mir122_count, (3, 100), (1, 100)),
        "utr5_accessibility": structure_scores.get("utr5_accessibility", {}).get("status", "GREY"),
        "manufacturability": _traffic_light(-mfg_violations, (-3, 0), (-999, 0)),
    }

    return {
        "sequence_info": {
            "total_length": len(parsed.raw),
            "utr5_length": len(parsed.utr5),
            "cds_length": len(parsed.cds),
            "utr3_length": len(parsed.utr3),
            "num_codons": len(parsed.codons),
        },
        "codon_scores": codon_scores,
        "mirna_scores": mirna_scores,
        "structure_scores": structure_scores,
        "manufacturing_scores": manufacturing_scores,
        "summary": summary,
    }


def format_report(report: dict) -> str:
    """Format a report dict as a readable markdown string."""
    lines = []
    info = report["sequence_info"]
    summary = report["summary"]
    codons = report["codon_scores"]
    mirna = report["mirna_scores"]
    structure = report["structure_scores"]
    mfg = report["manufacturing_scores"]

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

    # Codon details
    lines.append("## 1. Codon Quality")
    lines.append(f"- **CAI:** {codons['cai']}")
    gc = codons["gc_content"]
    lines.append(f"- **GC content:** overall {gc['overall']}% | CDS {gc['cds']}%")
    lines.append(f"- **Liver selectivity score:** {codons['liver_selectivity']}")
    if "target_selectivity" in codons:
        lines.append(f"- **Target selectivity score:** {codons['target_selectivity']}")
    if "selectivity_ratio" in codons:
        lines.append(f"- **Selectivity ratio (target/liver):** {codons['selectivity_ratio']}")
    clusters = codons.get("rare_codon_clusters", [])
    lines.append(f"- **Rare codon clusters:** {len(clusters)} found")
    lines.append("")

    # miRNA details
    lines.append("## 2. miRNA Detargeting")
    for mirna_name, details in mirna.get("detargeting", {}).items():
        lines.append(f"- **{mirna_name}:** {details['utr3_sites']} sites in 3'UTR "
                     f"({details['total_sites']} total)")
        if details.get("inter_site_spacing"):
            lines.append(f"  - Spacing: {details['inter_site_spacing']}nt")
            lines.append(f"  - Adequate spacing: {'Yes' if details['adequate_spacing'] else 'No'}")
    for warning in mirna.get("warnings", []):
        lines.append(f"- **WARNING:** {warning['message']}")
    lines.append("")

    # Structure details
    lines.append("## 3. Structure")
    utr5 = structure.get("utr5_accessibility", {})
    if utr5.get("mfe") is not None:
        lines.append(f"- **5'UTR MFE:** {utr5['mfe']} kcal/mol ({utr5['status']})")
    gmfe = structure.get("global_mfe", {})
    lines.append(f"- **Global MFE:** {gmfe.get('mfe', 'N/A')} kcal/mol "
                 f"({gmfe.get('mfe_per_nt', 'N/A')} per nt)")
    sites = structure.get("mirna_site_accessibility", [])
    if sites:
        accessible = sum(1 for s in sites if s["accessible"])
        lines.append(f"- **miRNA site accessibility:** {accessible}/{len(sites)} sites accessible")
    lines.append("")

    # Manufacturing details
    lines.append("## 4. Manufacturability")
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


def _cai_bar(value: float, width: int = 10) -> Text:
    filled = round(value * width)
    bar = Text()
    bar.append("▓" * filled, style="bold green" if value >= 0.8 else "bold yellow" if value >= 0.6 else "bold red")
    bar.append("░" * (width - filled), style="dim")
    return bar


def _metric_label(metric: str) -> str:
    labels = {
        "codon_quality": "Codon quality",
        "gc_content": "GC content",
        "mir122_detargeting": "miR-122 detargeting",
        "utr5_accessibility": "5'UTR accessibility",
        "manufacturability": "Manufacturability",
    }
    return labels.get(metric, metric)


def _metric_value(metric: str, report: dict) -> str:
    if metric == "codon_quality":
        return f"CAI {report['codon_scores']['cai']}"
    if metric == "gc_content":
        return f"{report['codon_scores']['gc_content']['cds']}%"
    if metric == "mir122_detargeting":
        sites = report["mirna_scores"]["detargeting"].get("miR-122-5p", {}).get("utr3_sites", 0)
        return f"{sites} site{'s' if sites != 1 else ''}"
    if metric == "utr5_accessibility":
        mfe = report["structure_scores"]["utr5_accessibility"].get("mfe")
        return f"{mfe} kcal" if mfe is not None else "N/A"
    if metric == "manufacturability":
        v = report["manufacturing_scores"]["total_violations"]
        return f"{v} issue{'s' if v != 1 else ''}"
    return ""


def _metric_hint(metric: str) -> str:
    hints = {
        "codon_quality": ">0.8 for green",
        "gc_content": "40-60% ideal",
        "mir122_detargeting": "need 3+ in 3'UTR",
        "utr5_accessibility": "> -20 for green",
        "manufacturability": "0 = green",
    }
    return hints.get(metric, "")


def print_report(console: Console, report: dict, label: str | None = None) -> None:
    """Print a single sequence report using Rich formatting."""
    info = report["sequence_info"]
    summary = report["summary"]
    codons = report["codon_scores"]
    mirna = report["mirna_scores"]
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
            _metric_value(metric, report),
            _styled_status(status),
            _metric_hint(metric),
        )

    console.print()
    console.print(table)
    console.print()

    # ── Detail sections (only for non-GREEN metrics)
    has_details = False

    # Codon quality detail
    if summary.get("codon_quality") != "GREEN" or summary.get("gc_content") != "GREEN":
        has_details = True
        console.print(Rule("Codon Quality", style="dim"))
        cai = codons["cai"]
        bar = _cai_bar(cai)
        line = Text("  CAI  ")
        line.append(f"{cai}  ", style="bold")
        line.append_text(bar)
        console.print(line)

        gc = codons["gc_content"]
        gc_parts = [f"overall {gc['overall']}%", f"CDS {gc['cds']}%"]
        if gc.get("utr5") is not None:
            gc_parts.append(f"5'UTR {gc['utr5']}%")
        if gc.get("utr3") is not None:
            gc_parts.append(f"3'UTR {gc['utr3']}%")
        console.print(f"  GC   {'  '.join(gc_parts)}")
        console.print(f"  Liver selectivity  {codons['liver_selectivity']:+.4f}")
        if "target_selectivity" in codons:
            console.print(f"  Target selectivity  {codons['target_selectivity']:+.4f}")
        if "selectivity_ratio" in codons:
            console.print(f"  Selectivity ratio  {codons['selectivity_ratio']:+.4f}")
        clusters = codons.get("rare_codon_clusters", [])
        if clusters:
            console.print(f"  Rare codon clusters  {len(clusters)} found", style="yellow")
        console.print()

    # miRNA detail
    mir122 = mirna.get("detargeting", {}).get("miR-122-5p", {})
    warnings = mirna.get("warnings", [])
    if summary.get("mir122_detargeting") != "GREEN" or warnings:
        has_details = True
        console.print(Rule("miRNA Detargeting", style="dim"))
        utr3 = mir122.get("utr3_sites", 0)
        total = mir122.get("total_sites", 0)
        console.print(f"  miR-122-5p  {utr3} site{'s' if utr3 != 1 else ''} in 3'UTR ({total} total)")
        if mir122.get("inter_site_spacing"):
            spacing = mir122["inter_site_spacing"]
            adequate = mir122.get("adequate_spacing", False)
            style = "" if adequate else "yellow"
            console.print(f"  Spacing  {spacing}nt  {'adequate' if adequate else 'too close'}", style=style)
        for w in warnings:
            console.print(f"  {w['mirna']}  {w['sites_found']} accidental site{'s' if w['sites_found'] != 1 else ''}", style="red")
        console.print()

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
        sites = structure.get("mirna_site_accessibility", [])
        if sites:
            accessible = sum(1 for s in sites if s["accessible"])
            console.print(f"  miRNA site accessibility  {accessible}/{len(sites)} accessible")
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
            enzymes = [v["enzyme"] for v in rs["violations"]]
            unique = list(dict.fromkeys(enzymes))
            positions = [f"{v['enzyme']}@{v['position']}" for v in rs["violations"]]
            console.print(f"  Restriction sites  {len(rs['violations'])} found ({', '.join(positions)})", style="yellow")
        else:
            console.print("  Restriction sites  pass", style="green")
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


def print_batch_report(console: Console, results: list[dict]) -> None:
    """Print a ranked comparison table for multiple candidates."""
    console.print()

    table = Table(title="Candidate Ranking", show_header=True, header_style="bold", padding=(0, 1))
    table.add_column("Rank", justify="right", style="bold")
    table.add_column("Candidate")
    table.add_column("CAI", justify="right")
    table.add_column("GC", justify="center")
    table.add_column("miR-122", justify="center")
    table.add_column("5'UTR", justify="center")
    table.add_column("Mfg", justify="center")
    table.add_column("Score", justify="right", style="bold")

    for i, r in enumerate(results, 1):
        report = r["report"]
        fitness = r["fitness"]
        summary = report["summary"]

        score = fitness["overall"]
        score_style = "green" if score >= 0.7 else "yellow" if score >= 0.4 else "red"

        table.add_row(
            str(i),
            r["label"],
            str(report["codon_scores"]["cai"]),
            _styled_status(summary["gc_content"]),
            _metric_value("mir122_detargeting", report),
            _styled_status(summary["utr5_accessibility"]),
            _styled_status(summary["manufacturability"]),
            Text(f"{score:.2f}", style=score_style),
        )

    console.print(table)
    console.print()
