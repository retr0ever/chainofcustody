"""Normalised fitness scoring and suggestion engine for candidate ranking."""

DEFAULT_WEIGHTS = {
    "codon_quality": 0.20,
    "gc_content": 0.10,
    "mir122_detargeting": 0.30,
    "utr5_accessibility": 0.15,
    "manufacturability": 0.25,
}


def _normalise_cai(report: dict) -> float:
    """CAI is already 0-1."""
    return report["codon_scores"]["cai"]


def _normalise_gc(report: dict) -> float:
    """1.0 if 40-60%, linear falloff to 0 at 20% or 80%."""
    gc = report["codon_scores"]["gc_content"]["cds"]
    if 40 <= gc <= 60:
        return 1.0
    elif gc < 40:
        return max(0.0, (gc - 20) / 20)
    else:
        return max(0.0, (80 - gc) / 20)


def _normalise_mir122(report: dict) -> float:
    """min(utr3_sites / 3, 1.0) — 3+ sites in 3'UTR = perfect."""
    sites = report["mirna_scores"].get("detargeting", {}).get("miR-122-5p", {}).get("utr3_sites", 0)
    return min(sites / 3, 1.0)


def _normalise_utr5(report: dict) -> float:
    """1.0 if MFE > -20, linear to 0 at -50."""
    mfe = report["structure_scores"].get("utr5_accessibility", {}).get("mfe")
    if mfe is None:
        return 0.5  # no data — neutral
    if mfe > -20:
        return 1.0
    elif mfe < -50:
        return 0.0
    return (mfe + 50) / 30


def _normalise_manufacturing(report: dict) -> float:
    """max(0, 1 - violations / 10)."""
    violations = report["manufacturing_scores"]["total_violations"]
    return max(0.0, 1 - violations / 10)


NORMALISERS = {
    "codon_quality": _normalise_cai,
    "gc_content": _normalise_gc,
    "mir122_detargeting": _normalise_mir122,
    "utr5_accessibility": _normalise_utr5,
    "manufacturability": _normalise_manufacturing,
}


def compute_fitness(report: dict, weights: dict[str, float] | None = None) -> dict:
    """
    Compute normalised per-metric scores and weighted overall score.

    Returns:
        {
            "scores": {metric: {"value": float, "weight": float, "weighted": float, "status": str}},
            "overall": float,
            "suggestions": [{"metric": str, "priority": str, "action": str}],
        }
    """
    weights = weights or DEFAULT_WEIGHTS
    summary = report["summary"]

    scores = {}
    for metric, normaliser in NORMALISERS.items():
        value = round(normaliser(report), 3)
        w = weights.get(metric, 0)
        scores[metric] = {
            "value": value,
            "weight": w,
            "weighted": round(value * w, 4),
            "status": summary.get(metric, "GREY"),
        }

    overall = sum(s["weighted"] for s in scores.values())

    suggestions = _generate_suggestions(report, scores)

    return {
        "scores": scores,
        "overall": round(overall, 4),
        "suggestions": suggestions,
    }


def _priority_from_weight(weight: float) -> str:
    if weight >= 0.25:
        return "high"
    elif weight >= 0.15:
        return "medium"
    return "low"


def _generate_suggestions(report: dict, scores: dict) -> list[dict]:
    """Generate actionable suggestions for non-GREEN metrics."""
    suggestions = []

    for metric, info in scores.items():
        if info["status"] == "GREEN":
            continue

        priority = _priority_from_weight(info["weight"])
        action = _suggestion_for(metric, report)
        if action:
            suggestions.append({"metric": metric, "priority": priority, "action": action})

    # Sort by priority: high > medium > low
    order = {"high": 0, "medium": 1, "low": 2}
    suggestions.sort(key=lambda s: order.get(s["priority"], 3))
    return suggestions


def _suggestion_for(metric: str, report: dict) -> str | None:
    if metric == "codon_quality":
        cai = report["codon_scores"]["cai"]
        return f"Optimise synonymous codons for higher human CAI (current: {cai:.2f}, target: >0.8)"

    if metric == "gc_content":
        gc = report["codon_scores"]["gc_content"]["cds"]
        return f"Adjust GC content in CDS (current: {gc:.1f}%, target: 40-60%)"

    if metric == "mir122_detargeting":
        sites = report["mirna_scores"]["detargeting"].get("miR-122-5p", {}).get("utr3_sites", 0)
        need = max(0, 3 - sites)
        return f"Add {need}+ miR-122-5p seed sites (CACTCC) to 3'UTR with >=8nt spacing"

    if metric == "utr5_accessibility":
        mfe = report["structure_scores"]["utr5_accessibility"].get("mfe")
        if mfe is not None:
            return f"Reduce 5'UTR secondary structure (current MFE: {mfe:.1f}, target: > -20 kcal/mol)"
        return None

    if metric == "manufacturability":
        mfg = report["manufacturing_scores"]
        parts = []
        gc_v = len(mfg["gc_windows"]["violations"])
        hp_v = len(mfg["homopolymers"]["violations"])
        rs_v = len(mfg["restriction_sites"]["violations"])
        if gc_v:
            parts.append(f"{gc_v} high-GC windows")
        if hp_v:
            parts.append(f"{hp_v} homopolymer runs")
        if rs_v:
            enzymes = [v["enzyme"] for v in mfg["restriction_sites"]["violations"]]
            unique = list(dict.fromkeys(enzymes))
            parts.append(f"{rs_v} restriction sites ({', '.join(unique)})")
        return "Fix " + " + ".join(parts) if parts else None

    return None
