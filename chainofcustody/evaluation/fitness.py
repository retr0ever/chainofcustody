"""Normalised fitness scoring and suggestion engine for candidate ranking."""

DEFAULT_WEIGHTS = {
    "utr5_accessibility": 0.22,
    "manufacturability": 0.30,
    "stability": 0.35,
    "translation_efficiency": 0.13,
}


def _normalise_utr5(report: dict) -> float:
    """1.0 if MFE < -30, linear to 0 at 0."""
    mfe = report["structure_scores"].get("utr5_accessibility", {}).get("mfe")
    if mfe is None:
        return 0.5  # no data — neutral
    if mfe < -30:
        return 1.0
    elif mfe > 0:
        return 0.0
    return -mfe / 30


def _normalise_manufacturing(report: dict) -> float:
    """max(0, 1 - violations / 10)."""
    violations = report["manufacturing_scores"]["total_violations"]
    return max(0.0, 1 - violations / 10)


def _normalise_stability(report: dict) -> float:
    """Use the combined stability score directly (already 0-1)."""
    return report.get("stability_scores", {}).get("stability_score", 0.5)


def _normalise_te(report: dict) -> float:
    """1.0 if mean TE >= 2.0, linear to 0 at 0.5."""
    mean_te = report["ribonn_scores"]["mean_te"]
    if mean_te >= 2.0:
        return 1.0
    elif mean_te <= 0.5:
        return 0.0
    return (mean_te - 0.5) / 1.5


NORMALISERS = {
    "utr5_accessibility": _normalise_utr5,
    "manufacturability": _normalise_manufacturing,
    "stability": _normalise_stability,
    "translation_efficiency": _normalise_te,
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
    if weight >= 0.35:
        return "high"
    elif weight >= 0.25:
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
    order = {"high": 0, "medium": 1, "low": 2}
    suggestions.sort(key=lambda s: order.get(s["priority"], 3))
    return suggestions


def _suggestion_for(metric: str, report: dict) -> str | None:
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

    if metric == "stability":
        stab = report.get("stability_scores", {})
        parts = []
        if stab.get("gc3", 0) < 0.5:
            parts.append(f"increase GC3 wobble content (current: {stab['gc3']:.1%})")
        if stab.get("mfe_per_nt", 0) > -0.3:
            parts.append("increase thermodynamic stability")
        return "Improve stability: " + ", ".join(parts) if parts else "Improve overall mRNA stability"

    if metric == "translation_efficiency":
        mean_te = report["ribonn_scores"]["mean_te"]
        return f"Optimise 5'UTR and codon usage for higher translation efficiency (current mean TE: {mean_te:.2f}, target: >= 1.5)"

    return None
