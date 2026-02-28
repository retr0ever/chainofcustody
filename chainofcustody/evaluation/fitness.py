"""Normalised fitness scoring and suggestion engine for candidate ranking."""

import math

DEFAULT_WEIGHTS = {
    "utr5_accessibility": 0.15,
    "manufacturability": 0.30,
    "stability": 0.20,
    "specificity": 0.35,
}


def _sigmoid(x: float, midpoint: float, k: float) -> float:
    """Logistic sigmoid: 1 / (1 + exp(-k * (x - midpoint))).

    Returns values in (0, 1) and never fully saturates, providing gradient
    even for extreme inputs far outside the expected range.
    k > 0  → higher x gives higher score
    k < 0  → higher x gives lower score (e.g. for violation counts)
    """
    return 1.0 / (1.0 + math.exp(-k * (x - midpoint)))


def _normalise_utr5(report: dict) -> float:
    """Sigmoid on MFE/nt of the 5'UTR (higher / less negative → more accessible → 1.0).

    Midpoint at -0.2 kcal/mol/nt; transitions from ~0.9 at -0.05 to ~0.1 at -0.35.
    """
    mfe_per_nt = report["structure_scores"].get("utr5_accessibility", {}).get("mfe_per_nt")
    if mfe_per_nt is None:
        return 0.5  # no data — neutral
    return _sigmoid(mfe_per_nt, midpoint=-0.2, k=15)


def _normalise_manufacturing(report: dict) -> float:
    """Sigmoid on 5'UTR-only violation count (fewer violations → 1.0).

    Midpoint at 1 violation; ~0.88 at 0 violations, ~0.12 at 2+ violations.
    """
    mfg = report["manufacturing_scores"]
    violations = mfg.get("utr5_violations", mfg["total_violations"])
    return _sigmoid(violations, midpoint=1.0, k=-2.0)


def _normalise_stability(report: dict) -> float:
    """Sigmoid on the combined stability score (higher → 1.0).

    Midpoint at 0.6; transitions from ~0.31 at 0.5 to ~0.92 at 0.9.
    """
    score = report.get("stability_scores", {}).get("stability_score", 0.5)
    return _sigmoid(score, midpoint=0.6, k=8.0)


def _normalise_te(report: dict) -> float:
    """Sigmoid on absolute target-tissue TE (higher → 1.0).

    RiboNN TE values observed across cell types typically range 0.1–2.5, with
    the CDS dominating the prediction and the 5'UTR contributing only ~1% of
    the variance between candidates in a single run.  The metric is therefore
    most useful for **final ranking** rather than as an optimisation gradient.

    Midpoint at 1.0; k=6 → ~3× more gradient per TE unit compared to the
    previous (midpoint=1.2, k=3) setting, giving NSGA-III a stronger
    directional signal in the practically achievable 5'UTR-tuning range
    (typically 0.8–1.6 TE units).  Score ~0.5 at TE=1.0, ~0.27 at TE=0.8,
    ~0.73 at TE=1.2, ~0.88 at TE=1.4.
    """
    ribonn = report["ribonn_scores"]
    target_te = ribonn.get("target_te", ribonn.get("mean_te", 0.0))
    return _sigmoid(target_te, midpoint=1.0, k=6.0)


NORMALISERS = {
    "utr5_accessibility": _normalise_utr5,
    "manufacturability": _normalise_manufacturing,
    "stability": _normalise_stability,
    "specificity": _normalise_te,
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
        d = report["structure_scores"]["utr5_accessibility"]
        val = d.get("mfe_per_nt")
        if val is not None:
            return f"Reduce 5'UTR secondary structure (current MFE/nt: {val:.3f}, target: >= -0.1 kcal/mol/nt)"
        return None

    if metric == "manufacturability":
        mfg = report["manufacturing_scores"]
        parts = []
        uorf_v = mfg.get("uorfs", {}).get("count", 0)
        gc_v = len(mfg["gc_windows"]["violations"])
        hp_v = len(mfg["homopolymers"]["violations"])
        rs_v = len(mfg["restriction_sites"]["violations"])
        if uorf_v:
            parts.append(f"{uorf_v} upstream AUG(s) (uORFs reduce main-ORF translation)")
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

    if metric == "specificity":
        ribonn = report["ribonn_scores"]
        target = ribonn.get("target_cell_type", "target")
        target_te = ribonn.get("target_te", ribonn.get("mean_te", 0.0))
        return (
            f"Improve {target} TE = {target_te:.2f} (target: >= 1.5)"
        )

    return None
