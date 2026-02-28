from __future__ import annotations

import time
import traceback
from pathlib import Path
from typing import Any, Dict, List, Optional


def _now_run_name() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _to_rna(seq: str) -> str:
    return (seq or "").strip().upper().replace("T", "U")


def _base_payload(run_dir: str) -> Dict[str, Any]:
    return {
        "ok": False,
        "error": "",
        "traceback": "",
        "run_dir": run_dir,
        "lengths": {"cds": 0, "utr3": 0},
        "history": [],
        "best": {"label": "none", "overall": None},
        "mirna_sites": [],
        "plots": {
            "full": {"files": {"png": ""}},
            "utr5": {"files": {"png": ""}},
            "utr3": {"files": {"png": ""}},
        },
    }


def _try_get_mirna_sites(target_cell_type: str) -> List[Dict[str, Any]]:
    """
    Best-effort adapter for three_prime/filtering_on_target.py.
    Replace this with your real function call when you're awake.
    """
    try:
        mod = __import__("chainofcustody.three_prime.filtering_on_target", fromlist=["*"])
    except Exception:
        return []

    for fn_name in (
        "get_mirna_binding_sites",
        "mirna_binding_sites",
        "list_mirna_sites",
        "compute_mirna_sites",
    ):
        fn = getattr(mod, fn_name, None)
        if callable(fn):
            try:
                sites = fn(target_cell_type)
                if sites is None:
                    return []
                if isinstance(sites, list):
                    return [s if isinstance(s, dict) else {"site": str(s)} for s in sites]
                return [{"site": str(sites)}]
            except Exception:
                return []

    for attr in ("SITES", "MIRNA_SITES", "mirna_sites"):
        val = getattr(mod, attr, None)
        if isinstance(val, list):
            return [s if isinstance(s, dict) else {"site": str(s)} for s in val]

    return []


def optimize_and_plot(
    gene: str = "POU5F1",
    target_cell_type: str = "Dendritic_cell",
    utr5_min: int = 10,
    utr5_max: int = 100,
    pop_size: int = 128,
    n_gen: int = 10,
    out_dir: str = "runs",
    run_name: Optional[str] = None,
    plot_temperature_c: float = 37.0,
) -> Dict[str, Any]:
    """
    Dashboard-safe wrapper:
    - Always returns a predictable dict.
    - Generates secondary structure PNGs via evaluation/plot_secondary_structure.py
    """
    if run_name is None or not str(run_name).strip():
        run_name = _now_run_name()

    run_dir = Path(out_dir) / run_name
    plot_dir = run_dir / "plots"
    _ensure_dir(plot_dir)

    payload = _base_payload(str(run_dir))

    if utr5_min > utr5_max:
        payload["error"] = f"utr5_min ({utr5_min}) must be <= utr5_max ({utr5_max})."
        return payload

    # Hardcode mutation rate for now since you removed it from the UI
    MUTATION_RATE = 0.01

    try:
        from chainofcustody.cds import GeneNotFoundError, get_canonical_cds
        from chainofcustody.three_prime import generate_utr3
        from chainofcustody.optimization import KOZAK, METRIC_NAMES, mRNASequence, SequenceProblem, run, score_parsed
        from chainofcustody.evaluation.fitness import compute_fitness

        # Add project root to path to reach dashboard/ folder
        import sys
        root = Path(__file__).resolve().parent.parent.parent
        if str(root) not in sys.path:
            sys.path.insert(0, str(root))
        from dashboard.plot_secondary_structure import predict_and_plot_full_and_utrs

        try:
            cds = _to_rna(get_canonical_cds(gene))
        except GeneNotFoundError as exc:
            payload["error"] = f"Gene not found: {exc}"
            return payload

        utr3 = _to_rna(generate_utr3(target_cell_type))

        # miRNA binding sites list (best-effort)
        payload["mirna_sites"] = _try_get_mirna_sites(target_cell_type)

        # Run GA (no progress; no workers; no seed)
        X, F, history = run(
            utr5_min=utr5_min,
            utr5_max=utr5_max,
            cds=cds,
            utr3=utr3,
            pop_size=pop_size,
            n_gen=n_gen,
            mutation_rate=MUTATION_RATE,
            seed=None,
            n_workers=None,
            progress=None,
            progress_task=None,
        )

        problem = SequenceProblem(utr5_min=utr5_min, utr5_max=utr5_max, cds=cds, utr3=utr3)
        sequences = problem.decode(X)

        scored: List[Dict[str, Any]] = []
        for i, seq in enumerate(sequences):
            utr5_len = int(X[i][0])
            utr5 = seq[: utr5_len + len(KOZAK)]
            parsed = mRNASequence(utr5=utr5, cds=cds, utr3=utr3)

            try:
                report = score_parsed(parsed)
                fitness = compute_fitness(report)
            except Exception:
                continue

            entry: Dict[str, Any] = {
                "label": f"candidate_{i+1}",
                "overall": float(fitness.get("overall")),
                "sequence": seq,
                "utr5": utr5,
            }
            for m in METRIC_NAMES:
                if m in report:
                    entry[m] = report[m]
            scored.append(entry)

        if not scored:
            payload["error"] = "No sequences could be scored."
            return payload

        scored.sort(key=lambda r: r["overall"], reverse=True)
        best = scored[0]

        # === THIS is what makes the dashboard show the plots ===
        # It writes PNGs into run_dir/plots and returns their file paths.
        plots = predict_and_plot_full_and_utrs(
            utr5=best["utr5"],
            cds=cds,
            utr3=utr3,
            out_dir=str(plot_dir),
            base_prefix=f"{gene}_{best['label']}",
            temperature_c=plot_temperature_c,
            show=False,
        )

        payload.update(
            {
                "ok": True,
                "error": "",
                "traceback": "",
                "lengths": {"cds": len(cds), "utr3": len(utr3)},
                "history": history,
                "best": best,
                "plots": plots,
            }
        )
        return payload

    except Exception as exc:
        payload["ok"] = False
        payload["error"] = f"Python error: {exc}"
        payload["traceback"] = traceback.format_exc()
        return payload