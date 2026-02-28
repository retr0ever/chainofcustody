"""Metric 4: Translation efficiency prediction via RiboNN (Sanofi).

RiboNN is a deep-CNN tool that predicts per-tissue translation efficiency from
mRNA sequences. It is bundled as a git submodule at ``vendor/RiboNN`` and its
Python dependencies are installed in the project venv.

Pretrained weights are downloaded automatically on first run.

Reference: Karollus et al., Nature Biotechnology (2024).
Repo: https://github.com/Sanofi-Public/RiboNN

Performance notes
-----------------
All 50 models (10 CV folds × 5 top models each) are loaded once into GPU
memory by :class:`RiboNNPredictor` and reused for every subsequent call.
Use :func:`score_ribonn_batch` during optimisation to score a whole
population in a single GPU pass; :func:`score_ribonn` is a thin wrapper for
single-sequence use (CLI / evaluation).
"""

from __future__ import annotations

import contextlib
import os
import sys
import tempfile
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import torch

from chainofcustody.sequence import mRNASequence
from chainofcustody.progress import update_status

if TYPE_CHECKING:
    pass

# Path to the RiboNN submodule, relative to this file:
#   chainofcustody/evaluation/ribonn.py  →  parents[2] = repo root
_RIBONN_DIR = Path(__file__).parents[2] / "vendor" / "RiboNN"

_SPECIES = "human"
_TOP_K = 5
_BATCH_SIZE = 512


def _ensure_importable() -> None:
    """Add vendor/RiboNN to sys.path so src.* modules can be imported."""
    ribonn_str = str(_RIBONN_DIR)
    if ribonn_str not in sys.path:
        sys.path.insert(0, ribonn_str)


def _te_status(mean_te: float) -> str:
    if mean_te >= 1.5:
        return "GREEN"
    if mean_te >= 1.0:
        return "AMBER"
    return "RED"


class RiboNNPredictor:
    """Holds all CV models in GPU memory for fast repeated inference.

    Use the module-level :func:`get_predictor` to obtain the shared singleton
    rather than instantiating this class directly.

    Parameters
    ----------
    ribonn_dir:
        Path to the ``vendor/RiboNN`` submodule root.
    species:
        ``"human"`` or ``"mouse"``.
    top_k:
        Number of top-ranked models to load per CV fold.
    """

    def __init__(
        self,
        ribonn_dir: Path = _RIBONN_DIR,
        species: str = _SPECIES,
        top_k: int = _TOP_K,
    ) -> None:
        _ensure_importable()
        from src.model import RiboNN  # noqa: PLC0415
        from src.utils.helpers import extract_config  # noqa: PLC0415

        self._ribonn_dir = ribonn_dir
        self._species = species
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        run_df = pd.read_csv(ribonn_dir / "models" / species / "runs.csv")
        self._config = extract_config(run_df, run_df.run_id[0])
        self._config["species"] = species
        self._config["max_utr5_len"] = 1_381
        self._config["max_cds_utr3_len"] = 11_937

        # Build list of (fold, model_path) pairs, sorted by fold
        run_df_sorted = run_df.sort_values("metrics.val_r2", ascending=False)
        self._fold_models: list[tuple[int, list[torch.nn.Module]]] = []
        n_folds = len(run_df["params.test_fold"].unique())
        for fold_idx, fold in enumerate(np.sort(run_df["params.test_fold"].unique()), 1):
            update_status(
                f"RiboNN  loading models  fold {fold_idx}/{n_folds}  "
                f"(device: {self.device})"
            )
            fold_str = str(fold)
            sub_df = run_df_sorted.query(
                "`params.test_fold` == @fold_str or `params.test_fold` == @fold"
            ).head(top_k)

            models = []
            with contextlib.chdir(ribonn_dir):
                for run_id in sub_df.run_id:
                    model = RiboNN(**self._config)
                    state = torch.load(
                        f"models/{species}/{run_id}/state_dict.pth",
                        map_location=self.device,
                    )
                    model.load_state_dict(state)
                    model.to(self.device)
                    model.eval()
                    models.append(model)

            self._fold_models.append((int(fold), models))

        # Derive predicted column names from one fold
        self._predicted_cols = self._get_predicted_cols()

    def _get_predicted_cols(self) -> list[str]:
        """Return the list of predicted_TE_* column names for this species."""
        _ensure_importable()
        from src.predict import predict_using_models_trained_in_one_fold  # noqa: PLC0415  # noqa: F401

        # Column names are derived from the species constant inside predict.py.
        # Re-derive them here without running inference.
        human_cols = (
            "TE_108T,TE_12T,TE_A2780,TE_A549,TE_BJ,TE_BRx.142,TE_C643,TE_CRL.1634,"
            "TE_Calu.3,TE_Cybrid_Cells,TE_H1.hESC,TE_H1933,TE_H9.hESC,TE_HAP.1,"
            "TE_HCC_tumor,TE_HCC_adjancent_normal,TE_HCT116,TE_HEK293,TE_HEK293T,"
            "TE_HMECs,TE_HSB2,TE_HSPCs,TE_HeLa,TE_HeLa_S3,TE_HepG2,TE_Huh.7.5,"
            "TE_Huh7,TE_K562,TE_Kidney_normal_tissue,TE_LCL,TE_LuCaP.PDX,TE_MCF10A,"
            "TE_MCF10A.ER.Src,TE_MCF7,TE_MD55A3,TE_MDA.MB.231,TE_MM1.S,TE_MOLM.13,"
            "TE_Molt.3,TE_Mutu,TE_OSCC,TE_PANC1,TE_PATU.8902,TE_PC3,TE_PC9,"
            "TE_Primary_CD4._T.cells,TE_Primary_human_bronchial_epithelial_cells,"
            "TE_RD.CCL.136,TE_RPE.1,TE_SH.SY5Y,TE_SUM159PT,TE_SW480TetOnAPC,"
            "TE_T47D,TE_THP.1,TE_U.251,TE_U.343,TE_U2392,TE_U2OS,TE_Vero_6,"
            "TE_WI38,TE_WM902B,TE_WTC.11,TE_ZR75.1,TE_cardiac_fibroblasts,TE_ccRCC,"
            "TE_early_neurons,TE_fibroblast,TE_hESC,TE_human_brain_tumor,"
            "TE_iPSC.differentiated_dopamine_neurons,TE_megakaryocytes,TE_muscle_tissue,"
            "TE_neuronal_precursor_cells,TE_neurons,TE_normal_brain_tissue,"
            "TE_normal_prostate,TE_primary_macrophages,TE_skeletal_muscle"
        )
        mouse_cols = (
            "TE_3T3,TE_4T1,TE_A3-1_mESC,TE_B16_melanoma_cell,"
            "TE_Bone_marrow_derived_macrophage,TE_Bone_marrow_derived_primary_dendritic_cells,"
            "TE_Bone_marrow_derived_regulatory_dendritic_cells,TE_CD4_T-cells,TE_CGR8_mESC,"
            "TE_Cerebellum,TE_DRG_neuronal_culture,TE_Dorsal_section_of_lumbar_spinal_cord,"
            "TE_E14_mESC,TE_E14Tg2a_mESC,TE_ES_cell_derived_neurons,TE_Epiblast_like_cells,"
            "TE_Epidermal_basal_cells,TE_Fetal_cortex,TE_Flt3L-DC,TE_Forebrain,TE_HFSC,"
            "TE_Ileum,TE_J1_mESC,TE_KH2_mESC,"
            "TE_KRPC-A_cells_(from_gen._eng._murine_tumors),TE_MN1_cells,TE_Mammary_tissue,"
            "TE_Mouse_back_skins,TE_NSC,TE_Neuro2a,"
            "TE_Neurons_(DIV_8)_derived_from_CGR8_ES_cells,TE_P19,TE_Primary_Keratinocytes,"
            "TE_Primary_cortical_neurons,TE_Quadriceps_muscle,TE_R1_mESC,TE_R1/E_mESC,"
            "TE_RAW264.7,TE_RF8_mESC,"
            "TE_Spontaneously_Immortalized_Mouse_Keratinocyte_Cult,TE_Striatal_cells,"
            "TE_T-ALL,TE_brain,TE_dentate_gyrus,TE_duodenum,TE_embryonic_fibroblast,"
            "TE_embryonic_stem_cells,TE_epididymal_white_fat,TE_forelimbs,"
            "TE_gastrocnemius_tissue,TE_heart,TE_hippocampal,TE_interscapular_brown_fat,"
            "TE_kidney,TE_liver,TE_lung,TE_lymphoid_ba/f3_cells,TE_mouse_eye,"
            "TE_neural_tube,TE_neutrophils,TE_pancreas,TE_resting_state_T_cells,"
            "TE_skeletal_muscle,TE_skin_squamous_tumours_(skin_papilloma),"
            "TE_subcutaneous_white_fat,TE_testis,TE_tibialis_anterior_tissue,TE_v6.5_mESC"
        )
        cols = human_cols if self._species == "human" else mouse_cols
        return [c.replace("TE_", "predicted_TE_") for c in cols.split(",")]

    def predict_batch(self, sequences: list[mRNASequence]) -> list[dict]:
        """Score a batch of sequences.  Returns one result dict per sequence.

        Each dict has keys ``mean_te``, ``per_tissue``, ``status``, ``message``.
        Sequences that exceed RiboNN's length limits are returned with
        ``mean_te=0.0`` and ``status="RED"``.
        """
        _ensure_importable()
        from src.data import RiboNNDataModule  # noqa: PLC0415

        rows = [
            {
                "tx_id": str(i),
                "utr5_sequence": seq.utr5,
                "cds_sequence": seq.cds,
                "utr3_sequence": seq.utr3,
            }
            for i, seq in enumerate(sequences)
        ]
        input_df = pd.DataFrame(rows)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=True
        ) as fh:
            input_df.to_csv(fh, sep="\t", index=False)
            fh.flush()
            input_path = fh.name

            config = dict(self._config)
            config["tx_info_path"] = input_path
            config["num_workers"] = 0
            config["test_batch_size"] = _BATCH_SIZE
            config["remove_extreme_txs"] = False
            config["target_column_pattern"] = None

            with contextlib.chdir(self._ribonn_dir), \
                 warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="'pin_memory' argument is set as true but not supported on MPS",
                    category=UserWarning,
                )
                dm = RiboNNDataModule(config)
                dataloader = dm.predict_dataloader()

                all_fold_preds: list[np.ndarray] = []
                for _fold, models in self._fold_models:
                    fold_model_preds: list[np.ndarray] = []
                    for model in models:
                        with torch.no_grad():
                            batches = [
                                model(batch.to(self.device))
                                for batch in dataloader
                            ]
                        preds = torch.cat(batches, dim=0).cpu().numpy()
                        fold_model_preds.append(preds)
                    # Average across the top-k models in this fold
                    all_fold_preds.append(np.stack(fold_model_preds).mean(axis=0))

        # Average across folds: shape (n_sequences, n_tissues)
        mean_preds = np.stack(all_fold_preds).mean(axis=0)

        # Map back — dm.df may have been filtered; align by tx_id
        valid_ids = set(dm.df["tx_id"].astype(str))
        pred_idx = 0
        results: list[dict] = []
        for i in range(len(sequences)):
            if str(i) not in valid_ids:
                results.append(_null_result())
                continue
            tissue_preds = mean_preds[pred_idx]
            pred_idx += 1
            mean_te = float(tissue_preds.mean())
            per_tissue = {
                col.removeprefix("predicted_TE_"): round(float(v), 4)
                for col, v in zip(self._predicted_cols, tissue_preds)
            }
            results.append({
                "mean_te": round(mean_te, 4),
                "per_tissue": per_tissue or None,
                "status": _te_status(mean_te),
                "message": f"RiboNN predicted mean TE = {mean_te:.4f}",
            })

        return results


# ---------------------------------------------------------------------------
# Module-level singleton
# ---------------------------------------------------------------------------

_predictor: RiboNNPredictor | None = None


def get_predictor(ribonn_dir: Path = _RIBONN_DIR) -> RiboNNPredictor:
    """Return (or create) the module-level :class:`RiboNNPredictor` singleton."""
    global _predictor
    if _predictor is None:
        _predictor = RiboNNPredictor(ribonn_dir=ribonn_dir)
    return _predictor


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def score_ribonn(parsed: mRNASequence) -> dict:
    """Predict translation efficiency for a single sequence using RiboNN.

    Loads models on first call; subsequent calls reuse cached GPU models.

    Returns a dict with keys: ``mean_te``, ``per_tissue``, ``status``,
    ``message``.
    """
    return get_predictor().predict_batch([parsed])[0]


def score_ribonn_batch(sequences: list[mRNASequence]) -> list[dict]:
    """Predict translation efficiency for a list of sequences in one GPU pass.

    Significantly faster than calling :func:`score_ribonn` in a loop because
    the whole population is forwarded through the cached GPU models together.

    Returns one result dict per input sequence (same order).
    """
    return get_predictor().predict_batch(sequences)


def _null_result() -> dict:
    return {
        "mean_te": 0.0,
        "per_tissue": None,
        "status": "RED",
        "message": "Sequence excluded (exceeds RiboNN length limits)",
    }
