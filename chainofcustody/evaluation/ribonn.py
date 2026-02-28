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

The hot path in :meth:`RiboNNPredictor.predict_batch` bypasses the
``RiboNNDataModule`` / ``DataFrameDataset`` entirely.  The original
``DataFrameDataset.__getitem__`` builds one-hot tensors with Python for-loops
over individual nucleotide characters — O(sequence_length) Python iterations
per sample, completely unvectorized.  Instead we use
:func:`_encode_sequences_vectorized` which encodes the whole batch in one
numpy operation on a pre-allocated ``(N, C, L)`` array, then pins and
transfers it to GPU as a single ``torch.Tensor``.

Use :func:`score_ribonn_batch` during optimisation to score a whole
population in a single GPU pass; :func:`score_ribonn` is a thin wrapper for
single-sequence use (CLI / evaluation).
"""

from __future__ import annotations

import contextlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from chainofcustody.sequence import mRNASequence
from chainofcustody.progress import update_status

# Path to the RiboNN submodule, relative to this file:
#   chainofcustody/evaluation/ribonn.py  →  parents[2] = repo root
_RIBONN_DIR = Path(__file__).parents[2] / "vendor" / "RiboNN"

_SPECIES = "human"
_TOP_K = 5

# Sequence length limits used when training the human model (from config)
_MAX_UTR5_LEN = 1_381
_MAX_CDS_UTR3_LEN = 11_937
_PADDED_LEN = _MAX_UTR5_LEN + _MAX_CDS_UTR3_LEN  # 13 318

# RiboNN uses DNA internally (U→T at ingestion).
# Channel layout for the human model (pad_5_prime=True,
# split_utr5_cds_utr3_channels=False, label_codons=True, all others False):
#   ch 0-3: one-hot A/T/C/G, right-padded so start codon aligns at _MAX_UTR5_LEN
#   ch 4:   codon-start mask (1 at first nt of every CDS codon)
_N_CHANNELS = 5

# Map nucleotide → channel index (DNA alphabet; U treated as T)
_NT_INDEX: dict[str, int] = {"A": 0, "T": 1, "U": 1, "C": 2, "G": 3}

# Pre-built 256-entry LUT: ASCII byte value → channel index (0 for unknowns).
# Built once at import time; reused for every sequence in every epoch.
_NT_LUT: np.ndarray = np.zeros(256, dtype=np.int8)
for _ch, _idx in _NT_INDEX.items():
    _NT_LUT[ord(_ch)] = _idx


def _ensure_importable() -> None:
    """Add vendor/RiboNN to sys.path so src.* modules can be imported."""
    ribonn_str = str(_RIBONN_DIR)
    if ribonn_str not in sys.path:
        sys.path.insert(0, ribonn_str)


_HUMAN_TISSUE_NAMES: list[str] = [
    c.strip().removeprefix("TE_")
    for c in (
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
    ).split(",")
]


def get_valid_tissue_names(species: str = "human") -> list[str]:
    """Return the tissue names available for *species* without loading model weights."""
    if species != "human":
        raise ValueError(f"Only 'human' tissue names are available statically; got {species!r}")
    return list(_HUMAN_TISSUE_NAMES)


def _te_status(target_te: float, mean_off_target_te: float) -> str:
    selectivity = target_te - mean_off_target_te
    if target_te >= 1.5 and selectivity >= 0.5:
        return "GREEN"
    if target_te >= 1.0 and selectivity >= 0.0:
        return "AMBER"
    return "RED"


def _encode_sequences_vectorized(
    sequences: list[mRNASequence],
) -> tuple[torch.Tensor, list[bool]]:
    """Encode a batch of mRNA sequences into the RiboNN input tensor format.

    Builds a ``(N, 5, 13318)`` float32 tensor on CPU using vectorized numpy
    operations, then returns it as a pinned-memory tensor ready for GPU
    transfer.  Sequences that exceed the model's length limits are skipped;
    their slot is left as zeros and ``valid[i]`` is set to ``False``.

    Returns
    -------
    tensor : torch.Tensor  shape (N, 5, 13318), pinned
    valid  : list[bool]     True for sequences that were encoded
    """
    n = len(sequences)
    # Pre-allocate on CPU as a numpy array; fill in-place per sequence
    arr = np.zeros((n, _N_CHANNELS, _PADDED_LEN), dtype=np.float32)
    valid = [False] * n

    for i, seq in enumerate(sequences):
        # RiboNN expects DNA (U→T already handled via _NT_INDEX mapping)
        utr5 = seq.utr5
        cds = seq.cds
        utr3 = seq.utr3

        utr5_len = len(utr5)
        cds_len = len(cds)
        cds_utr3_len = cds_len + len(utr3)
        tx_len = utr5_len + cds_utr3_len

        if utr5_len > _MAX_UTR5_LEN or cds_utr3_len > _MAX_CDS_UTR3_LEN:
            continue  # leave valid[i] = False

        # Full transcript as a numpy char array
        tx = np.frombuffer((utr5 + cds + utr3).encode(), dtype=np.uint8)

        # One-hot: map each character to its channel index via a lookup table
        nt_channels = _NT_LUT[tx]  # shape (tx_len,)

        # Position in the padded tensor: UTR5 is right-aligned to _MAX_UTR5_LEN
        pad_offset = _MAX_UTR5_LEN - utr5_len  # start position in padded axis
        positions = np.arange(tx_len, dtype=np.int32) + pad_offset

        # Scatter one-hot values using advanced indexing
        arr[i, nt_channels, positions] = 1.0

        # Codon-start mask (channel 4): every 3rd position starting at CDS start
        cds_start = _MAX_UTR5_LEN  # aligned after padding
        codon_positions = np.arange(cds_start, cds_start + cds_len - 3 + 1, 3)
        arr[i, 4, codon_positions] = 1.0

        valid[i] = True

    # Pin memory for faster CPU→GPU transfer
    tensor = torch.from_numpy(arr).pin_memory()
    return tensor, valid


class RiboNNPredictor:
    """Holds all CV models in GPU memory for fast repeated inference.

    Use the module-level :func:`get_predictor` to obtain the shared singleton
    rather than instantiating this class directly.
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

        # Enable TF32 on Ampere+ GPUs — uses tensor cores for matmuls at no API cost.
        torch.set_float32_matmul_precision("high")

        run_df = pd.read_csv(ribonn_dir / "models" / species / "runs.csv")
        config = extract_config(run_df, run_df.run_id[0])
        config["species"] = species
        config["max_utr5_len"] = _MAX_UTR5_LEN
        config["max_cds_utr3_len"] = _MAX_CDS_UTR3_LEN

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
                    model = RiboNN(**config)
                    state = torch.load(
                        f"models/{species}/{run_id}/state_dict.pth",
                        map_location=self.device,
                    )
                    model.load_state_dict(state)
                    model.to(self.device)
                    model.eval()
                    models.append(model)

            self._fold_models.append((int(fold), models))

        self._predicted_cols = self._get_predicted_cols()
        update_status("RiboNN  ready")

    def _get_predicted_cols(self) -> list[str]:
        if self._species == "human":
            names = _HUMAN_TISSUE_NAMES
        else:
            mouse_names = [
                c.strip().removeprefix("TE_")
                for c in (
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
                ).split(",")
            ]
            names = mouse_names
        return [f"predicted_TE_{n}" for n in names]

    def predict_batch(
        self,
        sequences: list[mRNASequence],
        target_cell_type: str = "megakaryocytes",
    ) -> list[dict]:
        """Score a batch of sequences in one GPU pass.

        Encodes the whole batch with vectorized numpy (bypassing the slow
        per-character Python loop in ``DataFrameDataset.__getitem__``), then
        runs a single forward pass per model with the full batch on GPU.

        Args:
            sequences: Batch of sequences to score.
            target_cell_type: Tissue name that should have *high* TE.  All
                other tissues are treated as off-target (low TE desired).
                Must be a key in the per-tissue prediction dict.

        Returns one result dict per input sequence.
        """
        n = len(sequences)

        # Vectorized CPU encoding → pinned tensor
        batch_tensor, valid = _encode_sequences_vectorized(sequences)
        # Move the whole batch to GPU once
        batch_gpu = batch_tensor.to(self.device, non_blocking=True)

        # --- Run all 50 models (10 folds × 5 top-k) ---
        all_fold_preds: list[np.ndarray] = []
        for _fold, models in self._fold_models:
            fold_model_preds: list[np.ndarray] = []
            for model in models:
                with torch.no_grad():
                    out = model(batch_gpu).cpu().numpy()  # (N, n_tissues)
                fold_model_preds.append(out)
            all_fold_preds.append(np.stack(fold_model_preds).mean(axis=0))

        # Average across folds: (N, n_tissues)
        mean_preds = np.stack(all_fold_preds).mean(axis=0)

        # Resolve the index of the target tissue once for the whole batch
        tissue_names = [col.removeprefix("predicted_TE_") for col in self._predicted_cols]
        if target_cell_type not in tissue_names:
            raise ValueError(
                f"Unknown target cell type {target_cell_type!r}. "
                f"Valid names: {', '.join(sorted(tissue_names))}"
            )
        target_idx = tissue_names.index(target_cell_type)

        results: list[dict] = []
        for i in range(n):
            if not valid[i]:
                results.append(_null_result(target_cell_type))
                continue
            tissue_preds = mean_preds[i]
            mean_te = float(tissue_preds.mean())
            target_te = float(tissue_preds[target_idx])
            off_target_mask = np.ones(len(tissue_preds), dtype=bool)
            off_target_mask[target_idx] = False
            mean_off_target_te = float(tissue_preds[off_target_mask].mean())
            per_tissue = {
                name: round(float(v), 4)
                for name, v in zip(tissue_names, tissue_preds)
            }
            results.append({
                "mean_te": round(mean_te, 4),
                "target_cell_type": target_cell_type,
                "target_te": round(target_te, 4),
                "mean_off_target_te": round(mean_off_target_te, 4),
                "per_tissue": per_tissue,
                "status": _te_status(target_te, mean_off_target_te),
                "message": (
                    f"RiboNN: {target_cell_type} TE = {target_te:.4f}, "
                    f"mean off-target = {mean_off_target_te:.4f}"
                ),
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

def score_ribonn(parsed: mRNASequence, target_cell_type: str = "megakaryocytes") -> dict:
    """Predict translation efficiency for a single sequence using RiboNN.

    Loads models on first call; subsequent calls reuse cached GPU models.
    """
    return get_predictor().predict_batch([parsed], target_cell_type=target_cell_type)[0]


def score_ribonn_batch(
    sequences: list[mRNASequence],
    target_cell_type: str = "megakaryocytes",
) -> list[dict]:
    """Predict translation efficiency for a list of sequences in one GPU pass.

    Significantly faster than calling :func:`score_ribonn` in a loop because
    the whole population is encoded and forwarded in a single vectorized
    operation.
    """
    return get_predictor().predict_batch(sequences, target_cell_type=target_cell_type)


def _null_result(target_cell_type: str = "megakaryocytes") -> dict:
    return {
        "mean_te": 0.0,
        "target_cell_type": target_cell_type,
        "target_te": 0.0,
        "mean_off_target_te": 0.0,
        "per_tissue": None,
        "status": "RED",
        "message": "Sequence excluded (exceeds RiboNN length limits)",
    }
