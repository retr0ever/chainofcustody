"""Gradient-based 5'UTR sequence design via RiboNN backpropagation.

Uses the RiboNN ensemble as a differentiable oracle: continuous nucleotide
logits are optimised by gradient ascent on the target-tissue TE prediction,
then discretized to produce warm-start chromosome rows for NSGA-III.

The entire approach is cell-type-agnostic — the only cell-type-specific input
is the index of the target tissue column in the RiboNN output, which is
already resolved by :class:`~chainofcustody.evaluation.ribonn.RiboNNPredictor`.

Algorithm
---------
1. Initialise a ``(utr5_len, 4)`` logit tensor as a ``nn.Parameter``.
2. For each gradient step:
   a. Apply ``softmax`` over the nucleotide axis → soft one-hot probabilities.
   b. Construct the full ``(1, 5, 13318)`` RiboNN input tensor by replacing
      the 5'UTR channels with the soft embedding (the CDS/3'UTR channels are
      fixed one-hot, as in normal inference).
   c. Average the forward pass across all 50 ensemble models.
   d. Backpropagate the negative target-tissue TE (minimise → gradient ascent).
3. Discretize: ``argmax`` over the nucleotide axis gives the best sequence.
4. Repeat from a fresh random start *n_restarts* times; return the top-*n_seeds*
   sequences ranked by their discretized TE prediction.
"""

from __future__ import annotations

import logging

import numpy as np
import torch
import torch.nn as nn

from chainofcustody.sequence import mRNASequence, KOZAK
from chainofcustody.evaluation.ribonn import get_predictor, score_ribonn

logger = logging.getLogger(__name__)

# ── Constants matching RiboNN's encoding (from evaluation/ribonn.py) ──────────
_MAX_UTR5_LEN = 1_381
_MAX_CDS_UTR3_LEN = 11_937
_PADDED_LEN = _MAX_UTR5_LEN + _MAX_CDS_UTR3_LEN  # 13 318
_N_CHANNELS = 5
_NT_INDEX: dict[str, int] = {"A": 0, "T": 1, "U": 1, "C": 2, "G": 3}

# Nucleotide encoding used by the optimizer (A=0, C=1, G=2, U=3)
_NUCLEOTIDES = np.array(["A", "C", "G", "U"])
# Map optimizer index → RiboNN channel index (DNA)
#   optimizer: A=0, C=1, G=2, U=3
#   ribonn:    A=0, T/U=1, C=2, G=3
_OPT_TO_RIBONN = [0, 2, 3, 1]  # A→0, C→2, G→3, U→1


def _build_fixed_cds_utr3_tensor(
    cds: str,
    utr3: str,
    device: torch.device,
) -> torch.Tensor:
    """Build the fixed (CDS + 3'UTR) portion of the RiboNN input.

    Returns a ``(1, 5, 13318)`` float32 tensor with:
    - Channels 0-3: one-hot nucleotide encoding for CDS+3'UTR at positions
      _MAX_UTR5_LEN.._PADDED_LEN (5'UTR region stays zero).
    - Channel 4: codon-start mask at every 3rd CDS position.
    """
    arr = np.zeros((1, _N_CHANNELS, _PADDED_LEN), dtype=np.float32)

    cds_utr3 = (cds + utr3).upper()
    cds_len = len(cds)
    cds_utr3_len = len(cds_utr3)

    if cds_utr3_len > _MAX_CDS_UTR3_LEN:
        cds_utr3 = cds_utr3[:_MAX_CDS_UTR3_LEN]
        cds_utr3_len = _MAX_CDS_UTR3_LEN

    tx_bytes = np.frombuffer(cds_utr3.encode(), dtype=np.uint8)
    lut = np.zeros(256, dtype=np.int8)
    for ch, idx in _NT_INDEX.items():
        lut[ord(ch)] = idx

    nt_channels = lut[tx_bytes]
    positions = np.arange(cds_utr3_len, dtype=np.int32) + _MAX_UTR5_LEN
    arr[0, nt_channels, positions] = 1.0

    codon_positions = np.arange(_MAX_UTR5_LEN, _MAX_UTR5_LEN + cds_len - 3 + 1, 3)
    arr[0, 4, codon_positions] = 1.0

    return torch.from_numpy(arr).to(device)


def _soft_utr5_to_ribonn_input(
    logits: torch.Tensor,          # (utr5_len, 4) — nucleotide logits in optimizer order
    fixed_cds_utr3: torch.Tensor,  # (1, 5, 13318) — CDS+3'UTR channels, no grad
    utr5_len: int,
) -> torch.Tensor:
    """Combine soft 5'UTR probabilities with fixed CDS/3'UTR to form full RiboNN input.

    The 5'UTR is right-aligned: position ``_MAX_UTR5_LEN - utr5_len`` is the
    first nucleotide of the 5'UTR (matching the hard-encoding in ribonn.py).

    Returns: ``(1, 5, 13318)`` float32 tensor with gradients through *logits*.
    """
    probs = torch.softmax(logits, dim=-1)  # (utr5_len, 4) in optimizer order

    # Re-order from optimizer (A=0,C=1,G=2,U=3) to RiboNN (A=0,T/U=1,C=2,G=3)
    reorder_idx = torch.tensor(_OPT_TO_RIBONN, device=logits.device)
    probs_ribonn = probs[:, reorder_idx]  # (utr5_len, 4)

    # Build a (1, 4, utr5_len) slice and place it into a copy of the fixed tensor
    utr5_channels = probs_ribonn.T.unsqueeze(0)  # (1, 4, utr5_len)

    # Clone fixed tensor and fill 5'UTR channels (no in-place on grad tensors)
    x = fixed_cds_utr3.clone()
    pad_start = _MAX_UTR5_LEN - utr5_len
    x[0, :4, pad_start:pad_start + utr5_len] = utr5_channels[0]

    return x


def _run_ensemble(
    x: torch.Tensor,
    fold_models: list[tuple[int, list[nn.Module]]],
) -> torch.Tensor:
    """Average forward pass over all ensemble models. Returns ``(1, n_tissues)``."""
    preds = []
    for _fold, models in fold_models:
        fold_preds = [model(x) for model in models]
        preds.append(torch.stack(fold_preds).mean(0))
    return torch.stack(preds).mean(0)


def generate_gradient_seeds(
    cds: str,
    utr3: str,
    target_cell_type: str,
    utr5_len: int = 100,
    n_steps: int = 200,
    n_seeds: int = 16,
    n_restarts: int = 4,
    lr: float = 0.05,
    utr5_max: int = 1000,
) -> list[np.ndarray]:
    """Design high-TE 5'UTR sequences by gradient ascent through RiboNN.

    Runs *n_restarts* independent gradient-ascent trajectories, each starting
    from a fresh random logit initialisation.  Returns up to *n_seeds*
    chromosome rows (``np.ndarray`` of shape ``utr5_max + 1``) ranked by their
    discretized TE prediction for *target_cell_type*.

    Args:
        cds: Fixed CDS sequence (RNA, uppercase).
        utr3: Fixed 3'UTR sequence (RNA, uppercase).
        target_cell_type: RiboNN tissue column name (e.g. ``"neurons"``).
        utr5_len: Length of the 5'UTR to design.  Clipped to [1, 1381].
        n_steps: Number of gradient-ascent steps per restart.
        n_seeds: Maximum number of chromosome rows to return.
        n_restarts: Number of independent restarts; best sequences are kept.
        lr: Learning rate for Adam.
        utr5_max: Maximum 5'UTR length supported by the optimizer problem
            (determines chromosome row length = utr5_max + 1).

    Returns:
        List of chromosome rows as integer ``np.ndarray`` of shape
        ``(utr5_max + 1,)``, sorted best-TE-first.  Empty list if RiboNN is
        not available.
    """
    utr5_len = int(np.clip(utr5_len, 1, min(_MAX_UTR5_LEN, utr5_max)))

    try:
        predictor = get_predictor()
    except Exception as exc:
        logger.warning("Could not load RiboNN predictor: %s", exc)
        return []

    device = predictor.device
    fold_models = predictor._fold_models

    # Resolve target tissue index
    tissue_names = [col.removeprefix("predicted_TE_") for col in predictor._predicted_cols]
    if target_cell_type not in tissue_names:
        logger.warning(
            "gradient_seed: unknown target_cell_type %r; valid: %s",
            target_cell_type,
            ", ".join(tissue_names),
        )
        return []
    target_idx = tissue_names.index(target_cell_type)

    # Pre-build fixed CDS+3'UTR tensor once (no grad needed)
    fixed = _build_fixed_cds_utr3_tensor(cds, utr3, device)
    fixed.requires_grad_(False)

    results: list[tuple[float, np.ndarray]] = []

    for restart in range(n_restarts):
        logits = nn.Parameter(torch.randn(utr5_len, 4, device=device))
        optimizer_gd = torch.optim.Adam([logits], lr=lr)

        for step in range(n_steps):
            optimizer_gd.zero_grad()
            x = _soft_utr5_to_ribonn_input(logits, fixed, utr5_len)
            pred = _run_ensemble(x, fold_models)  # (1, n_tissues)
            loss = -pred[0, target_idx]           # maximise target TE
            loss.backward()
            optimizer_gd.step()

        # Discretize: argmax over nucleotide axis → sequence in optimizer encoding
        with torch.no_grad():
            best_nt_indices = logits.argmax(dim=-1).cpu().numpy()  # (utr5_len,) in opt order

        # Evaluate the discretized sequence to get a comparable TE score
        utr5_str = "".join(_NUCLEOTIDES[best_nt_indices])
        seq = mRNASequence(utr5=utr5_str + KOZAK, cds=cds, utr3=utr3)
        try:
            ribonn = score_ribonn(seq, target_cell_type=target_cell_type)
            te = ribonn.get("target_te", 0.0)
        except Exception:
            te = float(-pred[0, target_idx].item())  # fallback: use soft TE

        # Build chromosome row
        row = np.zeros(utr5_max + 1, dtype=int)
        row[0] = utr5_len
        row[1:utr5_len + 1] = best_nt_indices
        results.append((te, row))

    results.sort(key=lambda x: x[0], reverse=True)
    return [row for _, row in results[:n_seeds]]
