# Evaluation Metrics

Four metrics are currently implemented in `chainofcustody/evaluation/`. All four contribute to a weighted overall fitness score (0–1, higher is better).

## Metric Overview

| Metric | Key | Weight | Module | GREEN threshold |
|---|---|---:|---|---|
| 5'UTR Accessibility | `utr5_accessibility` | 22% | `structure.py` | MFE < −30 kcal/mol |
| Manufacturability | `manufacturability` | 30% | `manufacturing.py` | 0 violations |
| Stability | `stability` | 35% | `stability.py` | Combined score ≥ 0.7 |
| Translation Efficiency | `translation_efficiency` | 13% | `ribonn.py` | Mean TE ≥ 1.5 |

---

## 1. 5'UTR Accessibility (`structure.py`)

**What it measures:** Whether the 5'UTR is accessible for ribosome loading. Strong secondary structure in the 5'UTR blocks the 43S pre-initiation complex.

**How it works:**
- Folds the 5'UTR + first 30 nt of CDS (the ribosome landing zone) with ViennaRNA
- For long sequences (> 2000 nt), uses a windowed fold (500 nt windows, 250 nt step) to avoid O(n³) blowup

**Traffic light:**
| Status | Condition |
|---|---|
| GREEN | MFE < −30 kcal/mol |
| AMBER | −30 ≤ MFE < −20 kcal/mol |
| RED | MFE ≥ −20 kcal/mol |

**Normalisation (0→1):** Linear from 0 at MFE = 0 to 1.0 at MFE ≤ −30.

**Additional outputs (not in fitness score):**
- `global_mfe` — full-sequence MFE and MFE/nt
- `mirna_site_accessibility` — per-site seed region pairing (only when site positions are provided)

---

## 2. Manufacturability (`manufacturing.py`)

**What it measures:** Whether the sequence can be reliably synthesised as DNA. Checks three independent failure modes.

### Sub-checks

| Check | Rule | Violates if… |
|---|---|---|
| GC windows | Sliding 50 nt window, GC must be 30–70% | Any window falls outside range |
| Homopolymers | Max run length 8 nt | Any single-nucleotide run > 8 nt |
| Restriction sites | 6 enzymes on forward + reverse strand | Recognition sequence present anywhere |

**Restriction enzymes checked:** BsaI (`GGUCUC`), BsmBI (`CGUCUC`), EcoRI (`GAAUUC`), BamHI (`GGAUCC`), HindIII (`AAGCUU`), NotI (`GCGGCCGC`).

**Traffic light:** Derived from violation count:
| Status | Condition |
|---|---|
| GREEN | 0 violations |
| AMBER | 1–3 violations |
| RED | > 3 violations |

**Normalisation (0→1):** `max(0, 1 − violations / 10)` — reaches 0 at 10+ violations.

---

## 3. Stability (`stability.py`)

**What it measures:** mRNA half-life in vivo. Combines three independent sub-metrics.

### Sub-metrics

| Sub-metric | Function | What it captures |
|---|---|---|
| GC3 wobble content | `compute_gc3` | G/C at 3rd codon position; higher GC3 → longer half-life in mammals |
| MFE per nucleotide | `compute_mfe_per_nt` | Thermodynamic stability of full sequence (kcal/mol/nt); more negative = more stable |
| AU-rich elements (AREs) | `count_au_rich_elements` | Count of `AUUUA` pentamers in 3'UTR; each recruits exosome-mediated degradation |

**Normalisation of sub-metrics:**
- **GC3:** 1.0 in the 0.5–0.7 optimal range; linearly penalised outside
- **MFE/nt:** 1.0 at ≤ −0.4 kcal/mol/nt; 0.0 at ≥ −0.1 kcal/mol/nt
- **ARE count:** `max(0, 1 − count / 3)` — 0 at 3+ elements

**Combined stability score:** `0.4 × GC3_norm + 0.4 × MFE/nt_norm + 0.2 × ARE_norm`

**Traffic light:**
| Status | Condition |
|---|---|
| GREEN | Combined score ≥ 0.7 |
| AMBER | 0.4 ≤ score < 0.7 |
| RED | score < 0.4 |

---

## 4. Translation Efficiency (`ribonn.py`)

**What it measures:** Predicted ribosome loading and translation rate across human tissues, using the [RiboNN deep-CNN model](https://github.com/Sanofi-Public/RiboNN) (Karollus et al., *Nature Biotechnology* 2024).

**How it works:**
- Writes the split sequence (5'UTR / CDS / 3'UTR) to RiboNN's tab-separated input format
- Invokes `python -m src.main --predict human` as a subprocess using the project venv
- Parses `results/human/prediction_output.txt` for `mean_predicted_TE` and per-tissue predictions

**Setup:** The RiboNN repository is included as a git submodule at `vendor/RiboNN`. After cloning, initialise it with:
```bash
git submodule update --init vendor/RiboNN
```
Pretrained weights are downloaded automatically on first run.

**Traffic light:**
| Status | Condition |
|---|---|
| GREEN | Mean TE ≥ 1.5 |
| AMBER | 1.0 ≤ mean TE < 1.5 |
| RED | Mean TE < 1.0 |
| GREY | RiboNN unavailable |

**Normalisation (0→1):** Linear from 0 at mean TE ≤ 0.5 to 1.0 at mean TE ≥ 2.0.

---

## Fitness Aggregation (`fitness.py`)

`compute_fitness(report)` returns:
```python
{
    "scores": {
        "<metric>": {"value": float, "weight": float, "weighted": float, "status": str}
    },
    "overall": float,          # weighted sum, 0–1
    "suggestions": [           # only for non-GREEN metrics, sorted high→low priority
        {"metric": str, "priority": "high"|"medium"|"low", "action": str}
    ]
}
```

Priority is derived from metric weight: ≥ 35% → high, ≥ 25% → medium, < 25% → low.  
Currently: Stability = high, Manufacturability = medium, 5'UTR and TE = low.

---

## Not yet implemented

- **Codon optimisation** (CAI, rare codon clusters) — removed; would go in a new `codons.py`
- **miRNA detargeting** (miR-122 liver silencing) — removed; would go in a new `mirna.py`
- **5'UTR generation** — no builder exists yet
- **Full mRNA assembly** (5'UTR + Kozak + CDS + 3'UTR concatenation) — handled by the optimiser but not as a standalone scorer
