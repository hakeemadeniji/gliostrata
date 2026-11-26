"""L2 — Verhaak subtype classification.

Stub: real implementation is an ensemble of random forest + 1D-CNN trained on
TCGA-GBM labeled samples, with class probabilities returned alongside the argmax.

The stub `classify_subtypes` below reads the `subtype` column from the input
cohort (which `generate_synthetic_cohort` already populates) and adds calibrated
class-probability columns. This makes downstream code testable without the
trained classifier.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

SUBTYPES = ["Classical", "Mesenchymal", "Proneural", "Neural"]


def classify_subtypes(cohort: pd.DataFrame, *, seed: int = 42) -> pd.DataFrame:
    """Classify each sample into a Verhaak subtype.

    Returns a DataFrame with columns: patient_id, subtype, plus class probabilities
    for each subtype (subtype_prob_Classical, etc.).
    """
    rng = np.random.default_rng(seed)

    if "subtype" not in cohort.columns:
        raise ValueError("Input cohort missing 'subtype' column. "
                         "Use generate_synthetic_cohort() or populate manually.")

    n = len(cohort)
    # Stub: trust the input subtype, but add noisy class probabilities that peak there.
    probs = np.zeros((n, len(SUBTYPES)))
    for i, st in enumerate(cohort["subtype"].values):
        # Base distribution centered on the (already-known) subtype
        idx = SUBTYPES.index(st)
        base = np.full(4, 0.10)
        base[idx] = 0.70
        noisy = base + rng.normal(0, 0.05, size=4)
        noisy = np.clip(noisy, 0.01, None)
        probs[i] = noisy / noisy.sum()

    out = pd.DataFrame({
        "patient_id": cohort["patient_id"].values,
        "subtype": cohort["subtype"].values,
    })
    for j, st in enumerate(SUBTYPES):
        out[f"subtype_prob_{st}"] = probs[:, j]
    out["subtype_classifier_version"] = "0.1.0-stub"
    return out


def train_subtype_classifier(
    training_features: pd.DataFrame,
    training_labels: pd.Series,
) -> dict:
    """STUB — train the Verhaak subtype classifier ensemble.

    Real implementation:
      1. Build a feature matrix from Verhaak signature genes.
      2. Train random forest + 1D-CNN with 5-fold stratified CV.
      3. Save model weights + class label encoder to models/subtype_v1/.
      4. Return a metadata dict with CV metrics.

    Returns metadata dict (empty model) for now.
    """
    return {
        "model_version": "0.1.0-stub",
        "trained": False,
        "cv_accuracy": None,
        "n_training_samples": len(training_features),
        "n_features": training_features.shape[1],
        "note": "Stub — no actual training performed. See docs/roadmap.md.",
    }
