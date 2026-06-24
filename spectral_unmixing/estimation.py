"""
Alpha estimation helpers for spectral bleed-through correction.

Author: Fabrizio Musacchio
Date: June 2026
"""

from __future__ import annotations

import numpy as np

MIN_MASK_VOXELS = 16


def _validate_percentile(name: str, value: float) -> float:
    value = float(value)
    if not 0.0 <= value <= 100.0:
        raise ValueError(f"{name} must be between 0 and 100. Got {value!r}.")
    return value


def estimate_alpha_from_volume(
    source,
    target,
    signal_percentile: float = 99.0,
    background_percentile: float = 1.0,
    min_mask_voxels: int = MIN_MASK_VOXELS,
) -> float:
    """
    Estimate bleed-through coefficient alpha from matching source and target volumes.

    Parameters
    ----------
    source, target : array-like
        Matching source and target volumes, typically with shape ``(Z, Y, X)``.
    signal_percentile : float, optional
        Percentile used to define the bright-source mask.
    background_percentile : float, optional
        Percentile used as a rough per-channel background estimate before masking.
    min_mask_voxels : int, optional
        Minimum number of voxels required in the source-signal mask.

    Returns
    -------
    float
        Estimated scalar bleed-through coefficient.

    Raises
    ------
    ValueError
        If the input shapes differ, the signal mask is too small, or the estimate
        is numerically invalid.
    """

    signal_percentile = _validate_percentile("signal_percentile", signal_percentile)
    background_percentile = _validate_percentile(
        "background_percentile", background_percentile
    )

    source_f = np.asarray(source, dtype=np.float32)
    target_f = np.asarray(target, dtype=np.float32)

    if source_f.shape != target_f.shape:
        raise ValueError(
            "source and target must have the same shape. "
            f"Got {source_f.shape!r} and {target_f.shape!r}."
        )
    if source_f.size == 0:
        raise ValueError("source and target must not be empty.")
    if min_mask_voxels < 1:
        raise ValueError(f"min_mask_voxels must be >= 1. Got {min_mask_voxels!r}.")

    source_bg = np.percentile(source_f, background_percentile)
    target_bg = np.percentile(target_f, background_percentile)

    source_corrected = np.clip(source_f - source_bg, a_min=0.0, a_max=None)
    target_corrected = np.clip(target_f - target_bg, a_min=0.0, a_max=None)

    signal_threshold = np.percentile(source_corrected, signal_percentile)
    mask = source_corrected > signal_threshold

    mask_voxels = int(np.count_nonzero(mask))
    if mask_voxels < min_mask_voxels:
        raise ValueError(
            "Source signal mask does not contain enough voxels for alpha estimation. "
            f"Found {mask_voxels}, need at least {min_mask_voxels}."
        )

    denominator = float(np.mean(source_corrected[mask]))
    if denominator <= 0.0:
        raise ValueError(
            "Mean source intensity inside the source-signal mask must be > 0."
        )

    numerator = float(np.mean(target_corrected[mask]))
    alpha = numerator / denominator

    if not np.isfinite(alpha):
        raise ValueError(f"Estimated alpha is not finite: {alpha!r}.")

    return float(alpha)
