"""
Alpha estimation helpers for spectral bleed-through correction.

Author: Fabrizio Musacchio
Date: June 2026
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy.optimize import minimize_scalar

MIN_MASK_VOXELS = 16
DEFAULT_ALPHA_MAX = 1.0
DEFAULT_MI_BINS = 64
DEFAULT_MAX_ALPHA_VOXELS = 500_000
DEFAULT_RANDOM_STATE = 0
EPSILON = 1e-12

SUPPORTED_ALPHA_ESTIMATION_METHODS = {
    "mean_ratio",
    "linear_fit",
    "corr_min",
    "mi_min",
}


def _validate_percentile(name: str, value: float | None) -> float | None:
    """Validate an optional percentile argument and return it as ``float``."""

    if value is None:
        return None
    value = float(value)
    if not 0.0 <= value <= 100.0:
        raise ValueError(f"{name} must be between 0 and 100. Got {value!r}.")
    return value


def _validate_alpha_estimation_method(method: str) -> str:
    """Normalize and validate a two-channel alpha-estimation method name."""

    method = str(method).strip().lower()
    if method not in SUPPORTED_ALPHA_ESTIMATION_METHODS:
        raise ValueError(
            f"method must be one of {sorted(SUPPORTED_ALPHA_ESTIMATION_METHODS)}. "
            f"Got {method!r}."
        )
    return method


def _validate_positive_int(name: str, value: int, *, minimum: int = 1) -> int:
    """Validate that an integer parameter satisfies a lower bound."""

    value = int(value)
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum}. Got {value!r}.")
    return value


def _validate_positive_float(name: str, value: float, *, strictly_positive: bool = True) -> float:
    """Validate that a floating-point parameter satisfies a lower bound."""

    value = float(value)
    if strictly_positive and value <= 0.0:
        raise ValueError(f"{name} must be > 0. Got {value!r}.")
    if not strictly_positive and value < 0.0:
        raise ValueError(f"{name} must be >= 0. Got {value!r}.")
    return value


def _safe_correlation(x: np.ndarray, y: np.ndarray) -> float:
    """Return a finite Pearson correlation when both vectors contain variance."""

    if x.size == 0 or y.size == 0:
        return float("nan")
    if np.std(x) <= EPSILON or np.std(y) <= EPSILON:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def mutual_information_1d(x, y, bins: int = DEFAULT_MI_BINS) -> float:
    """
    Estimate mutual information between two 1D intensity arrays using a 2D histogram.

    Parameters
    ----------
    x, y : array-like
        One-dimensional or flattenable intensity arrays of equal length.
    bins : int, optional
        Number of histogram bins used for the joint density estimate. Must be
        at least ``2``.

    Returns
    -------
    float
        Estimated mutual information in natural-log units.

    Notes
    -----
    This is a histogram-based estimator intended for relative comparisons
    during alpha optimization, not for bias-corrected information-theoretic
    inference.
    """

    bins = _validate_positive_int("mi_bins", bins, minimum=2)
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()

    if x.shape != y.shape:
        raise ValueError(
            f"x and y must have the same flattened shape. Got {x.shape!r} and {y.shape!r}."
        )
    if x.size == 0:
        raise ValueError("x and y must not be empty.")

    hist_2d, _, _ = np.histogram2d(x, y, bins=bins)
    total = float(np.sum(hist_2d))
    if total <= 0.0:
        return 0.0

    pxy = hist_2d / total
    px = np.sum(pxy, axis=1)
    py = np.sum(pxy, axis=0)
    px_py = px[:, None] * py[None, :]
    nz = pxy > 0

    return float(np.sum(pxy[nz] * np.log(pxy[nz] / px_py[nz])))


def _prepare_single_volume_for_alpha(
    volume,
    *,
    background_percentile: float = 1.0,
    preprocess_alpha_inputs: bool = True,
) -> tuple[np.ndarray, float]:
    """Prepare one image volume for alpha estimation and return its background."""

    volume = np.asarray(volume, dtype=np.float32)
    background_percentile = _validate_percentile(
        "background_percentile",
        background_percentile,
    )

    if not preprocess_alpha_inputs:
        return volume.copy(), 0.0

    background = float(np.percentile(volume, background_percentile))
    prepared = np.clip(volume - background, a_min=0.0, a_max=None)
    return prepared.astype(np.float32, copy=False), background


def prepare_source_target_for_alpha(
    source_volume,
    target_volume,
    *,
    background_percentile: float = 1.0,
    preprocess_alpha_inputs: bool = True,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """
    Convert, optionally background-correct, and clip source and target volumes.

    Parameters
    ----------
    source_volume, target_volume : array-like
        Matching source and target image volumes. Any matching shape is
        accepted, for example ``ZYX`` or flattened arrays.
    background_percentile : float, optional
        Low percentile used to estimate a rough background in each input.
    preprocess_alpha_inputs : bool, optional
        If ``True``, subtract the percentile-based background from each input
        and clip negative values to zero. If ``False``, the inputs are only
        converted to ``float32``.

    Returns
    -------
    tuple
        ``(source_prepared, target_prepared, source_background, target_background)``.

    Raises
    ------
    ValueError
        If the input volumes do not share the same shape or are empty.
    """

    source_f = np.asarray(source_volume, dtype=np.float32)
    target_f = np.asarray(target_volume, dtype=np.float32)

    if source_f.shape != target_f.shape:
        raise ValueError(
            "source and target must have the same shape. "
            f"Got {source_f.shape!r} and {target_f.shape!r}."
        )
    if source_f.size == 0:
        raise ValueError("source and target must not be empty.")

    source_prepared, source_background = _prepare_single_volume_for_alpha(
        source_f,
        background_percentile=background_percentile,
        preprocess_alpha_inputs=preprocess_alpha_inputs,
    )
    target_prepared, target_background = _prepare_single_volume_for_alpha(
        target_f,
        background_percentile=background_percentile,
        preprocess_alpha_inputs=preprocess_alpha_inputs,
    )
    return source_prepared, target_prepared, source_background, target_background


def make_alpha_mask(
    source,
    target=None,
    *,
    signal_percentile: float = 99.0,
    target_low_percentile: float | None = None,
    min_voxels: int = MIN_MASK_VOXELS,
) -> tuple[np.ndarray, dict]:
    """
    Create a robust mask for alpha estimation.

    Parameters
    ----------
    source : array-like
        Prepared source intensities used to define the bright-source mask.
    target : array-like or None, optional
        Optional prepared target intensities used to additionally constrain the
        mask to low-target voxels.
    signal_percentile : float, optional
        Source percentile above which voxels are considered signal-rich.
    target_low_percentile : float or None, optional
        If provided, only voxels at or below this target percentile are kept,
        unless that mask would become too small and the source-only fallback is
        used instead.
    min_voxels : int, optional
        Minimum number of voxels required for a valid estimation mask.

    Returns
    -------
    tuple
        ``(mask, details)`` where ``mask`` is a boolean array and ``details`` is
        a metadata dictionary describing thresholds, fallback behavior, and voxel
        counts.
    """

    source = np.asarray(source, dtype=np.float32)
    target = None if target is None else np.asarray(target, dtype=np.float32)
    signal_percentile = _validate_percentile("signal_percentile", signal_percentile)
    target_low_percentile = _validate_percentile(
        "target_low_percentile",
        target_low_percentile,
    )
    min_voxels = _validate_positive_int("min_voxels", min_voxels, minimum=1)

    if target is not None and source.shape != target.shape:
        raise ValueError(
            "source and target must have the same shape for mask creation. "
            f"Got {source.shape!r} and {target.shape!r}."
        )

    source_threshold = float(np.percentile(source, signal_percentile))
    mask_source_only = source >= source_threshold
    source_only_voxels = int(np.count_nonzero(mask_source_only))

    mask = mask_source_only
    target_threshold = None
    fallback_used = False
    strategy = "source_only"

    if target is not None and target_low_percentile is not None:
        target_threshold = float(np.percentile(target, target_low_percentile))
        mask_candidate = mask_source_only & (target <= target_threshold)
        candidate_voxels = int(np.count_nonzero(mask_candidate))
        if candidate_voxels >= min_voxels:
            mask = mask_candidate
            strategy = "source_and_low_target"
        elif source_only_voxels >= min_voxels:
            fallback_used = True
        else:
            raise ValueError(
                "Alpha mask is too small even after fallback to the source-only mask. "
                f"Candidate voxels={candidate_voxels}, source-only voxels={source_only_voxels}, "
                f"required min_voxels={min_voxels}."
            )
    elif source_only_voxels < min_voxels:
        raise ValueError(
            "Source signal mask does not contain enough voxels for alpha estimation. "
            f"Found {source_only_voxels}, need at least {min_voxels}."
        )

    mask_voxels = int(np.count_nonzero(mask))
    if mask_voxels < min_voxels:
        raise ValueError(
            "Final alpha mask does not contain enough voxels for estimation. "
            f"Found {mask_voxels}, need at least {min_voxels}."
        )

    details = {
        "signal_percentile": float(signal_percentile),
        "target_low_percentile": None
        if target_low_percentile is None
        else float(target_low_percentile),
        "source_threshold": source_threshold,
        "target_threshold": target_threshold,
        "mask_strategy": strategy,
        "mask_fallback_used": bool(fallback_used),
        "mask_voxel_count": mask_voxels,
        "source_only_mask_voxel_count": source_only_voxels,
        "min_mask_voxels": int(min_voxels),
    }
    return mask, details


def _subsample_masked_vectors(
    x: np.ndarray,
    y: np.ndarray,
    *,
    max_alpha_voxels: int | None = DEFAULT_MAX_ALPHA_VOXELS,
    random_state: int = DEFAULT_RANDOM_STATE,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Optionally subsample masked source/target vectors for faster estimation."""

    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()

    if x.shape != y.shape:
        raise ValueError(
            f"x and y must have the same shape. Got {x.shape!r} and {y.shape!r}."
        )

    original_count = int(x.size)
    if max_alpha_voxels is None or original_count <= int(max_alpha_voxels):
        return x, y, {
            "subsampled": False,
            "voxel_count_before_subsampling": original_count,
            "voxel_count_after_subsampling": original_count,
            "max_alpha_voxels": None if max_alpha_voxels is None else int(max_alpha_voxels),
            "random_state": int(random_state),
        }

    max_alpha_voxels = _validate_positive_int(
        "max_alpha_voxels",
        max_alpha_voxels,
        minimum=1,
    )
    rng = np.random.default_rng(int(random_state))
    indices = rng.choice(original_count, size=max_alpha_voxels, replace=False)
    details = {
        "subsampled": True,
        "voxel_count_before_subsampling": original_count,
        "voxel_count_after_subsampling": int(max_alpha_voxels),
        "max_alpha_voxels": int(max_alpha_voxels),
        "random_state": int(random_state),
    }
    return x[indices], y[indices], details


def _estimate_alpha_mean_ratio(x: np.ndarray, y: np.ndarray) -> tuple[float, dict]:
    """Estimate alpha from the ratio of masked mean target and source intensities."""

    denominator = float(np.mean(x))
    if denominator <= EPSILON:
        raise ValueError("Mean source intensity inside the estimation mask must be > 0.")
    alpha = max(float(np.mean(y)) / denominator, 0.0)
    return alpha, {
        "numerator_mean": float(np.mean(y)),
        "denominator_mean": denominator,
    }


def _estimate_alpha_linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, dict]:
    """Estimate alpha by least-squares fitting of ``y ≈ alpha * x`` without intercept."""

    denominator = float(np.sum(x * x))
    if denominator <= EPSILON:
        raise ValueError(
            "Sum of squared source intensities inside the estimation mask must be > 0."
        )
    numerator = float(np.sum(x * y))
    alpha = max(numerator / denominator, 0.0)
    return alpha, {
        "numerator_sum_xy": numerator,
        "denominator_sum_x2": denominator,
    }


def _estimate_alpha_corr_min(
    x: np.ndarray,
    y: np.ndarray,
    *,
    alpha_max: float,
) -> tuple[float, dict]:
    """Estimate alpha by minimizing residual Pearson correlation after correction."""

    alpha_max = _validate_positive_float("alpha_max", alpha_max)

    def objective(alpha_value: float) -> float:
        corrected = y - float(alpha_value) * x
        r = _safe_correlation(x, corrected)
        if not np.isfinite(r):
            return float("inf")
        return float(r * r)

    result = minimize_scalar(objective, bounds=(0.0, alpha_max), method="bounded")
    if not result.success or not np.isfinite(result.x):
        raise ValueError("Correlation-minimization alpha estimation failed.")

    alpha = max(float(result.x), 0.0)
    corrected = y - alpha * x
    r = _safe_correlation(x, corrected)
    return alpha, {
        "optimization_success": bool(result.success),
        "optimization_nfev": int(getattr(result, "nfev", -1)),
        "objective_value": float(result.fun),
        "post_correlation": None if not np.isfinite(r) else float(r),
    }


def _estimate_alpha_mi_min(
    x: np.ndarray,
    y: np.ndarray,
    *,
    alpha_max: float,
    mi_bins: int,
) -> tuple[float, dict]:
    """Estimate alpha by minimizing residual mutual information after correction."""

    alpha_max = _validate_positive_float("alpha_max", alpha_max)
    mi_bins = _validate_positive_int("mi_bins", mi_bins, minimum=2)

    def objective(alpha_value: float) -> float:
        corrected = y - float(alpha_value) * x
        mi = mutual_information_1d(x, corrected, bins=mi_bins)
        if not np.isfinite(mi):
            return float("inf")
        return float(mi)

    result = minimize_scalar(objective, bounds=(0.0, alpha_max), method="bounded")
    if not result.success or not np.isfinite(result.x):
        raise ValueError("Mutual-information-minimization alpha estimation failed.")

    alpha = max(float(result.x), 0.0)
    corrected = y - alpha * x
    mi = mutual_information_1d(x, corrected, bins=mi_bins)
    return alpha, {
        "optimization_success": bool(result.success),
        "optimization_nfev": int(getattr(result, "nfev", -1)),
        "objective_value": float(result.fun),
        "post_mutual_information": float(mi),
    }


def estimate_alpha_from_volume(
    source,
    target,
    signal_percentile: float = 99.0,
    background_percentile: float = 1.0,
    min_mask_voxels: int = MIN_MASK_VOXELS,
    *,
    method: str = "mean_ratio",
    target_low_percentile: float | None = None,
    preprocess_alpha_inputs: bool = True,
    alpha_max: float = DEFAULT_ALPHA_MAX,
    mi_bins: int = DEFAULT_MI_BINS,
    max_alpha_voxels: int | None = DEFAULT_MAX_ALPHA_VOXELS,
    random_state: int = DEFAULT_RANDOM_STATE,
    return_details: bool = False,
) -> float | tuple[float, dict]:
    """
    Estimate a bleed-through coefficient alpha from matching source and target volumes.

    Parameters
    ----------
    source, target : array-like
        Matching source and target image volumes. Typical microscopy input is a
        ``ZYX`` volume, but any matching shape is accepted.
    signal_percentile : float, optional
        Percentile used to define a bright-source signal mask.
    background_percentile : float, optional
        Low percentile used for optional percentile-based background subtraction
        before alpha estimation.
    min_mask_voxels : int, optional
        Minimum number of voxels required in the final alpha-estimation mask.
    method : {"mean_ratio", "linear_fit", "corr_min", "mi_min"}, optional
        Scalar alpha-estimation strategy.
    target_low_percentile : float or None, optional
        Optional extra target constraint for the estimation mask. If provided,
        voxels are restricted to bright source signal and low target signal
        whenever that yields enough voxels.
    preprocess_alpha_inputs : bool, optional
        If ``True``, apply percentile-based background subtraction and clip
        negative values to zero before alpha estimation. If ``False``, estimate
        alpha directly on the original intensities converted to ``float32``.
    alpha_max : float, optional
        Upper search bound for optimization-based methods ``"corr_min"`` and
        ``"mi_min"``.
    mi_bins : int, optional
        Number of histogram bins used by the mutual-information estimator.
    max_alpha_voxels : int or None, optional
        Optional cap on the number of voxels used after masking. If the mask is
        larger, voxels are subsampled without replacement.
    random_state : int, optional
        Random seed used for optional voxel subsampling.
    return_details : bool, optional
        If ``True``, also return a dictionary describing preprocessing,
        thresholds, mask size, subsampling, and method-specific diagnostics.

    Returns
    -------
    float or tuple
        Estimated alpha, or ``(alpha, details)`` if ``return_details=True``.

    Raises
    ------
    ValueError
        If the input data are incompatible, the estimation mask is too small,
        or the estimated alpha is invalid.
    """

    method = _validate_alpha_estimation_method(method)
    signal_percentile = _validate_percentile("signal_percentile", signal_percentile)
    background_percentile = _validate_percentile(
        "background_percentile",
        background_percentile,
    )
    target_low_percentile = _validate_percentile(
        "target_low_percentile",
        target_low_percentile,
    )
    min_mask_voxels = _validate_positive_int(
        "min_mask_voxels",
        min_mask_voxels,
        minimum=1,
    )
    alpha_max = _validate_positive_float("alpha_max", alpha_max)
    mi_bins = _validate_positive_int("mi_bins", mi_bins, minimum=2)

    source_prepared, target_prepared, source_background, target_background = (
        prepare_source_target_for_alpha(
            source,
            target,
            background_percentile=background_percentile,
            preprocess_alpha_inputs=preprocess_alpha_inputs,
        )
    )

    mask, mask_details = make_alpha_mask(
        source_prepared,
        target_prepared,
        signal_percentile=signal_percentile,
        target_low_percentile=target_low_percentile,
        min_voxels=min_mask_voxels,
    )

    x = np.asarray(source_prepared[mask], dtype=np.float64)
    y = np.asarray(target_prepared[mask], dtype=np.float64)
    x, y, subsampling_details = _subsample_masked_vectors(
        x,
        y,
        max_alpha_voxels=max_alpha_voxels,
        random_state=random_state,
    )

    if x.size < 2:
        raise ValueError(
            "Not enough masked voxels remain after optional subsampling for alpha estimation."
        )

    if method == "mean_ratio":
        alpha, method_details = _estimate_alpha_mean_ratio(x, y)
    elif method == "linear_fit":
        alpha, method_details = _estimate_alpha_linear_fit(x, y)
    elif method == "corr_min":
        alpha, method_details = _estimate_alpha_corr_min(
            x,
            y,
            alpha_max=alpha_max,
        )
    else:
        alpha, method_details = _estimate_alpha_mi_min(
            x,
            y,
            alpha_max=alpha_max,
            mi_bins=mi_bins,
        )

    if not np.isfinite(alpha):
        raise ValueError(f"Estimated alpha is not finite: {alpha!r}.")

    details = {
        "method": method,
        "alpha": float(alpha),
        "signal_percentile": float(signal_percentile),
        "target_low_percentile": None
        if target_low_percentile is None
        else float(target_low_percentile),
        "background_percentile": float(background_percentile),
        "preprocess_alpha_inputs": bool(preprocess_alpha_inputs),
        "source_background": float(source_background),
        "target_background": float(target_background),
        "alpha_max": float(alpha_max),
        "mi_bins": int(mi_bins),
        "mask_voxel_count": int(mask_details["mask_voxel_count"]),
        "mask_details": mask_details,
        "subsampling_details": subsampling_details,
        "estimation_voxel_count": int(x.size),
        "method_details": method_details,
    }

    if return_details:
        return float(alpha), details
    return float(alpha)


def estimate_picasso_unmixing_matrix_from_volume(
    channel_volumes,
    *,
    background_percentile: float = 1.0,
    preprocess_alpha_inputs: bool = True,
    mi_bins: int = DEFAULT_MI_BINS,
    alpha_max: float = DEFAULT_ALPHA_MAX,
    max_iter: int = 10,
    tolerance: float = 1e-4,
    max_alpha_voxels: int | None = DEFAULT_MAX_ALPHA_VOXELS,
    random_state: int = DEFAULT_RANDOM_STATE,
) -> tuple[np.ndarray, dict]:
    """
    Estimate a PICASSO-like blind unmixing matrix from multi-channel image data.

    Parameters
    ----------
    channel_volumes : array-like
        Multi-channel image data with channel as the first axis, for example
        ``(C, Z, Y, X)`` or ``(C, N)``.
    background_percentile : float, optional
        Low percentile used for optional per-channel background subtraction.
    preprocess_alpha_inputs : bool, optional
        If ``True``, apply percentile-based background subtraction and clipping
        before estimating the unmixing matrix.
    mi_bins : int, optional
        Number of histogram bins used by the mutual-information estimator.
    alpha_max : float, optional
        Upper bound for each pairwise subtraction coefficient.
    max_iter : int, optional
        Maximum number of pairwise update sweeps.
    tolerance : float, optional
        Convergence criterion on the largest pairwise coefficient update in one
        iteration.
    max_alpha_voxels : int or None, optional
        Optional cap on the number of voxels used for matrix estimation.
    random_state : int, optional
        Random seed used for optional subsampling.

    Returns
    -------
    tuple
        ``(matrix, details)`` with the estimated unmixing matrix and a metadata
        dictionary describing convergence and pairwise updates.

    Notes
    -----
    This implements an iterative, pairwise, mutual-information-minimizing blind
    unmixing routine. It is inspired by the PICASSO criterion, but it is not a
    deep-learning approach and not a full spectral reference-based method.
    """

    channels_f = np.asarray(channel_volumes, dtype=np.float32)
    if channels_f.ndim < 2:
        raise ValueError(
            "channel_volumes must have channel as the first axis and at least one "
            f"additional data axis. Got shape {channels_f.shape!r}."
        )

    n_channels = int(channels_f.shape[0])
    if n_channels < 2:
        raise ValueError("PICASSO-style blind unmixing requires at least two channels.")

    alpha_max = _validate_positive_float("alpha_max", alpha_max)
    mi_bins = _validate_positive_int("mi_bins", mi_bins, minimum=2)
    max_iter = _validate_positive_int("max_iter", max_iter, minimum=1)
    tolerance = _validate_positive_float("tolerance", tolerance)

    prepared_channels = np.empty_like(channels_f, dtype=np.float32)
    channel_backgrounds: list[float] = []
    for c in range(n_channels):
        prepared_channels[c], background = _prepare_single_volume_for_alpha(
            channels_f[c],
            background_percentile=background_percentile,
            preprocess_alpha_inputs=preprocess_alpha_inputs,
        )
        channel_backgrounds.append(float(background))

    flattened = prepared_channels.reshape(n_channels, -1).astype(np.float64, copy=False)
    if max_alpha_voxels is not None and flattened.shape[1] > int(max_alpha_voxels):
        max_alpha_voxels = _validate_positive_int(
            "max_alpha_voxels",
            max_alpha_voxels,
            minimum=1,
        )
        rng = np.random.default_rng(int(random_state))
        indices = rng.choice(flattened.shape[1], size=max_alpha_voxels, replace=False)
        flattened = flattened[:, indices]
        estimation_voxel_count = int(max_alpha_voxels)
        subsampled = True
    else:
        estimation_voxel_count = int(flattened.shape[1])
        subsampled = False

    F = flattened.copy()
    U = np.eye(n_channels, dtype=np.float64)
    pairwise_updates: list[dict] = []
    converged = False
    iterations_run = 0

    for iteration in range(max_iter):
        max_change = 0.0
        for i in range(n_channels):
            for j in range(n_channels):
                if i == j:
                    continue
                x = F[i]
                y = F[j]
                if np.std(x) <= EPSILON or np.std(y) <= EPSILON:
                    alpha = 0.0
                    mi_before = mutual_information_1d(x, y, bins=mi_bins)
                    mi_after = mi_before
                else:
                    alpha, method_details = _estimate_alpha_mi_min(
                        x,
                        y,
                        alpha_max=alpha_max,
                        mi_bins=mi_bins,
                    )
                    mi_before = mutual_information_1d(x, y, bins=mi_bins)
                    y_new = y - alpha * x
                    mi_after = mutual_information_1d(x, y_new, bins=mi_bins)
                    F[j] = y_new
                    U[j, :] = U[j, :] - alpha * U[i, :]
                    max_change = max(max_change, float(alpha))
                    pairwise_updates.append(
                        {
                            "iteration": int(iteration),
                            "source_index": int(i),
                            "target_index": int(j),
                            "alpha": float(alpha),
                            "mi_before": float(mi_before),
                            "mi_after": float(mi_after),
                            "optimization_nfev": int(
                                method_details.get("optimization_nfev", -1)
                            ),
                        }
                    )
        iterations_run = iteration + 1
        if max_change < tolerance:
            converged = True
            break

    details = {
        "background_percentile": float(background_percentile),
        "preprocess_alpha_inputs": bool(preprocess_alpha_inputs),
        "channel_backgrounds": channel_backgrounds,
        "alpha_max": float(alpha_max),
        "mi_bins": int(mi_bins),
        "max_iter": int(max_iter),
        "tolerance": float(tolerance),
        "iterations_run": int(iterations_run),
        "converged": bool(converged),
        "estimation_voxel_count": int(estimation_voxel_count),
        "subsampled": bool(subsampled),
        "max_alpha_voxels": None if max_alpha_voxels is None else int(max_alpha_voxels),
        "random_state": int(random_state),
        "pairwise_updates": pairwise_updates,
        "unmixing_matrix": U.tolist(),
    }
    return U.astype(np.float64, copy=False), details


__all__ = [
    "DEFAULT_ALPHA_MAX",
    "DEFAULT_MAX_ALPHA_VOXELS",
    "DEFAULT_MI_BINS",
    "DEFAULT_RANDOM_STATE",
    "MIN_MASK_VOXELS",
    "SUPPORTED_ALPHA_ESTIMATION_METHODS",
    "estimate_alpha_from_volume",
    "estimate_picasso_unmixing_matrix_from_volume",
    "make_alpha_mask",
    "mutual_information_1d",
    "prepare_source_target_for_alpha",
]
