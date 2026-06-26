"""
Internal PICASSO-style unmixing implementations.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy.optimize import minimize

from .estimation import (
    DEFAULT_ALPHA_MAX,
    DEFAULT_MAX_ALPHA_VOXELS,
    DEFAULT_MI_BINS,
    DEFAULT_RANDOM_STATE,
    mutual_information_1d,
)
# %% CONSTANTS
SUPPORTED_PICASSO_IMPLEMENTATIONS = {
    "matlab_3c",
    "matlab_n",
    "source_sink_n",
}
DEFAULT_PICASSO_IMPLEMENTATION = "matlab_3c"
DEFAULT_PICASSO_STEP_SIZE = 0.2
DEFAULT_PICASSO_QN = 100
DEFAULT_PICASSO_BIN_FACTOR = 16
DEFAULT_PICASSO_ALPHA_CLIP = 0.5
DEFAULT_PICASSO_NEGATIVITY_THRESHOLD = 0.9e-3
DEFAULT_PICASSO_CLIP_EVERY_N_ITERATIONS = 50

EPSILON = 1e-12

# %% INTERNAL HELPERS
def validate_picasso_implementation(implementation: str) -> str:
    """Normalize and validate the selected PICASSO implementation name."""

    implementation = str(implementation).strip().lower()
    if implementation not in SUPPORTED_PICASSO_IMPLEMENTATIONS:
        raise ValueError(
            "implementation must be one of "
            f"{sorted(SUPPORTED_PICASSO_IMPLEMENTATIONS)}. Got {implementation!r}."
        )
    return implementation


def validate_source_sink_matrix(source_sink_matrix, *, n_channels: int) -> np.ndarray:
    """
    Validate a square source-sink mixing matrix.

    The matrix uses the napari-plugin-inspired convention:

    - diagonal entry ``1``: the corresponding channel is the sink to be corrected
    - off-diagonal entry ``-1``: the corresponding row channel contributes to that sink
    - off-diagonal entry ``0``: no modeled spillover relation
    """

    matrix = np.asarray(source_sink_matrix, dtype=np.int8)
    if matrix.shape != (n_channels, n_channels):
        raise ValueError(
            "source_sink_matrix must have shape "
            f"({n_channels}, {n_channels}). Got {matrix.shape!r}."
        )
    if not np.all(np.diag(matrix) == 1):
        raise ValueError("source_sink_matrix must have ones on the diagonal.")
    off_diagonal_mask = ~np.eye(n_channels, dtype=bool)
    off_diagonal = matrix[off_diagonal_mask]
    if not np.all(np.isin(off_diagonal, (-1, 0))):
        raise ValueError(
            "Off-diagonal entries of source_sink_matrix must be -1 or 0."
        )
    return matrix


def default_source_sink_matrix(n_channels: int) -> np.ndarray:
    """Return an all-to-all square source-sink matrix with identity sinks."""

    matrix = -np.ones((n_channels, n_channels), dtype=np.int8)
    np.fill_diagonal(matrix, 1)
    return matrix


def build_source_sink_matrix(
    selected_channels: Sequence[int],
    *,
    sink_channels: Sequence[int] | None = None,
    neutral_channels: Sequence[int] | None = None,
) -> np.ndarray:
    """
    Build a source-sink matrix from channel-role lists.

    Parameters
    ----------
    selected_channels : sequence of int
        The actual channel indices included in the current unmixing run.
    sink_channels : sequence of int or None, optional
        Actual channel indices, drawn from ``selected_channels``, that should be
        corrected as sinks. Every non-neutral selected channel may contribute to
        these sinks except the sink channel itself.
    neutral_channels : sequence of int or None, optional
        Actual channel indices, drawn from ``selected_channels``, that should be
        excluded from source-sink modeling. Neutral channels are kept in the
        output but are neither corrected as sinks nor used as sources.

    Returns
    -------
    numpy.ndarray
        Square ``int8`` source-sink matrix following the package convention:
        diagonal ``1``, modeled source-to-sink relations ``-1``, and all other
        entries ``0``.

    Notes
    -----
    This helper is intended as a more readable alternative to writing
    ``source_sink_matrix`` manually. It operates on the user-facing channel
    indices rather than on channel positions inside the selected subset.
    """

    selected_channels = [int(channel) for channel in selected_channels]
    if len(selected_channels) < 2:
        raise ValueError(
            "selected_channels must contain at least two channels to build a "
            "source-sink matrix."
        )

    selected_set = set(selected_channels)
    neutral_channels = [] if neutral_channels is None else [int(channel) for channel in neutral_channels]
    neutral_set = set(neutral_channels)
    if not neutral_set.issubset(selected_set):
        missing = sorted(neutral_set.difference(selected_set))
        raise ValueError(
            "neutral_channels must be a subset of selected_channels. "
            f"Missing from selected_channels: {missing!r}."
        )

    if sink_channels is None:
        sink_channels = [
            channel for channel in selected_channels if channel not in neutral_set
        ]
    else:
        sink_channels = [int(channel) for channel in sink_channels]

    sink_set = set(sink_channels)
    if not sink_set:
        raise ValueError("sink_channels must contain at least one non-neutral channel.")
    if not sink_set.issubset(selected_set):
        missing = sorted(sink_set.difference(selected_set))
        raise ValueError(
            "sink_channels must be a subset of selected_channels. "
            f"Missing from selected_channels: {missing!r}."
        )
    if sink_set.intersection(neutral_set):
        overlap = sorted(sink_set.intersection(neutral_set))
        raise ValueError(
            "sink_channels and neutral_channels must be disjoint. "
            f"Overlapping channels: {overlap!r}."
        )

    channel_to_position = {channel: index for index, channel in enumerate(selected_channels)}
    matrix = np.zeros((len(selected_channels), len(selected_channels)), dtype=np.int8)
    np.fill_diagonal(matrix, 1)

    source_channels = [
        channel for channel in selected_channels if channel not in neutral_set
    ]
    for source_channel in source_channels:
        source_position = channel_to_position[source_channel]
        for sink_channel in sink_channels:
            if source_channel == sink_channel:
                continue
            sink_position = channel_to_position[sink_channel]
            matrix[source_position, sink_position] = -1

    return matrix


def _normalize_positive_int(name: str, value: int, *, minimum: int = 1) -> int:
    value = int(value)
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum}. Got {value!r}.")
    return value


def _normalize_positive_float(name: str, value: float, *, minimum: float = 0.0) -> float:
    value = float(value)
    if value <= minimum:
        relation = ">" if minimum == 0.0 else f"> {minimum}"
        raise ValueError(f"{name} must be {relation}. Got {value!r}.")
    return value


def _prepare_matlab_channels(
    channel_volumes: np.ndarray,
    *,
    background_percentile: float,
) -> tuple[np.ndarray, list[float], list[float]]:
    """Apply the MATLAB PICASSO preprocessing: per-channel max normalization and background removal."""

    prepared = np.empty_like(channel_volumes, dtype=np.float32)
    normalization_factors: list[float] = []
    background_values: list[float] = []

    for c in range(channel_volumes.shape[0]):
        current = np.asarray(channel_volumes[c], dtype=np.float32)
        norm_factor = float(np.max(current))
        if norm_factor > 0.0:
            current = current / norm_factor
        background = float(np.percentile(current, background_percentile))
        current = np.maximum(current - background, 0.0)
        prepared[c] = current.astype(np.float32, copy=False)
        normalization_factors.append(norm_factor)
        background_values.append(background)

    return prepared, normalization_factors, background_values


def _prepare_source_sink_channels(
    channel_volumes: np.ndarray,
    *,
    background_percentile: float,
) -> tuple[np.ndarray, list[float], float]:
    """Normalize all selected channels by a shared maximum and estimate one background per channel."""

    channel_volumes = np.asarray(channel_volumes, dtype=np.float32)
    global_max = float(np.max(channel_volumes))
    if global_max > 0.0:
        normalized = channel_volumes / global_max
    else:
        normalized = channel_volumes.copy()

    backgrounds: list[float] = []
    for c in range(normalized.shape[0]):
        backgrounds.append(float(np.percentile(normalized[c], background_percentile)))

    return normalized.astype(np.float32, copy=False), backgrounds, global_max


def _matlab_check_negativity(values: np.ndarray) -> float:
    """Compute the MATLAB-style negativity ratio used by PICASSO."""

    values = np.asarray(values, dtype=np.float64).ravel()
    total_abs = float(np.sum(np.abs(values)))
    if total_abs <= EPSILON:
        return 0.0
    negative_abs = float(np.sum(np.abs(values[values < 0])))
    return negative_abs / total_abs


def _matlab_pixel_bin_plane(plane: np.ndarray, bin_factor: int) -> np.ndarray:
    """Replicate the MATLAB ``pixelBin`` helper for one 2D plane."""

    if plane.ndim != 2:
        raise ValueError(f"Expected a 2D plane. Got shape {plane.shape!r}.")

    bin_factor = _normalize_positive_int("pixel_bin_size", bin_factor, minimum=1)
    height = (plane.shape[0] // bin_factor) * bin_factor
    width = (plane.shape[1] // bin_factor) * bin_factor
    cropped = np.asarray(plane[:height, :width], dtype=np.float64)
    if cropped.size == 0:
        return np.zeros((0, 0), dtype=np.float64)

    reshaped = cropped.reshape(height // bin_factor, bin_factor, width // bin_factor, bin_factor)
    return reshaped.sum(axis=(1, 3))


def _matlab_bin_volume(volume: np.ndarray, bin_factor: int) -> np.ndarray:
    """Apply MATLAB-style 2D pixel binning slice-wise and flatten the result."""

    volume = np.asarray(volume, dtype=np.float64)
    if volume.ndim == 2:
        return _matlab_pixel_bin_plane(volume, bin_factor).ravel()
    if volume.ndim < 2:
        raise ValueError(f"Expected at least 2D input. Got shape {volume.shape!r}.")

    flattened_planes: list[np.ndarray] = []
    for index in np.ndindex(volume.shape[:-2]):
        flattened_planes.append(_matlab_pixel_bin_plane(volume[index], bin_factor).ravel())
    if not flattened_planes:
        return np.zeros((0,), dtype=np.float64)
    return np.concatenate(flattened_planes, axis=0)


def _matlab_mutual_information(im1: np.ndarray, im2: np.ndarray, qn: int) -> float:
    """Replicate the MATLAB histogram-based mutual information calculation."""

    qn = _normalize_positive_int("qN", qn, minimum=1)
    im1 = np.asarray(im1, dtype=np.float64).ravel()
    im2 = np.asarray(im2, dtype=np.float64).ravel()
    if im1.shape != im2.shape:
        raise ValueError(
            f"im1 and im2 must have the same shape. Got {im1.shape!r} and {im2.shape!r}."
        )
    if im1.size == 0:
        return 0.0

    scale1 = float(np.max(np.abs(im1)))
    scale2 = float(np.max(np.abs(im2)))
    if scale1 <= EPSILON:
        q1 = np.zeros_like(im1, dtype=np.int8)
    else:
        q1 = (qn * im1 / scale1).astype(np.int8, copy=False)
    if scale2 <= EPSILON:
        q2 = np.zeros_like(im2, dtype=np.int8)
    else:
        q2 = (qn * im2 / scale2).astype(np.int8, copy=False)

    _, inv1 = np.unique(q1, return_inverse=True)
    _, inv2 = np.unique(q2, return_inverse=True)
    joint_hist = np.zeros((int(inv1.max()) + 1, int(inv2.max()) + 1), dtype=np.int64)
    np.add.at(joint_hist, (inv1, inv2), 1)

    joint_prob = joint_hist.astype(np.float64) / float(im1.size)
    joint_nz = joint_prob[joint_prob > 0]
    joint_entropy = -float(np.sum(joint_nz * np.log2(joint_nz)))

    hist1 = np.sum(joint_hist, axis=1).astype(np.float64)
    hist2 = np.sum(joint_hist, axis=0).astype(np.float64)
    prob1 = hist1[hist1 > 0]
    prob2 = hist2[hist2 > 0]
    prob1 = prob1 / np.sum(prob1)
    prob2 = prob2 / np.sum(prob2)
    entropy1 = -float(np.sum(prob1 * np.log2(prob1)))
    entropy2 = -float(np.sum(prob2 * np.log2(prob2)))
    return entropy1 + entropy2 - joint_entropy


def _estimate_matlab_pair_alpha(
    sink_volume: np.ndarray,
    source_volume: np.ndarray,
    *,
    pixel_bin_size: int,
    qn: int,
    alpha_clip: float,
    negativity_threshold: float,
) -> tuple[float, dict]:
    """Estimate one MATLAB-style pairwise PICASSO coefficient."""

    sink_binned = _matlab_bin_volume(sink_volume, pixel_bin_size)
    source_binned = _matlab_bin_volume(source_volume, pixel_bin_size)
    if sink_binned.shape != source_binned.shape:
        raise ValueError(
            "Binned sink and source vectors must match. "
            f"Got {sink_binned.shape!r} and {source_binned.shape!r}."
        )
    if sink_binned.size == 0:
        return 0.0, {"alpha0": 0.0, "negativity_ratio": 0.0}

    norm = np.sqrt(np.dot(sink_binned, sink_binned) * np.dot(source_binned, source_binned))
    alpha0 = float(np.dot(sink_binned, source_binned) / norm) if norm > EPSILON else 0.0

    def objective(alpha_array: np.ndarray) -> float:
        alpha_value = float(alpha_array[0])
        return _matlab_mutual_information(
            sink_binned - alpha_value * source_binned,
            source_binned,
            qn,
        )

    result = minimize(
        objective,
        np.array([0.1 * alpha0], dtype=np.float64),
        method="Nelder-Mead",
        options={"maxiter": 200, "xatol": 1e-4, "fatol": 1e-4},
    )
    alpha = float(result.x[0]) if np.isfinite(result.x[0]) else 0.0
    alpha = float(np.clip(alpha, -alpha_clip, alpha_clip))
    negativity_ratio = _matlab_check_negativity(sink_binned - alpha * source_binned)
    if negativity_ratio > negativity_threshold:
        alpha = alpha / 10.0

    return alpha, {
        "alpha0": float(alpha0),
        "optimization_success": bool(result.success),
        "optimization_nit": int(getattr(result, "nit", -1)),
        "optimization_nfev": int(getattr(result, "nfev", -1)),
        "negativity_ratio": float(negativity_ratio),
        "pixel_bin_size": int(pixel_bin_size),
        "qN": int(qn),
    }


def _apply_sequential_update(channel_volumes: np.ndarray, incremental_matrix: np.ndarray) -> np.ndarray:
    """Apply one MATLAB-style in-place row update sweep to all channels."""

    current = np.asarray(channel_volumes, dtype=np.float32)
    updated = current.copy()
    n_channels = current.shape[0]
    for row in range(n_channels):
        new_row = np.zeros_like(current[row], dtype=np.float32)
        for col in range(n_channels):
            source = updated[col] if col < row else current[col]
            new_row = new_row + float(incremental_matrix[row, col]) * source
        updated[row] = new_row
    return updated


def run_picasso_matlab_like(
    channel_volumes,
    *,
    background_percentile: float,
    max_iter: int,
    step_size: float,
    qn: int,
    pixel_bin_size: int,
    alpha_clip: float,
    negativity_threshold: float,
    clip_every_n_iterations: int,
    require_three_channels: bool,
) -> tuple[np.ndarray, dict]:
    """Run the MATLAB-like PICASSO iteration on one multi-channel volume."""

    channel_volumes = np.asarray(channel_volumes, dtype=np.float32)
    if channel_volumes.ndim < 3:
        raise ValueError(
            "channel_volumes must have channel first and at least 2 spatial dimensions. "
            f"Got shape {channel_volumes.shape!r}."
        )
    n_channels = int(channel_volumes.shape[0])
    if require_three_channels and n_channels != 3:
        raise ValueError(
            f"The matlab_3c implementation requires exactly 3 channels. Got {n_channels}."
        )

    max_iter = _normalize_positive_int("max_iter", max_iter, minimum=1)
    step_size = _normalize_positive_float("step_size", step_size)
    qn = _normalize_positive_int("qN", qn, minimum=1)
    pixel_bin_size = _normalize_positive_int("pixel_bin_size", pixel_bin_size, minimum=1)
    alpha_clip = _normalize_positive_float("alpha_clip", alpha_clip)
    clip_every_n_iterations = _normalize_positive_int(
        "clip_every_n_iterations",
        clip_every_n_iterations,
        minimum=1,
    )

    prepared, normalization_factors, background_values = _prepare_matlab_channels(
        channel_volumes,
        background_percentile=background_percentile,
    )
    current = prepared.copy()
    cumulative_matrix = np.eye(n_channels, dtype=np.float64)
    incremental_matrices: list[list[list[float]]] = []
    raw_alpha_matrices: list[list[list[float]]] = []
    pairwise_alpha_details: list[list[dict]] = []
    max_abs_alpha_by_iteration: list[float] = []
    negativity_by_iteration: list[list[float]] = []

    for iteration in range(max_iter):
        alpha_matrix = np.zeros((n_channels, n_channels), dtype=np.float64)
        iteration_details: list[dict] = []
        for sink in range(n_channels):
            for source in range(n_channels):
                if sink == source:
                    continue
                alpha, details = _estimate_matlab_pair_alpha(
                    current[sink],
                    current[source],
                    pixel_bin_size=pixel_bin_size,
                    qn=qn,
                    alpha_clip=alpha_clip,
                    negativity_threshold=negativity_threshold,
                )
                alpha_matrix[sink, source] = -float(alpha)
                iteration_details.append(
                    {
                        "sink_index": int(sink),
                        "source_index": int(source),
                        "alpha": float(alpha),
                        **details,
                    }
                )

        incremental_matrix = np.eye(n_channels, dtype=np.float64) + step_size * alpha_matrix
        cumulative_matrix = incremental_matrix @ cumulative_matrix
        current = _apply_sequential_update(current, incremental_matrix)
        if (iteration + 1) % clip_every_n_iterations == 1:
            current = np.maximum(current, 0.0)

        incremental_matrices.append(incremental_matrix.tolist())
        raw_alpha_matrices.append(alpha_matrix.tolist())
        pairwise_alpha_details.append(iteration_details)
        max_abs_alpha_by_iteration.append(float(np.max(np.abs(alpha_matrix))))
        negativity_by_iteration.append(
            [float(_matlab_check_negativity(current[c])) for c in range(n_channels)]
        )

    final_volumes = np.maximum(current, 0.0).astype(np.float32, copy=False)
    details = {
        "background_percentile": float(background_percentile),
        "normalization_factors": [float(value) for value in normalization_factors],
        "background_values": [float(value) for value in background_values],
        "max_iter": int(max_iter),
        "step_size": float(step_size),
        "qN": int(qn),
        "pixel_bin_size": int(pixel_bin_size),
        "alpha_clip": float(alpha_clip),
        "negativity_threshold": float(negativity_threshold),
        "clip_every_n_iterations": int(clip_every_n_iterations),
        "iterations_run": int(max_iter),
        "incremental_matrices": incremental_matrices,
        "raw_alpha_matrices": raw_alpha_matrices,
        "pairwise_alpha_details": pairwise_alpha_details,
        "max_abs_alpha_by_iteration": max_abs_alpha_by_iteration,
        "negativity_by_iteration": negativity_by_iteration,
        "unmixing_matrix": cumulative_matrix.tolist(),
    }
    return final_volumes, details


def apply_matlab_incremental_sequence(
    channel_volumes,
    *,
    background_percentile: float,
    incremental_matrices: Sequence[Sequence[Sequence[float]]],
    clip_every_n_iterations: int,
) -> tuple[np.ndarray, dict]:
    """Apply a previously learned MATLAB-like update sequence to a new volume."""

    prepared, normalization_factors, background_values = _prepare_matlab_channels(
        np.asarray(channel_volumes, dtype=np.float32),
        background_percentile=background_percentile,
    )
    current = prepared.copy()
    clip_every_n_iterations = _normalize_positive_int(
        "clip_every_n_iterations",
        clip_every_n_iterations,
        minimum=1,
    )

    for iteration, matrix in enumerate(incremental_matrices):
        current = _apply_sequential_update(current, np.asarray(matrix, dtype=np.float64))
        if (iteration + 1) % clip_every_n_iterations == 1:
            current = np.maximum(current, 0.0)

    return np.maximum(current, 0.0).astype(np.float32, copy=False), {
        "background_values": [float(value) for value in background_values],
        "normalization_factors": [float(value) for value in normalization_factors],
        "iterations_run": int(len(incremental_matrices)),
    }


def _subsample_pair(
    source_vector: np.ndarray,
    sink_vector: np.ndarray,
    *,
    max_alpha_voxels: int | None,
    random_state: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Optionally subsample one source-sink vector pair for optimization."""

    source_vector = np.asarray(source_vector, dtype=np.float64).ravel()
    sink_vector = np.asarray(sink_vector, dtype=np.float64).ravel()
    if source_vector.shape != sink_vector.shape:
        raise ValueError(
            "source_vector and sink_vector must have the same shape. "
            f"Got {source_vector.shape!r} and {sink_vector.shape!r}."
        )
    if max_alpha_voxels is None or source_vector.size <= int(max_alpha_voxels):
        return source_vector, sink_vector

    max_alpha_voxels = _normalize_positive_int("max_alpha_voxels", max_alpha_voxels, minimum=1)
    rng = np.random.default_rng(int(random_state))
    indices = rng.choice(source_vector.size, size=max_alpha_voxels, replace=False)
    return source_vector[indices], sink_vector[indices]


def _estimate_source_sink_alpha(
    source_prepared: np.ndarray,
    sink_current: np.ndarray,
    *,
    mi_bins: int,
    alpha_max: float,
    max_alpha_voxels: int | None,
    random_state: int,
) -> tuple[float, dict]:
    """Estimate one plugin-like source-to-sink coefficient by mutual-information minimization."""

    source_vector, sink_vector = _subsample_pair(
        source_prepared,
        sink_current,
        max_alpha_voxels=max_alpha_voxels,
        random_state=random_state,
    )
    if source_vector.size == 0 or np.max(source_vector) <= EPSILON:
        return 0.0, {"optimization_success": True, "optimization_nfev": 0}

    def objective(alpha_value: np.ndarray) -> float:
        corrected = sink_vector - float(alpha_value[0]) * source_vector
        return mutual_information_1d(source_vector, corrected, bins=mi_bins)

    result = minimize(
        objective,
        np.array([0.1], dtype=np.float64),
        method="Nelder-Mead",
        options={"maxiter": 200, "xatol": 1e-4, "fatol": 1e-4},
    )
    alpha = float(result.x[0]) if np.isfinite(result.x[0]) else 0.0
    alpha = float(np.clip(alpha, 0.0, alpha_max))
    return alpha, {
        "optimization_success": bool(result.success),
        "optimization_nit": int(getattr(result, "nit", -1)),
        "optimization_nfev": int(getattr(result, "nfev", -1)),
    }


def run_source_sink_unmixing(
    channel_volumes,
    *,
    source_sink_matrix: np.ndarray,
    background_percentile: float,
    mi_bins: int,
    alpha_max: float,
    max_alpha_voxels: int | None,
    random_state: int,
) -> tuple[np.ndarray, dict]:
    """Run a napari-plugin-inspired source-sink multi-channel unmixing step."""

    normalized, backgrounds, global_max = _prepare_source_sink_channels(
        np.asarray(channel_volumes, dtype=np.float32),
        background_percentile=background_percentile,
    )
    n_channels = int(normalized.shape[0])
    source_sink_matrix = validate_source_sink_matrix(source_sink_matrix, n_channels=n_channels)
    mi_bins = _normalize_positive_int("mi_bins", mi_bins, minimum=2)
    alpha_max = _normalize_positive_float("alpha_max", alpha_max)

    corrected = normalized.copy()
    alpha_parameters = np.zeros((n_channels, n_channels), dtype=np.float64)
    background_parameters = np.zeros((n_channels, n_channels), dtype=np.float64)
    np.fill_diagonal(alpha_parameters, 1.0)
    pairwise_details: list[dict] = []

    prepared_sources = np.empty_like(normalized, dtype=np.float32)
    for source_index in range(n_channels):
        prepared_sources[source_index] = np.maximum(
            normalized[source_index] - backgrounds[source_index],
            0.0,
        )

    for sink_index in range(n_channels):
        sink_current = normalized[sink_index].copy()
        for source_index in range(n_channels):
            if source_index == sink_index:
                continue
            if source_sink_matrix[source_index, sink_index] != -1:
                continue
            alpha, details = _estimate_source_sink_alpha(
                prepared_sources[source_index],
                sink_current,
                mi_bins=mi_bins,
                alpha_max=alpha_max,
                max_alpha_voxels=max_alpha_voxels,
                random_state=random_state + source_index * 1000 + sink_index,
            )
            sink_current = sink_current - alpha * prepared_sources[source_index]
            alpha_parameters[source_index, sink_index] = float(alpha)
            background_parameters[source_index, sink_index] = float(backgrounds[source_index])
            pairwise_details.append(
                {
                    "sink_index": int(sink_index),
                    "source_index": int(source_index),
                    "alpha": float(alpha),
                    "background": float(backgrounds[source_index]),
                    **details,
                }
            )
        corrected[sink_index] = np.maximum(sink_current, 0.0)

    if global_max > 0.0:
        corrected = corrected * global_max

    details = {
        "background_percentile": float(background_percentile),
        "global_max": float(global_max),
        "background_values": [float(value) for value in backgrounds],
        "mi_bins": int(mi_bins),
        "alpha_max": float(alpha_max),
        "max_alpha_voxels": None if max_alpha_voxels is None else int(max_alpha_voxels),
        "random_state": int(random_state),
        "source_sink_matrix": source_sink_matrix.astype(int).tolist(),
        "alpha_parameters": alpha_parameters.tolist(),
        "background_parameters": background_parameters.tolist(),
        "pairwise_alpha_details": pairwise_details,
    }
    return corrected.astype(np.float32, copy=False), details


def apply_source_sink_parameters(
    channel_volumes,
    *,
    source_sink_matrix: np.ndarray,
    alpha_parameters,
    background_values,
) -> np.ndarray:
    """Apply previously learned source-sink parameters to a new multi-channel volume."""

    channel_volumes = np.asarray(channel_volumes, dtype=np.float32)
    source_sink_matrix = validate_source_sink_matrix(
        source_sink_matrix,
        n_channels=int(channel_volumes.shape[0]),
    )
    alpha_parameters = np.asarray(alpha_parameters, dtype=np.float64)
    background_values = np.asarray(background_values, dtype=np.float64)

    global_max = float(np.max(channel_volumes))
    if global_max > 0.0:
        normalized = channel_volumes / global_max
    else:
        normalized = channel_volumes.copy()

    prepared_sources = np.empty_like(normalized, dtype=np.float32)
    for source_index in range(normalized.shape[0]):
        prepared_sources[source_index] = np.maximum(
            normalized[source_index] - float(background_values[source_index]),
            0.0,
        )

    corrected = normalized.copy()
    for sink_index in range(normalized.shape[0]):
        sink_current = normalized[sink_index].copy()
        for source_index in range(normalized.shape[0]):
            if source_index == sink_index:
                continue
            if source_sink_matrix[source_index, sink_index] != -1:
                continue
            sink_current = sink_current - float(alpha_parameters[source_index, sink_index]) * prepared_sources[source_index]
        corrected[sink_index] = np.maximum(sink_current, 0.0)

    if global_max > 0.0:
        corrected = corrected * global_max
    return corrected.astype(np.float32, copy=False)
# %% END
