"""
Core spectral bleed-through correction routines.

Author: Fabrizio Musacchio
Date: June 2026
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np

from .estimation import (
    DEFAULT_ALPHA_MAX,
    DEFAULT_MAX_ALPHA_VOXELS,
    DEFAULT_MI_BINS,
    DEFAULT_RANDOM_STATE,
    MIN_MASK_VOXELS,
    SUPPORTED_ALPHA_ESTIMATION_METHODS,
    estimate_alpha_from_volume,
    estimate_picasso_unmixing_matrix_from_volume,
)
from .io import CANONICAL_AXIS_ORDER, load_stack_with_omio, write_stack_with_omio

ALPHA_MODES = {"fixed", "reference_t", "per_t"}
SUPPORTED_UNMIX_METHODS = {"manual", *SUPPORTED_ALPHA_ESTIMATION_METHODS}
SUPPORTED_PICASSO_METHODS = {"picasso"}
PICASSO_ALPHA_MODES = {"reference_t", "per_t"}


def report_path_from_output_path(output_path: str | Path) -> Path:
    """Return the JSON sidecar path used for reproducibility metadata."""

    output_path = Path(output_path)
    return output_path.with_suffix(output_path.suffix + ".json")


def _print_verbose(verbose: bool, message: str) -> None:
    if verbose:
        print(message)


def _validate_channel_index(name: str, channel: int, channel_count: int) -> int:
    channel = int(channel)
    if not 0 <= channel < channel_count:
        raise ValueError(
            f"{name} must be between 0 and {channel_count - 1}. Got {channel!r}."
        )
    return channel


def _validate_alpha_mode(alpha_mode: str) -> str:
    alpha_mode = str(alpha_mode).strip().lower()
    if alpha_mode not in ALPHA_MODES:
        raise ValueError(
            f"alpha_mode must be one of {sorted(ALPHA_MODES)}. Got {alpha_mode!r}."
        )
    return alpha_mode


def _validate_picasso_alpha_mode(alpha_mode: str) -> str:
    alpha_mode = _validate_alpha_mode(alpha_mode)
    if alpha_mode not in PICASSO_ALPHA_MODES:
        raise ValueError(
            "unmix_picasso currently supports only alpha_mode='reference_t' or "
            f"'per_t'. Got {alpha_mode!r}."
        )
    return alpha_mode


def _validate_unmix_method(method: str) -> str:
    method = str(method).strip().lower()
    if method not in SUPPORTED_UNMIX_METHODS:
        raise ValueError(
            f"method must be one of {sorted(SUPPORTED_UNMIX_METHODS)}. Got {method!r}."
        )
    return method


def _validate_picasso_method(method: str) -> str:
    method = str(method).strip().lower()
    if method not in SUPPORTED_PICASSO_METHODS:
        raise ValueError(
            f"method must be one of {sorted(SUPPORTED_PICASSO_METHODS)}. Got {method!r}."
        )
    return method


def _validate_alpha_value(alpha: float) -> float:
    alpha = float(alpha)
    if not np.isfinite(alpha) or alpha < 0.0:
        raise ValueError(
            f"alpha must be finite and >= 0 for manual/fixed unmixing. Got {alpha!r}."
        )
    return alpha


def _validate_common_estimation_parameters(
    *,
    alpha_max: float,
    mi_bins: int,
    min_mask_voxels: int,
) -> tuple[float, int, int]:
    alpha_max = float(alpha_max)
    if alpha_max <= 0.0:
        raise ValueError(f"alpha_max must be > 0. Got {alpha_max!r}.")
    mi_bins = int(mi_bins)
    if mi_bins < 2:
        raise ValueError(f"mi_bins must be >= 2. Got {mi_bins!r}.")
    min_mask_voxels = int(min_mask_voxels)
    if min_mask_voxels < 1:
        raise ValueError(f"min_mask_voxels must be >= 1. Got {min_mask_voxels!r}.")
    return alpha_max, mi_bins, min_mask_voxels


def _validate_channels(channels, channel_count: int) -> list[int]:
    if channels is None:
        return list(range(channel_count))

    normalized = [_validate_channel_index("channel", channel, channel_count) for channel in channels]
    if len(normalized) < 2:
        raise ValueError("channels must contain at least two distinct channel indices.")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"channels must not contain duplicates. Got {normalized!r}.")
    return normalized


def _estimate_reference_alpha(
    stack: np.ndarray,
    *,
    alpha_reference_t: int,
    source_channel: int,
    target_channel: int,
    signal_percentile: float,
    target_low_percentile: float | None,
    background_percentile: float,
    method: str,
    preprocess_alpha_inputs: bool,
    alpha_max: float,
    mi_bins: int,
    max_alpha_voxels: int | None,
    random_state: int,
    min_mask_voxels: int,
) -> tuple[float, dict]:
    time_count = int(stack.shape[0])
    if not 0 <= int(alpha_reference_t) < time_count:
        raise ValueError(
            f"alpha_reference_t must be between 0 and {time_count - 1}. "
            f"Got {alpha_reference_t!r}."
        )

    source = stack[int(alpha_reference_t), :, source_channel, :, :]
    target = stack[int(alpha_reference_t), :, target_channel, :, :]
    return estimate_alpha_from_volume(
        source=source,
        target=target,
        signal_percentile=signal_percentile,
        target_low_percentile=target_low_percentile,
        background_percentile=background_percentile,
        min_mask_voxels=min_mask_voxels,
        method=method,
        preprocess_alpha_inputs=preprocess_alpha_inputs,
        alpha_max=alpha_max,
        mi_bins=mi_bins,
        max_alpha_voxels=max_alpha_voxels,
        random_state=random_state,
        return_details=True,
    )


def _estimate_per_t_alphas(
    stack: np.ndarray,
    *,
    source_channel: int,
    target_channel: int,
    signal_percentile: float,
    target_low_percentile: float | None,
    background_percentile: float,
    method: str,
    preprocess_alpha_inputs: bool,
    alpha_max: float,
    mi_bins: int,
    max_alpha_voxels: int | None,
    random_state: int,
    min_mask_voxels: int,
) -> tuple[np.ndarray, list[dict]]:
    alpha_values = np.empty(stack.shape[0], dtype=np.float32)
    details_by_t: list[dict] = []
    for t in range(stack.shape[0]):
        alpha_value, details = estimate_alpha_from_volume(
            source=stack[t, :, source_channel, :, :],
            target=stack[t, :, target_channel, :, :],
            signal_percentile=signal_percentile,
            target_low_percentile=target_low_percentile,
            background_percentile=background_percentile,
            min_mask_voxels=min_mask_voxels,
            method=method,
            preprocess_alpha_inputs=preprocess_alpha_inputs,
            alpha_max=alpha_max,
            mi_bins=mi_bins,
            max_alpha_voxels=max_alpha_voxels,
            random_state=random_state + int(t),
            return_details=True,
        )
        alpha_values[t] = alpha_value
        details_by_t.append({"t": int(t), **details})
    return alpha_values, details_by_t


def _cast_output_stack(stack: np.ndarray, output_dtype: str | np.dtype) -> np.ndarray:
    dtype = np.dtype(output_dtype)
    if np.issubdtype(dtype, np.integer):
        info = np.iinfo(dtype)
        stack = np.clip(stack, info.min, info.max)
    return stack.astype(dtype, copy=False)


def _write_report(report: dict, report_path: Path) -> Path:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    return report_path


def _derive_picasso_output_path(input_path: Path) -> Path:
    return input_path.with_name(f"{input_path.stem}_picasso.tif")


def _move_selected_channels_last(stack: np.ndarray, channels: list[int]) -> np.ndarray:
    selected = np.take(stack, channels, axis=2)
    return np.moveaxis(selected, 2, -1)


def _apply_reference_unmixing_matrix(
    stack: np.ndarray,
    *,
    channels: list[int],
    matrix: np.ndarray,
) -> np.ndarray:
    moved = _move_selected_channels_last(stack, channels)
    unmixed = np.einsum("ij,tzyxj->tzyxi", matrix, moved, optimize=True)
    return np.moveaxis(unmixed, -1, 2)


def _apply_per_t_unmixing_matrices(
    stack: np.ndarray,
    *,
    channels: list[int],
    matrices: np.ndarray,
) -> np.ndarray:
    moved = _move_selected_channels_last(stack, channels)
    unmixed = np.einsum("tij,tzyxj->tzyxi", matrices, moved, optimize=True)
    return np.moveaxis(unmixed, -1, 2)


def _estimate_reference_picasso_matrix(
    stack: np.ndarray,
    *,
    channels: list[int],
    alpha_reference_t: int,
    background_percentile: float,
    preprocess_alpha_inputs: bool,
    mi_bins: int,
    alpha_max: float,
    max_iter: int,
    tolerance: float,
    max_alpha_voxels: int | None,
    random_state: int,
) -> tuple[np.ndarray, dict]:
    time_count = int(stack.shape[0])
    if not 0 <= int(alpha_reference_t) < time_count:
        raise ValueError(
            f"alpha_reference_t must be between 0 and {time_count - 1}. "
            f"Got {alpha_reference_t!r}."
        )

    channel_volumes = np.moveaxis(
        np.take(stack[int(alpha_reference_t), :, :, :, :], channels, axis=1),
        1,
        0,
    )
    matrix, details = estimate_picasso_unmixing_matrix_from_volume(
        channel_volumes,
        background_percentile=background_percentile,
        preprocess_alpha_inputs=preprocess_alpha_inputs,
        mi_bins=mi_bins,
        alpha_max=alpha_max,
        max_iter=max_iter,
        tolerance=tolerance,
        max_alpha_voxels=max_alpha_voxels,
        random_state=random_state,
    )
    return matrix, {"t": int(alpha_reference_t), **details}


def _estimate_per_t_picasso_matrices(
    stack: np.ndarray,
    *,
    channels: list[int],
    background_percentile: float,
    preprocess_alpha_inputs: bool,
    mi_bins: int,
    alpha_max: float,
    max_iter: int,
    tolerance: float,
    max_alpha_voxels: int | None,
    random_state: int,
) -> tuple[np.ndarray, list[dict]]:
    matrices = np.empty((stack.shape[0], len(channels), len(channels)), dtype=np.float64)
    details_by_t: list[dict] = []
    for t in range(stack.shape[0]):
        channel_volumes = np.moveaxis(
            np.take(stack[t, :, :, :, :], channels, axis=1),
            1,
            0,
        )
        matrix, details = estimate_picasso_unmixing_matrix_from_volume(
            channel_volumes,
            background_percentile=background_percentile,
            preprocess_alpha_inputs=preprocess_alpha_inputs,
            mi_bins=mi_bins,
            alpha_max=alpha_max,
            max_iter=max_iter,
            tolerance=tolerance,
            max_alpha_voxels=max_alpha_voxels,
            random_state=random_state + int(t),
        )
        matrices[t] = matrix
        details_by_t.append({"t": int(t), **details})
    return matrices, details_by_t


def unmix(
    input_path,
    output_path,
    alpha=None,
    alpha_mode="fixed",
    alpha_reference_t=0,
    source_channel=0,
    target_channel=1,
    signal_percentile=99.0,
    target_low_percentile=None,
    background_percentile=1.0,
    clip_negative=True,
    output_dtype="float32",
    verbose=True,
    method="mean_ratio",
    preprocess_alpha_inputs=True,
    alpha_max=DEFAULT_ALPHA_MAX,
    mi_bins=DEFAULT_MI_BINS,
    max_alpha_voxels=DEFAULT_MAX_ALPHA_VOXELS,
    random_state=DEFAULT_RANDOM_STATE,
    min_mask_voxels=MIN_MASK_VOXELS,
) -> Path:
    """
    Remove bleed-through from one source channel into one target channel in a TZCYX stack.
    """

    input_path = Path(input_path)
    output_path = Path(output_path)
    report_path = report_path_from_output_path(output_path)

    if input_path.resolve() == output_path.resolve():
        raise ValueError(
            "Refusing to overwrite the input file. Please choose a different output_path."
        )

    alpha_mode = _validate_alpha_mode(alpha_mode)
    method = _validate_unmix_method(method)
    alpha_max, mi_bins, min_mask_voxels = _validate_common_estimation_parameters(
        alpha_max=alpha_max,
        mi_bins=mi_bins,
        min_mask_voxels=min_mask_voxels,
    )

    if method == "manual" and alpha_mode != "fixed":
        raise ValueError("method='manual' is only valid with alpha_mode='fixed'.")

    _print_verbose(verbose, f"Reading stack with OMIO: {input_path}")
    stack, metadata = load_stack_with_omio(input_path)
    channel_count = int(stack.shape[2])
    time_count = int(stack.shape[0])
    z_count = int(stack.shape[1])

    source_channel = _validate_channel_index(
        "source_channel",
        source_channel,
        channel_count,
    )
    target_channel = _validate_channel_index(
        "target_channel",
        target_channel,
        channel_count,
    )
    if source_channel == target_channel:
        raise ValueError("source_channel and target_channel must be different.")

    _print_verbose(
        verbose,
        (
            f"Loaded stack with shape {tuple(int(v) for v in stack.shape)} in "
            f"{CANONICAL_AXIS_ORDER} order. T={'multiple' if time_count > 1 else 'single'} "
            f"({time_count}), Z={'multiple' if z_count > 1 else 'single'} ({z_count})."
        ),
    )

    alpha_scalar: float | None = None
    alpha_values: np.ndarray | None = None
    alpha_details: dict | None = None
    alpha_details_by_t: list[dict] | None = None
    alpha_source = "estimated"
    method_effective = method

    if alpha_mode == "fixed":
        if alpha is None:
            raise ValueError("alpha must be provided when alpha_mode='fixed'.")
        alpha_scalar = _validate_alpha_value(alpha)
        alpha_source = "user_provided"
        method_effective = "manual"
        _print_verbose(
            verbose,
            (
                f"Using user-provided alpha={alpha_scalar:.6f} for source_channel="
                f"{source_channel} -> target_channel={target_channel}."
            ),
        )
    elif alpha_mode == "reference_t":
        alpha_scalar, alpha_details = _estimate_reference_alpha(
            stack,
            alpha_reference_t=int(alpha_reference_t),
            source_channel=source_channel,
            target_channel=target_channel,
            signal_percentile=signal_percentile,
            target_low_percentile=target_low_percentile,
            background_percentile=background_percentile,
            method=method,
            preprocess_alpha_inputs=bool(preprocess_alpha_inputs),
            alpha_max=alpha_max,
            mi_bins=mi_bins,
            max_alpha_voxels=max_alpha_voxels,
            random_state=int(random_state),
            min_mask_voxels=min_mask_voxels,
        )
        _print_verbose(
            verbose,
            (
                f"Estimated reference alpha={alpha_scalar:.6f} from t={int(alpha_reference_t)} "
                f"with method='{method}'."
            ),
        )
    else:
        alpha_values, alpha_details_by_t = _estimate_per_t_alphas(
            stack,
            source_channel=source_channel,
            target_channel=target_channel,
            signal_percentile=signal_percentile,
            target_low_percentile=target_low_percentile,
            background_percentile=background_percentile,
            method=method,
            preprocess_alpha_inputs=bool(preprocess_alpha_inputs),
            alpha_max=alpha_max,
            mi_bins=mi_bins,
            max_alpha_voxels=max_alpha_voxels,
            random_state=int(random_state),
            min_mask_voxels=min_mask_voxels,
        )
        _print_verbose(
            verbose,
            "Estimated per-time-point alpha values: "
            + ", ".join(f"{float(value):.6f}" for value in np.asarray(alpha_values)),
        )

    working_stack = stack.astype(np.float32, copy=True)
    source_view = working_stack[:, :, source_channel, :, :]
    target_view = working_stack[:, :, target_channel, :, :]

    if alpha_mode in {"fixed", "reference_t"}:
        corrected_target = target_view - float(alpha_scalar) * source_view
        if clip_negative:
            corrected_target = np.maximum(corrected_target, 0.0)
        working_stack[:, :, target_channel, :, :] = corrected_target
    else:
        for t in range(working_stack.shape[0]):
            corrected_target = target_view[t] - float(alpha_values[t]) * source_view[t]
            if clip_negative:
                corrected_target = np.maximum(corrected_target, 0.0)
            working_stack[t, :, target_channel, :, :] = corrected_target

    _print_verbose(verbose, f"Writing corrected stack to: {output_path}")
    output_stack = _cast_output_stack(working_stack, output_dtype)
    actual_output_path = write_stack_with_omio(output_path, output_stack, metadata)

    report = {
        "input_path": str(input_path),
        "output_path": str(actual_output_path),
        "report_path": str(report_path),
        "alpha_mode": alpha_mode,
        "method": method,
        "method_effective": method_effective,
        "alpha_source": alpha_source,
        "alpha": None if alpha_scalar is None else float(alpha_scalar),
        "alpha_values": None
        if alpha_values is None
        else [float(value) for value in np.asarray(alpha_values)],
        "alpha_by_t": None
        if alpha_values is None
        else [float(value) for value in np.asarray(alpha_values)],
        "source_channel": int(source_channel),
        "target_channel": int(target_channel),
        "signal_percentile": float(signal_percentile),
        "target_low_percentile": None
        if target_low_percentile is None
        else float(target_low_percentile),
        "background_percentile": float(background_percentile),
        "preprocess_alpha_inputs": bool(preprocess_alpha_inputs),
        "alpha_max": float(alpha_max),
        "mi_bins": int(mi_bins),
        "max_alpha_voxels": None if max_alpha_voxels is None else int(max_alpha_voxels),
        "random_state": int(random_state),
        "min_mask_voxels": int(min_mask_voxels),
        "mask_voxel_count": None
        if alpha_details is None
        else int(alpha_details["mask_voxel_count"]),
        "mask_voxel_count_by_t": None
        if alpha_details_by_t is None
        else [int(item["mask_voxel_count"]) for item in alpha_details_by_t],
        "alpha_estimation": alpha_details,
        "alpha_estimation_by_t": alpha_details_by_t,
        "input_shape": tuple(int(v) for v in stack.shape),
        "axis_order": CANONICAL_AXIS_ORDER,
        "size_t": time_count,
        "size_z": z_count,
        "has_multiple_t": bool(time_count > 1),
        "has_multiple_z": bool(z_count > 1),
        "clip_negative": bool(clip_negative),
        "output_dtype": str(np.dtype(output_dtype)),
    }
    actual_report_path = _write_report(report, report_path)
    _print_verbose(verbose, f"Wrote processing report to: {actual_report_path}")
    return actual_output_path


def unmix_picasso(
    input_path,
    output_path=None,
    channels=None,
    *,
    method="picasso",
    alpha_mode="reference_t",
    alpha_reference_t=0,
    background_percentile=1.0,
    preprocess_alpha_inputs=True,
    mi_bins=DEFAULT_MI_BINS,
    alpha_max=DEFAULT_ALPHA_MAX,
    max_iter=10,
    tolerance=1e-4,
    max_alpha_voxels=DEFAULT_MAX_ALPHA_VOXELS,
    random_state=DEFAULT_RANDOM_STATE,
    clip_negative=True,
    output_dtype="float32",
    verbose=True,
) -> Path:
    """
    Perform PICASSO-like iterative multi-channel blind unmixing.

    Notes
    -----
    This function assumes that the number of measured channels equals the number
    of fluorophores to reconstruct.
    """

    input_path = Path(input_path)
    output_path = _derive_picasso_output_path(input_path) if output_path is None else Path(output_path)
    report_path = report_path_from_output_path(output_path)

    if input_path.resolve() == output_path.resolve():
        raise ValueError(
            "Refusing to overwrite the input file. Please choose a different output_path."
        )

    method = _validate_picasso_method(method)
    alpha_mode = _validate_picasso_alpha_mode(alpha_mode)
    alpha_max, mi_bins, _ = _validate_common_estimation_parameters(
        alpha_max=alpha_max,
        mi_bins=mi_bins,
        min_mask_voxels=1,
    )

    _print_verbose(verbose, f"Reading stack with OMIO: {input_path}")
    stack, metadata = load_stack_with_omio(input_path)
    time_count = int(stack.shape[0])
    channels = _validate_channels(channels, int(stack.shape[2]))

    if alpha_mode == "reference_t":
        matrix, matrix_details = _estimate_reference_picasso_matrix(
            stack,
            channels=channels,
            alpha_reference_t=int(alpha_reference_t),
            background_percentile=background_percentile,
            preprocess_alpha_inputs=bool(preprocess_alpha_inputs),
            mi_bins=mi_bins,
            alpha_max=alpha_max,
            max_iter=int(max_iter),
            tolerance=float(tolerance),
            max_alpha_voxels=max_alpha_voxels,
            random_state=int(random_state),
        )
        transformed_selected = _apply_reference_unmixing_matrix(
            stack.astype(np.float32, copy=False),
            channels=channels,
            matrix=matrix,
        )
        matrices_by_t = None
        matrix_details_by_t = None
    else:
        matrices_by_t, matrix_details_by_t = _estimate_per_t_picasso_matrices(
            stack,
            channels=channels,
            background_percentile=background_percentile,
            preprocess_alpha_inputs=bool(preprocess_alpha_inputs),
            mi_bins=mi_bins,
            alpha_max=alpha_max,
            max_iter=int(max_iter),
            tolerance=float(tolerance),
            max_alpha_voxels=max_alpha_voxels,
            random_state=int(random_state),
        )
        transformed_selected = _apply_per_t_unmixing_matrices(
            stack.astype(np.float32, copy=False),
            channels=channels,
            matrices=matrices_by_t,
        )
        matrix = None
        matrix_details = None

    if clip_negative:
        transformed_selected = np.maximum(transformed_selected, 0.0)

    working_stack = stack.astype(np.float32, copy=True)
    working_stack[:, :, channels, :, :] = transformed_selected

    _print_verbose(verbose, f"Writing blind-unmixed stack to: {output_path}")
    output_stack = _cast_output_stack(working_stack, output_dtype)
    actual_output_path = write_stack_with_omio(output_path, output_stack, metadata)

    report = {
        "input_path": str(input_path),
        "output_path": str(actual_output_path),
        "report_path": str(report_path),
        "method": method,
        "alpha_mode": alpha_mode,
        "alpha_reference_t": int(alpha_reference_t),
        "channels": [int(channel) for channel in channels],
        "background_percentile": float(background_percentile),
        "preprocess_alpha_inputs": bool(preprocess_alpha_inputs),
        "alpha_max": float(alpha_max),
        "mi_bins": int(mi_bins),
        "max_iter": int(max_iter),
        "tolerance": float(tolerance),
        "max_alpha_voxels": None if max_alpha_voxels is None else int(max_alpha_voxels),
        "random_state": int(random_state),
        "clip_negative": bool(clip_negative),
        "unmixing_matrix": None if matrix is None else matrix.tolist(),
        "unmixing_matrix_by_t": None
        if matrices_by_t is None
        else matrices_by_t.tolist(),
        "iterations_run": None
        if matrix_details is None
        else int(matrix_details["iterations_run"]),
        "iterations_run_by_t": None
        if matrix_details_by_t is None
        else [int(item["iterations_run"]) for item in matrix_details_by_t],
        "converged": None
        if matrix_details is None
        else bool(matrix_details["converged"]),
        "converged_by_t": None
        if matrix_details_by_t is None
        else [bool(item["converged"]) for item in matrix_details_by_t],
        "picasso_estimation": matrix_details,
        "picasso_estimation_by_t": matrix_details_by_t,
        "input_shape": tuple(int(v) for v in stack.shape),
        "axis_order": CANONICAL_AXIS_ORDER,
        "size_t": int(stack.shape[0]),
        "size_z": int(stack.shape[1]),
        "output_dtype": str(np.dtype(output_dtype)),
    }
    actual_report_path = _write_report(report, report_path)
    _print_verbose(verbose, f"Wrote processing report to: {actual_report_path}")
    return actual_output_path


def unmix_ch0_from_ch1(*args, **kwargs) -> Path:
    """
    Backward-compatible wrapper for older code paths.

    Deprecated in favor of :func:`unmix`.
    """

    warnings.warn(
        "unmix_ch0_from_ch1 is deprecated; use unmix with source_channel and "
        "target_channel instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return unmix(*args, **kwargs)
