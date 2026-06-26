"""
Core spectral bleed-through correction routines.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
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
from .picasso_impl import (
    build_source_sink_matrix,
    DEFAULT_PICASSO_ALPHA_CLIP,
    DEFAULT_PICASSO_BIN_FACTOR,
    DEFAULT_PICASSO_CLIP_EVERY_N_ITERATIONS,
    DEFAULT_PICASSO_IMPLEMENTATION,
    DEFAULT_PICASSO_NEGATIVITY_THRESHOLD,
    DEFAULT_PICASSO_QN,
    DEFAULT_PICASSO_STEP_SIZE,
    default_source_sink_matrix,
    apply_source_sink_parameters,
    run_picasso_matlab_like,
    run_source_sink_unmixing,
    apply_matlab_incremental_sequence,
    validate_picasso_implementation,
    validate_source_sink_matrix,
)
# %% CONSTANTS
ALPHA_MODES = {"fixed", "reference_t", "per_t"}
SUPPORTED_UNMIX_METHODS = {"manual", *SUPPORTED_ALPHA_ESTIMATION_METHODS}
SUPPORTED_PICASSO_METHODS = {"picasso"}
PICASSO_ALPHA_MODES = {"reference_t", "per_t"}

# %% INTERNAL HELPERS
def report_path_from_output_path(output_path: str | Path) -> Path:
    """Return the JSON sidecar path used for reproducibility metadata."""

    output_path = Path(output_path)
    return output_path.with_suffix(output_path.suffix + ".json")


def _print_verbose(verbose: bool, message: str) -> None:
    """Print a status message only when verbose mode is enabled."""

    if verbose:
        print(message)


def _validate_channel_index(name: str, channel: int, channel_count: int) -> int:
    """Validate a channel index against the available channel count."""

    channel = int(channel)
    if not 0 <= channel < channel_count:
        raise ValueError(
            f"{name} must be between 0 and {channel_count - 1}. Got {channel!r}."
        )
    return channel


def _validate_alpha_mode(alpha_mode: str) -> str:
    """Normalize and validate a scalar-alpha application mode."""

    alpha_mode = str(alpha_mode).strip().lower()
    if alpha_mode not in ALPHA_MODES:
        raise ValueError(
            f"alpha_mode must be one of {sorted(ALPHA_MODES)}. Got {alpha_mode!r}."
        )
    return alpha_mode


def _validate_picasso_alpha_mode(alpha_mode: str) -> str:
    """Validate the subset of alpha modes currently supported by ``unmix_picasso``."""

    alpha_mode = _validate_alpha_mode(alpha_mode)
    if alpha_mode not in PICASSO_ALPHA_MODES:
        raise ValueError(
            "unmix_picasso currently supports only alpha_mode='reference_t' or "
            f"'per_t'. Got {alpha_mode!r}."
        )
    return alpha_mode


def _validate_unmix_method(method: str) -> str:
    """Normalize and validate a two-channel unmixing method name."""

    method = str(method).strip().lower()
    if method not in SUPPORTED_UNMIX_METHODS:
        raise ValueError(
            f"method must be one of {sorted(SUPPORTED_UNMIX_METHODS)}. Got {method!r}."
        )
    return method


def _validate_picasso_method(method: str) -> str:
    """Normalize and validate the multi-channel blind-unmixing method name."""

    method = str(method).strip().lower()
    if method not in SUPPORTED_PICASSO_METHODS:
        raise ValueError(
            f"method must be one of {sorted(SUPPORTED_PICASSO_METHODS)}. Got {method!r}."
        )
    return method


def _validate_alpha_value(alpha: float) -> float:
    """Validate a user-provided fixed alpha value."""

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
    """Validate alpha-estimation parameters shared across multiple workflows."""

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
    """Normalize and validate the channel subset used for PICASSO-style unmixing."""

    if channels is None:
        return list(range(channel_count))

    normalized = [_validate_channel_index("channel", channel, channel_count) for channel in channels]
    if len(normalized) < 2:
        raise ValueError("channels must contain at least two distinct channel indices.")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"channels must not contain duplicates. Got {normalized!r}.")
    return normalized


def _resolve_reverse_parameter(forward_value, reverse_value):
    """Return a reverse-direction override or fall back to the forward value."""

    return forward_value if reverse_value is None else reverse_value


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
    """Estimate one scalar alpha from all Z slices of a chosen reference time point."""

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
    """Estimate one scalar alpha per time point using all Z slices of each time point."""

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


def _estimate_directional_alpha(
    stack: np.ndarray,
    *,
    alpha_mode: str,
    alpha,
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
    verbose: bool,
    direction_label: str,
) -> tuple[float | None, np.ndarray | None, dict | None, list[dict] | None, str, str]:
    """Estimate or validate one directional bleed-through coefficient workflow."""

    alpha_scalar: float | None = None
    alpha_values: np.ndarray | None = None
    alpha_details: dict | None = None
    alpha_details_by_t: list[dict] | None = None
    alpha_source = "estimated"
    method_effective = method

    if alpha_mode == "fixed":
        if alpha is None:
            raise ValueError(
                f"alpha must be provided for {direction_label} correction when alpha_mode='fixed'."
            )
        alpha_scalar = _validate_alpha_value(alpha)
        alpha_source = "user_provided"
        method_effective = "manual"
        _print_verbose(
            verbose,
            (
                f"Using user-provided {direction_label} alpha={alpha_scalar:.6f} for "
                f"source_channel={source_channel} -> target_channel={target_channel}."
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
                f"Estimated {direction_label} reference alpha={alpha_scalar:.6f} from "
                f"t={int(alpha_reference_t)} with method='{method}'."
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
            f"Estimated {direction_label} per-time-point alpha values: "
            + ", ".join(f"{float(value):.6f}" for value in np.asarray(alpha_values)),
        )

    return (
        alpha_scalar,
        alpha_values,
        alpha_details,
        alpha_details_by_t,
        alpha_source,
        method_effective,
    )


def _validate_bidirectional_determinant(
    alpha_forward,
    alpha_reverse,
    *,
    context: str,
) -> float:
    """Validate the determinant of the 2x2 bidirectional mixing matrix."""

    determinant = 1.0 - float(alpha_forward) * float(alpha_reverse)
    if determinant <= 0.0 or not np.isfinite(determinant):
        raise ValueError(
            f"Invalid bidirectional unmixing matrix for {context}: "
            f"determinant = 1 - alpha_forward * alpha_reverse = {determinant!r}. "
            "Please choose coefficients with alpha_forward * alpha_reverse < 1."
        )
    return float(determinant)


def _apply_bidirectional_unmixing(
    source_measured: np.ndarray,
    target_measured: np.ndarray,
    *,
    alpha_forward: float,
    alpha_reverse: float,
    clip_negative: bool,
    context: str,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Unmix a two-channel bidirectional linear mixture by 2x2 matrix inversion."""

    determinant = _validate_bidirectional_determinant(
        alpha_forward,
        alpha_reverse,
        context=context,
    )
    corrected_source = (
        source_measured - float(alpha_reverse) * target_measured
    ) / determinant
    corrected_target = (
        target_measured - float(alpha_forward) * source_measured
    ) / determinant

    if clip_negative:
        corrected_source = np.maximum(corrected_source, 0.0)
        corrected_target = np.maximum(corrected_target, 0.0)

    return corrected_source, corrected_target, determinant


def _cast_output_stack(stack: np.ndarray, output_dtype: str | np.dtype) -> np.ndarray:
    """Cast an output stack while clipping integer targets to their valid range."""

    dtype = np.dtype(output_dtype)
    if np.issubdtype(dtype, np.integer):
        info = np.iinfo(dtype)
        stack = np.clip(stack, info.min, info.max)
    return stack.astype(dtype, copy=False)


def _write_report(report: dict, report_path: Path) -> Path:
    """Write a JSON sidecar report for a processing run."""

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    return report_path


def _derive_picasso_output_path(input_path: Path) -> Path:
    """Derive a default output path for PICASSO-style unmixing results."""

    return input_path.with_name(f"{input_path.stem}_picasso.tif")


def _move_selected_channels_last(stack: np.ndarray, channels: list[int]) -> np.ndarray:
    """Extract selected channels and move channel to the trailing axis."""

    selected = np.take(stack, channels, axis=2)
    return np.moveaxis(selected, 2, -1)


def _apply_reference_unmixing_matrix(
    stack: np.ndarray,
    *,
    channels: list[int],
    matrix: np.ndarray,
) -> np.ndarray:
    """Apply one blind-unmixing matrix to the selected channels of an entire stack."""

    moved = _move_selected_channels_last(stack, channels)
    unmixed = np.einsum("ij,tzyxj->tzyxi", matrix, moved, optimize=True)
    return np.moveaxis(unmixed, -1, 2)


def _apply_per_t_unmixing_matrices(
    stack: np.ndarray,
    *,
    channels: list[int],
    matrices: np.ndarray,
) -> np.ndarray:
    """Apply one blind-unmixing matrix per time point to selected stack channels."""

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
    """Estimate one blind-unmixing matrix from a chosen reference time point."""

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
    """Estimate one blind-unmixing matrix per time point."""

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
    bidirectional=False,
    alpha_reverse=None,
    method_reverse=None,
    preprocess_alpha_inputs=True,
    alpha_max=DEFAULT_ALPHA_MAX,
    signal_percentile_reverse=None,
    background_percentile_reverse=None,
    target_low_percentile_reverse=None,
    alpha_max_reverse=None,
    mi_bins=DEFAULT_MI_BINS,
    mi_bins_reverse=None,
    max_alpha_voxels=DEFAULT_MAX_ALPHA_VOXELS,
    random_state=DEFAULT_RANDOM_STATE,
    min_mask_voxels=MIN_MASK_VOXELS,
) -> Path:
    """
    Remove bleed-through from one source channel into one target channel in a TZCYX stack.

    Parameters
    ----------
    input_path : str or Path
        Path to the input TIFF stack.
    output_path : str or Path
        Path to the output TIFF stack. A JSON sidecar report with the same name
        plus ``.json`` is written alongside it.
    alpha : float or None, optional
        User-provided bleed-through coefficient. Required when
        ``alpha_mode="fixed"``.
    alpha_mode : {"fixed", "reference_t", "per_t"}, optional
        Strategy that determines whether alpha is taken directly from ``alpha``,
        estimated once from a reference time point, or estimated separately for
        each time point.
    alpha_reference_t : int, optional
        Reference time point used when ``alpha_mode="reference_t"``.
    source_channel : int, optional
        Channel whose signal bleeds into the target channel.
    target_channel : int, optional
        Channel from which the source contribution should be removed.
    signal_percentile : float, optional
        Percentile used to define a bright-source signal mask for automatic
        alpha estimation.
    target_low_percentile : float or None, optional
        Optional low-target constraint for the alpha-estimation mask.
    background_percentile : float, optional
        Low percentile used for optional percentile-based background subtraction
        during alpha estimation.
    clip_negative : bool, optional
        If ``True``, clip negative corrected target values to zero.
    output_dtype : str or numpy.dtype, optional
        Data type used for the written output stack.
    verbose : bool, optional
        If ``True``, print processing progress and estimated coefficients.
    method : {"manual", "mean_ratio", "linear_fit", "corr_min", "mi_min"}, optional
        Method used to obtain alpha. ``"manual"`` is meaningful only together
        with ``alpha_mode="fixed"``; the other methods estimate alpha from the
        data.
    bidirectional : bool, optional
        If ``True``, estimate or use coefficients for both
        ``source_channel -> target_channel`` and
        ``target_channel -> source_channel`` and solve the resulting 2x2 linear
        mixing model by matrix inversion.
    alpha_reverse : float or None, optional
        Optional fixed reverse-direction bleed-through coefficient, i.e. from
        ``target_channel`` back into ``source_channel``. If ``None`` and
        ``bidirectional=True``, the forward ``alpha`` value is reused.
    method_reverse : {"manual", "mean_ratio", "linear_fit", "corr_min", "mi_min"} or None, optional
        Optional reverse-direction alpha-estimation method. If ``None`` and
        ``bidirectional=True``, the forward ``method`` value is reused.
    preprocess_alpha_inputs : bool, optional
        If ``True``, apply percentile-based background subtraction and clipping
        before automatic alpha estimation.
    alpha_max : float, optional
        Upper search bound for optimization-based alpha-estimation methods.
    signal_percentile_reverse : float or None, optional
        Optional reverse-direction source-mask percentile. Falls back to
        ``signal_percentile`` when ``None``.
    background_percentile_reverse : float or None, optional
        Optional reverse-direction background percentile. Falls back to
        ``background_percentile`` when ``None``.
    target_low_percentile_reverse : float or None, optional
        Optional reverse-direction low-target percentile. Falls back to
        ``target_low_percentile`` when ``None``.
    alpha_max_reverse : float or None, optional
        Optional reverse-direction optimization bound. Falls back to
        ``alpha_max`` when ``None``.
    mi_bins : int, optional
        Number of histogram bins used by the mutual-information estimator.
    mi_bins_reverse : int or None, optional
        Optional reverse-direction histogram-bin count used by ``mi_min``.
        Falls back to ``mi_bins`` when ``None``.
    max_alpha_voxels : int or None, optional
        Optional cap on the number of voxels used for alpha estimation after
        masking.
    random_state : int, optional
        Random seed used for optional voxel subsampling during alpha estimation.
    min_mask_voxels : int, optional
        Minimum number of voxels required for a valid alpha-estimation mask.

    Returns
    -------
    Path
        Actual path of the written TIFF stack.

    Raises
    ------
    ValueError
        If the input configuration is inconsistent or would overwrite the input.

    Notes
    -----
    Only ``target_channel`` is modified. ``source_channel`` remains unchanged in
    the output stack. Automatic alpha estimation is performed on prepared data,
    but the final subtraction is applied to the original stack intensities cast
    to ``float32``.

    If ``bidirectional=True``, both selected channels are corrected jointly by
    inverting the 2x2 linear mixing model. In that mode, ``clip_negative``
    applies to both corrected channels.
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
    bidirectional = bool(bidirectional)
    method_reverse_resolved = None if not bidirectional else _validate_unmix_method(
        _resolve_reverse_parameter(method, method_reverse)
    )
    signal_percentile_reverse_resolved = _resolve_reverse_parameter(
        signal_percentile,
        signal_percentile_reverse,
    )
    background_percentile_reverse_resolved = _resolve_reverse_parameter(
        background_percentile,
        background_percentile_reverse,
    )
    target_low_percentile_reverse_resolved = _resolve_reverse_parameter(
        target_low_percentile,
        target_low_percentile_reverse,
    )
    alpha_max_reverse_resolved = _resolve_reverse_parameter(
        alpha_max,
        alpha_max_reverse,
    )
    mi_bins_reverse_resolved = _resolve_reverse_parameter(
        mi_bins,
        mi_bins_reverse,
    )
    if bidirectional:
        alpha_max_reverse_resolved, mi_bins_reverse_resolved, _ = _validate_common_estimation_parameters(
            alpha_max=alpha_max_reverse_resolved,
            mi_bins=mi_bins_reverse_resolved,
            min_mask_voxels=min_mask_voxels,
        )

    if method == "manual" and alpha_mode != "fixed":
        raise ValueError("method='manual' is only valid with alpha_mode='fixed'.")
    if bidirectional and method_reverse_resolved == "manual" and alpha_mode != "fixed":
        raise ValueError("method_reverse='manual' is only valid with alpha_mode='fixed'.")

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
    (
        alpha_scalar,
        alpha_values,
        alpha_details,
        alpha_details_by_t,
        alpha_source,
        method_effective,
    ) = _estimate_directional_alpha(
        stack,
        alpha_mode=alpha_mode,
        alpha=alpha,
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
        verbose=bool(verbose),
        direction_label="forward",
    )

    alpha_reverse_scalar: float | None = None
    alpha_reverse_values: np.ndarray | None = None
    alpha_reverse_details: dict | None = None
    alpha_reverse_details_by_t: list[dict] | None = None
    alpha_reverse_source: str | None = None
    method_reverse_effective: str | None = None

    if bidirectional:
        reverse_alpha_input = _resolve_reverse_parameter(alpha, alpha_reverse)
        (
            alpha_reverse_scalar,
            alpha_reverse_values,
            alpha_reverse_details,
            alpha_reverse_details_by_t,
            alpha_reverse_source,
            method_reverse_effective,
        ) = _estimate_directional_alpha(
            stack,
            alpha_mode=alpha_mode,
            alpha=reverse_alpha_input,
            alpha_reference_t=int(alpha_reference_t),
            source_channel=target_channel,
            target_channel=source_channel,
            signal_percentile=signal_percentile_reverse_resolved,
            target_low_percentile=target_low_percentile_reverse_resolved,
            background_percentile=background_percentile_reverse_resolved,
            method=method_reverse_resolved,
            preprocess_alpha_inputs=bool(preprocess_alpha_inputs),
            alpha_max=alpha_max_reverse_resolved,
            mi_bins=mi_bins_reverse_resolved,
            max_alpha_voxels=max_alpha_voxels,
            random_state=int(random_state) + 10_000,
            min_mask_voxels=min_mask_voxels,
            verbose=bool(verbose),
            direction_label="reverse",
        )

    working_stack = stack.astype(np.float32, copy=True)
    source_view = working_stack[:, :, source_channel, :, :]
    target_view = working_stack[:, :, target_channel, :, :]
    bidirectional_determinants: list[float] | None = None

    if not bidirectional and alpha_mode in {"fixed", "reference_t"}:
        corrected_target = target_view - float(alpha_scalar) * source_view
        if clip_negative:
            corrected_target = np.maximum(corrected_target, 0.0)
        working_stack[:, :, target_channel, :, :] = corrected_target
    elif not bidirectional:
        for t in range(working_stack.shape[0]):
            corrected_target = target_view[t] - float(alpha_values[t]) * source_view[t]
            if clip_negative:
                corrected_target = np.maximum(corrected_target, 0.0)
            working_stack[t, :, target_channel, :, :] = corrected_target
    elif alpha_mode in {"fixed", "reference_t"}:
        corrected_source, corrected_target, determinant = _apply_bidirectional_unmixing(
            source_view,
            target_view,
            alpha_forward=float(alpha_scalar),
            alpha_reverse=float(alpha_reverse_scalar),
            clip_negative=bool(clip_negative),
            context="global bidirectional unmixing",
        )
        working_stack[:, :, source_channel, :, :] = corrected_source
        working_stack[:, :, target_channel, :, :] = corrected_target
        bidirectional_determinants = [float(determinant)]
    else:
        bidirectional_determinants = []
        for t in range(working_stack.shape[0]):
            corrected_source, corrected_target, determinant = _apply_bidirectional_unmixing(
                source_view[t],
                target_view[t],
                alpha_forward=float(alpha_values[t]),
                alpha_reverse=float(alpha_reverse_values[t]),
                clip_negative=bool(clip_negative),
                context=f"bidirectional unmixing at t={t}",
            )
            working_stack[t, :, source_channel, :, :] = corrected_source
            working_stack[t, :, target_channel, :, :] = corrected_target
            bidirectional_determinants.append(float(determinant))

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
        "bidirectional": bool(bidirectional),
        "alpha_reverse": None
        if alpha_reverse_scalar is None
        else float(alpha_reverse_scalar),
        "alpha_reverse_values": None
        if alpha_reverse_values is None
        else [float(value) for value in np.asarray(alpha_reverse_values)],
        "alpha_reverse_by_t": None
        if alpha_reverse_values is None
        else [float(value) for value in np.asarray(alpha_reverse_values)],
        "alpha_reverse_inherited_from_forward": bool(bidirectional and alpha_reverse is None),
        "method_reverse": None if not bidirectional else method_reverse_resolved,
        "method_reverse_effective": method_reverse_effective,
        "alpha_reverse_source": alpha_reverse_source,
        "source_channel": int(source_channel),
        "target_channel": int(target_channel),
        "signal_percentile": float(signal_percentile),
        "target_low_percentile": None
        if target_low_percentile is None
        else float(target_low_percentile),
        "background_percentile": float(background_percentile),
        "signal_percentile_reverse": None
        if not bidirectional
        else float(signal_percentile_reverse_resolved),
        "target_low_percentile_reverse": None
        if not bidirectional or target_low_percentile_reverse_resolved is None
        else float(target_low_percentile_reverse_resolved),
        "background_percentile_reverse": None
        if not bidirectional
        else float(background_percentile_reverse_resolved),
        "preprocess_alpha_inputs": bool(preprocess_alpha_inputs),
        "alpha_max": float(alpha_max),
        "alpha_max_reverse": None if not bidirectional else float(alpha_max_reverse_resolved),
        "mi_bins": int(mi_bins),
        "mi_bins_reverse": None if not bidirectional else int(mi_bins_reverse_resolved),
        "max_alpha_voxels": None if max_alpha_voxels is None else int(max_alpha_voxels),
        "random_state": int(random_state),
        "min_mask_voxels": int(min_mask_voxels),
        "mask_voxel_count": None
        if alpha_details is None
        else int(alpha_details["mask_voxel_count"]),
        "mask_voxel_count_reverse": None
        if alpha_reverse_details is None
        else int(alpha_reverse_details["mask_voxel_count"]),
        "mask_voxel_count_by_t": None
        if alpha_details_by_t is None
        else [int(item["mask_voxel_count"]) for item in alpha_details_by_t],
        "mask_voxel_count_reverse_by_t": None
        if alpha_reverse_details_by_t is None
        else [int(item["mask_voxel_count"]) for item in alpha_reverse_details_by_t],
        "alpha_estimation": alpha_details,
        "alpha_estimation_by_t": alpha_details_by_t,
        "alpha_estimation_reverse": alpha_reverse_details,
        "alpha_estimation_reverse_by_t": alpha_reverse_details_by_t,
        "bidirectional_mixing_matrix_determinant": None
        if bidirectional_determinants is None
        else (
            float(bidirectional_determinants[0])
            if len(bidirectional_determinants) == 1
            else None
        ),
        "bidirectional_mixing_matrix_determinant_by_t": None
        if bidirectional_determinants is None or len(bidirectional_determinants) == 1
        else [float(value) for value in bidirectional_determinants],
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
    implementation=DEFAULT_PICASSO_IMPLEMENTATION,
    alpha_mode="reference_t",
    alpha_reference_t=0,
    source_sink_matrix=None,
    sink_channels=None,
    neutral_channels=None,
    background_percentile=1.0,
    preprocess_alpha_inputs=True,
    mi_bins=DEFAULT_MI_BINS,
    alpha_max=DEFAULT_ALPHA_MAX,
    max_iter=200,
    tolerance=1e-4,
    max_alpha_voxels=DEFAULT_MAX_ALPHA_VOXELS,
    random_state=DEFAULT_RANDOM_STATE,
    step_size=DEFAULT_PICASSO_STEP_SIZE,
    qn=DEFAULT_PICASSO_QN,
    pixel_bin_size=DEFAULT_PICASSO_BIN_FACTOR,
    alpha_clip=DEFAULT_PICASSO_ALPHA_CLIP,
    negativity_threshold=DEFAULT_PICASSO_NEGATIVITY_THRESHOLD,
    clip_every_n_iterations=DEFAULT_PICASSO_CLIP_EVERY_N_ITERATIONS,
    clip_negative=True,
    output_dtype="float32",
    verbose=True,
) -> Path:
    """
    Perform PICASSO-family multi-channel blind unmixing.

    Parameters
    ----------
    input_path : str or Path
        Path to the input TIFF stack.
    output_path : str or Path or None, optional
        Output TIFF path. If ``None``, a filename ending in ``"_picasso.tif"``
        is created next to the input.
    channels : sequence of int or None, optional
        Channel indices to include in blind unmixing. If ``None``, all channels
        are used.
    method : {"picasso"}, optional
        Method label reserved for the PICASSO-like workflow.
    implementation : {"matlab_3c", "matlab_n", "source_sink_n"}, optional
        Concrete implementation variant:
        ``"matlab_3c"`` for a close Python port of the original MATLAB 3-channel
        workflow, ``"matlab_n"`` for an N-channel generalization of that
        workflow, and ``"source_sink_n"`` for a source-sink multi-channel model
        inspired by the napari plugin.
    alpha_mode : {"reference_t", "per_t"}, optional
        Whether to estimate one unmixing matrix from a reference time point or
        one matrix per time point.
    alpha_reference_t : int, optional
        Reference time point used when ``alpha_mode="reference_t"``.
    source_sink_matrix : array-like or None, optional
        Square matrix used only by ``implementation="source_sink_n"``. The
        selected channels define both rows and columns. The diagonal must be
        ``1``. Off-diagonal values must be ``-1`` for modeled source-to-sink
        spillover and ``0`` otherwise. If ``None``, an all-to-all matrix with
        diagonal ``1`` and off-diagonal ``-1`` is used.
    sink_channels : sequence of int or None, optional
        Convenience alternative to ``source_sink_matrix`` for
        ``implementation="source_sink_n"``. Provide actual channel indices from
        ``channels`` that should be corrected as sinks. All selected channels
        except those listed in ``neutral_channels`` may contribute to these
        sinks.
    neutral_channels : sequence of int or None, optional
        Convenience option for ``implementation="source_sink_n"``. Provide
        actual channel indices from ``channels`` that should remain neutral,
        meaning they are neither corrected as sinks nor used as sources when a
        source-sink matrix is auto-generated.
    background_percentile : float, optional
        Low percentile used for per-channel background subtraction in the
        selected PICASSO implementation.
    preprocess_alpha_inputs : bool, optional
        Backward-compatibility parameter retained from earlier PICASSO-like
        implementations. The current PICASSO-family implementations define
        their own internal preprocessing and therefore do not switch behavior
        based on this flag. The value is nevertheless recorded in the JSON
        sidecar report for reproducibility.
    mi_bins : int, optional
        Number of histogram bins used by the mutual-information estimator in the
        source-sink implementation.
    alpha_max : float, optional
        Upper bound for source-to-sink coefficients in the source-sink
        implementation.
    max_iter : int, optional
        Maximum number of iterations for the selected implementation.
    tolerance : float, optional
        Reserved for compatibility with the previous API. The MATLAB-like
        implementations follow a fixed iteration count rather than tolerance-
        based early stopping.
    max_alpha_voxels : int or None, optional
        Optional cap on the number of voxels used for coefficient estimation in
        the source-sink implementation.
    random_state : int, optional
        Random seed used for optional subsampling during matrix estimation.
    step_size : float, optional
        Incremental update step size used by the MATLAB-style implementations.
    qn : int, optional
        MATLAB-style mutual-information quantization parameter.
    pixel_bin_size : int, optional
        MATLAB-style 2D binning factor applied before mutual-information
        evaluation.
    alpha_clip : float, optional
        Absolute clipping bound applied to each MATLAB-style pairwise alpha
        estimate before the incremental update matrix is constructed.
    negativity_threshold : float, optional
        MATLAB-style negativity ratio threshold above which a pairwise alpha is
        divided by ten.
    clip_every_n_iterations : int, optional
        Frequency of positivity enforcement during the MATLAB-style iterations.
    clip_negative : bool, optional
        If ``True``, clip negative unmixed intensities to zero before writing.
    output_dtype : str or numpy.dtype, optional
        Data type used for the written output stack.
    verbose : bool, optional
        If ``True``, print processing progress and output paths.

    Returns
    -------
    Path
        Actual path of the written TIFF stack.

    Notes
    -----
    ``implementation="matlab_3c"`` is the default and aims to match the
    published MATLAB 3-channel workflow as closely as possible.

    ``implementation="matlab_n"`` is an explicit N-channel generalization of
    that MATLAB workflow.

    ``implementation="source_sink_n"`` uses a direct source-sink subtraction
    model inspired by the napari plugin structure. It is not a neural-MINE port
    of the plugin, but it follows the same idea of user-defined source/sink
    relations.
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
    implementation = validate_picasso_implementation(implementation)
    alpha_max, mi_bins, _ = _validate_common_estimation_parameters(
        alpha_max=alpha_max,
        mi_bins=mi_bins,
        min_mask_voxels=1,
    )

    _print_verbose(verbose, f"Reading stack with OMIO: {input_path}")
    stack, metadata = load_stack_with_omio(input_path)
    time_count = int(stack.shape[0])
    channels = _validate_channels(channels, int(stack.shape[2]))
    selected_stack = np.take(stack, channels, axis=2).astype(np.float32, copy=False)

    if implementation == "matlab_3c" and len(channels) != 3:
        raise ValueError(
            "implementation='matlab_3c' requires exactly 3 selected channels. "
            f"Got {len(channels)} channels: {channels!r}."
        )

    if source_sink_matrix is not None and implementation != "source_sink_n":
        raise ValueError(
            "source_sink_matrix is only used when implementation='source_sink_n'."
        )
    if (
        sink_channels is not None or neutral_channels is not None
    ) and implementation != "source_sink_n":
        raise ValueError(
            "sink_channels and neutral_channels are only used when "
            "implementation='source_sink_n'."
        )
    if implementation == "source_sink_n":
        if source_sink_matrix is not None and (
            sink_channels is not None or neutral_channels is not None
        ):
            raise ValueError(
                "Use either source_sink_matrix or sink_channels/neutral_channels, "
                "not both at the same time."
            )
        if source_sink_matrix is None and (
            sink_channels is not None or neutral_channels is not None
        ):
            source_sink_matrix = build_source_sink_matrix(
                channels,
                sink_channels=sink_channels,
                neutral_channels=neutral_channels,
            )
        elif source_sink_matrix is None:
            source_sink_matrix = default_source_sink_matrix(len(channels))
        source_sink_matrix = validate_source_sink_matrix(
            source_sink_matrix,
            n_channels=len(channels),
        )

    if not bool(preprocess_alpha_inputs):
        _print_verbose(
            verbose,
            (
                "Note: preprocess_alpha_inputs is retained for backward "
                "compatibility but is not used to switch preprocessing inside "
                f"implementation='{implementation}'."
            ),
        )

    _print_verbose(
        verbose,
        (
            f"Starting PICASSO unmixing with implementation='{implementation}', "
            f"method='{method}', alpha_mode='{alpha_mode}', channels={channels}."
        ),
    )

    transformed_selected = np.empty_like(selected_stack, dtype=np.float32)
    matrix = None
    matrices_by_t = None
    matrix_details = None
    matrix_details_by_t = None

    if alpha_mode == "reference_t":
        if not 0 <= int(alpha_reference_t) < time_count:
            raise ValueError(
                f"alpha_reference_t must be between 0 and {time_count - 1}. "
                f"Got {alpha_reference_t!r}."
            )
        reference_channel_volumes = np.moveaxis(selected_stack[int(alpha_reference_t)], 1, 0)

        if implementation in {"matlab_3c", "matlab_n"}:
            reference_unmixed, matrix_details = run_picasso_matlab_like(
                reference_channel_volumes,
                background_percentile=float(background_percentile),
                max_iter=int(max_iter),
                step_size=float(step_size),
                qn=int(qn),
                pixel_bin_size=int(pixel_bin_size),
                alpha_clip=float(alpha_clip),
                negativity_threshold=float(negativity_threshold),
                clip_every_n_iterations=int(clip_every_n_iterations),
                require_three_channels=(implementation == "matlab_3c"),
            )
            matrix = np.asarray(matrix_details["unmixing_matrix"], dtype=np.float64)
            transformed_selected[int(alpha_reference_t)] = np.moveaxis(reference_unmixed, 0, 1)
            application_details_by_t: list[dict] = []
            for t in range(time_count):
                if t == int(alpha_reference_t):
                    application_details_by_t.append({"t": int(t), "mode": "estimated"})
                    continue
                transformed_t, apply_details = apply_matlab_incremental_sequence(
                    np.moveaxis(selected_stack[t], 1, 0),
                    background_percentile=float(background_percentile),
                    incremental_matrices=matrix_details["incremental_matrices"],
                    clip_every_n_iterations=int(clip_every_n_iterations),
                )
                transformed_selected[t] = np.moveaxis(transformed_t, 0, 1)
                application_details_by_t.append({"t": int(t), **apply_details})
            matrix_details["application_details_by_t"] = application_details_by_t
        else:
            reference_unmixed, matrix_details = run_source_sink_unmixing(
                reference_channel_volumes,
                source_sink_matrix=source_sink_matrix,
                background_percentile=float(background_percentile),
                mi_bins=int(mi_bins),
                alpha_max=float(alpha_max),
                max_alpha_voxels=max_alpha_voxels,
                random_state=int(random_state),
            )
            transformed_selected[int(alpha_reference_t)] = np.moveaxis(reference_unmixed, 0, 1)
            matrix = np.asarray(matrix_details["alpha_parameters"], dtype=np.float64)
            application_details_by_t = []
            for t in range(time_count):
                if t == int(alpha_reference_t):
                    application_details_by_t.append({"t": int(t), "mode": "estimated"})
                    continue
                transformed_t = apply_source_sink_parameters(
                    np.moveaxis(selected_stack[t], 1, 0),
                    source_sink_matrix=source_sink_matrix,
                    alpha_parameters=matrix_details["alpha_parameters"],
                    background_values=matrix_details["background_values"],
                )
                transformed_selected[t] = np.moveaxis(transformed_t, 0, 1)
                application_details_by_t.append({"t": int(t), "mode": "applied"})
            matrix_details["application_details_by_t"] = application_details_by_t
    else:
        matrix_details_by_t = []
        if implementation in {"matlab_3c", "matlab_n"}:
            matrices_by_t = np.empty((time_count, len(channels), len(channels)), dtype=np.float64)
        else:
            matrices_by_t = np.empty((time_count, len(channels), len(channels)), dtype=np.float64)

        for t in range(time_count):
            channel_volumes_t = np.moveaxis(selected_stack[t], 1, 0)
            if implementation in {"matlab_3c", "matlab_n"}:
                transformed_t, details_t = run_picasso_matlab_like(
                    channel_volumes_t,
                    background_percentile=float(background_percentile),
                    max_iter=int(max_iter),
                    step_size=float(step_size),
                    qn=int(qn),
                    pixel_bin_size=int(pixel_bin_size),
                    alpha_clip=float(alpha_clip),
                    negativity_threshold=float(negativity_threshold),
                    clip_every_n_iterations=int(clip_every_n_iterations),
                    require_three_channels=(implementation == "matlab_3c"),
                )
                matrices_by_t[t] = np.asarray(details_t["unmixing_matrix"], dtype=np.float64)
            else:
                transformed_t, details_t = run_source_sink_unmixing(
                    channel_volumes_t,
                    source_sink_matrix=source_sink_matrix,
                    background_percentile=float(background_percentile),
                    mi_bins=int(mi_bins),
                    alpha_max=float(alpha_max),
                    max_alpha_voxels=max_alpha_voxels,
                    random_state=int(random_state) + int(t),
                )
                matrices_by_t[t] = np.asarray(details_t["alpha_parameters"], dtype=np.float64)

            transformed_selected[t] = np.moveaxis(transformed_t, 0, 1)
            matrix_details_by_t.append({"t": int(t), **details_t})

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
        "implementation": implementation,
        "alpha_mode": alpha_mode,
        "alpha_reference_t": int(alpha_reference_t),
        "channels": [int(channel) for channel in channels],
        "source_sink_matrix": None
        if source_sink_matrix is None
        else np.asarray(source_sink_matrix, dtype=int).tolist(),
        "sink_channels": None
        if sink_channels is None
        else [int(channel) for channel in sink_channels],
        "neutral_channels": None
        if neutral_channels is None
        else [int(channel) for channel in neutral_channels],
        "background_percentile": float(background_percentile),
        "preprocess_alpha_inputs": bool(preprocess_alpha_inputs),
        "alpha_max": float(alpha_max),
        "mi_bins": int(mi_bins),
        "max_iter": int(max_iter),
        "tolerance": float(tolerance),
        "max_alpha_voxels": None if max_alpha_voxels is None else int(max_alpha_voxels),
        "random_state": int(random_state),
        "step_size": float(step_size),
        "qN": int(qn),
        "pixel_bin_size": int(pixel_bin_size),
        "alpha_clip": float(alpha_clip),
        "negativity_threshold": float(negativity_threshold),
        "clip_every_n_iterations": int(clip_every_n_iterations),
        "clip_negative": bool(clip_negative),
        "unmixing_matrix": None if matrix is None else matrix.tolist(),
        "unmixing_matrix_by_t": None
        if matrices_by_t is None
        else matrices_by_t.tolist(),
        "iterations_run": None
        if matrix_details is None
        else (
            None
            if matrix_details.get("iterations_run") is None
            else int(matrix_details["iterations_run"])
        ),
        "iterations_run_by_t": None
        if matrix_details_by_t is None
        else [
            None if item.get("iterations_run") is None else int(item["iterations_run"])
            for item in matrix_details_by_t
        ],
        "converged": None,
        "converged_by_t": None,
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
# %% END
