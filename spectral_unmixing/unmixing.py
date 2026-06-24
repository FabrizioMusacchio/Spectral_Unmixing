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

from .estimation import estimate_alpha_from_volume
from .io import CANONICAL_AXIS_ORDER, load_stack_with_omio, write_stack_with_omio

ALPHA_MODES = {"fixed", "reference_t", "per_t"}


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
    if alpha_mode not in ALPHA_MODES:
        raise ValueError(
            f"alpha_mode must be one of {sorted(ALPHA_MODES)}. Got {alpha_mode!r}."
        )
    return alpha_mode


def _estimate_reference_alpha(
    stack: np.ndarray,
    alpha_reference_t: int,
    source_channel: int,
    target_channel: int,
    signal_percentile: float,
    background_percentile: float,
) -> float:
    time_count = stack.shape[0]
    if not 0 <= alpha_reference_t < time_count:
        raise ValueError(
            f"alpha_reference_t must be between 0 and {time_count - 1}. "
            f"Got {alpha_reference_t!r}."
        )

    source = stack[alpha_reference_t, :, source_channel, :, :]
    target = stack[alpha_reference_t, :, target_channel, :, :]
    return estimate_alpha_from_volume(
        source=source,
        target=target,
        signal_percentile=signal_percentile,
        background_percentile=background_percentile,
    )


def _estimate_per_t_alphas(
    stack: np.ndarray,
    source_channel: int,
    target_channel: int,
    signal_percentile: float,
    background_percentile: float,
) -> np.ndarray:
    alpha_values = np.empty(stack.shape[0], dtype=np.float32)
    for t in range(stack.shape[0]):
        alpha_values[t] = estimate_alpha_from_volume(
            source=stack[t, :, source_channel, :, :],
            target=stack[t, :, target_channel, :, :],
            signal_percentile=signal_percentile,
            background_percentile=background_percentile,
        )
    return alpha_values


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


def unmix(
    input_path,
    output_path,
    alpha=None,
    alpha_mode="fixed",
    alpha_reference_t=0,
    source_channel=0,
    target_channel=1,
    signal_percentile=99.0,
    background_percentile=1.0,
    clip_negative=True,
    output_dtype="float32",
    verbose=True,
) -> Path:
    """
    Remove bleed-through from one source channel into one target channel in a TZCYX stack.

    Parameters
    ----------
    input_path : str or Path
        Path to the input TIF stack.
    output_path : str or Path
        Path where the corrected TIF stack should be written.
    alpha : float or None
        Fixed bleed-through coefficient. Required when ``alpha_mode == "fixed"``.
        When ``alpha_mode`` is ``"reference_t"`` or ``"per_t"``, alpha is estimated
        from the data.
    alpha_mode : {"fixed", "reference_t", "per_t"}
        Determines how alpha is obtained.
    alpha_reference_t : int
        Time point used for alpha estimation when ``alpha_mode == "reference_t"``.
    source_channel : int
        Channel that bleeds into the target channel.
    target_channel : int
        Channel from which the source contribution should be removed.
    signal_percentile : float
        Percentile threshold used to define the bright source-signal mask.
    background_percentile : float
        Percentile used for rough per-channel background subtraction before alpha
        estimation.
    clip_negative : bool
        If True, clip corrected target-channel values below zero to zero.
    output_dtype : str
        Output dtype for the saved stack. ``"float32"`` is recommended.
    verbose : bool, optional
        If True, print relevant runtime information to the terminal.

    Returns
    -------
    Path
        Path to the written corrected TIFF stack.

    Warnings
    --------
    ``alpha_mode="reference_t"`` assumes that the bleed-through factor is stable
    across time.

    ``alpha_mode="per_t"`` can compensate for slow intensity changes but may also
    introduce time-dependent artifacts if the source and target channels change
    biologically.

    A fixed alpha measured from a proper single-label control recording is
    scientifically preferable.
    """

    input_path = Path(input_path)
    output_path = Path(output_path)
    report_path = report_path_from_output_path(output_path)

    if input_path.resolve() == output_path.resolve():
        raise ValueError(
            "Refusing to overwrite the input file. Please choose a different output_path."
        )

    alpha_mode = _validate_alpha_mode(alpha_mode)
    _print_verbose(verbose, f"Reading stack with OMIO: {input_path}")

    stack, metadata = load_stack_with_omio(input_path)
    channel_count = stack.shape[2]
    time_count = int(stack.shape[0])
    z_count = int(stack.shape[1])

    source_channel = _validate_channel_index(
        "source_channel", source_channel, channel_count
    )
    target_channel = _validate_channel_index(
        "target_channel", target_channel, channel_count
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

    if alpha_mode == "fixed":
        if alpha is None:
            raise ValueError("alpha must be provided when alpha_mode='fixed'.")
        alpha_scalar = float(alpha)
        alpha_values = None
        _print_verbose(
            verbose,
            (
                f"Using fixed alpha={alpha_scalar:.6f} for source_channel="
                f"{source_channel} -> target_channel={target_channel}."
            ),
        )
    elif alpha_mode == "reference_t":
        alpha_scalar = _estimate_reference_alpha(
            stack=stack,
            alpha_reference_t=int(alpha_reference_t),
            source_channel=source_channel,
            target_channel=target_channel,
            signal_percentile=signal_percentile,
            background_percentile=background_percentile,
        )
        alpha_values = None
        _print_verbose(
            verbose,
            (
                f"Estimated reference alpha={alpha_scalar:.6f} from t="
                f"{int(alpha_reference_t)} across all {z_count} z-plane(s)."
            ),
        )
    else:
        alpha_scalar = None
        alpha_values = _estimate_per_t_alphas(
            stack=stack,
            source_channel=source_channel,
            target_channel=target_channel,
            signal_percentile=signal_percentile,
            background_percentile=background_percentile,
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
        "alpha": None if alpha_scalar is None else float(alpha_scalar),
        "alpha_values": None
        if alpha_values is None
        else [float(value) for value in np.asarray(alpha_values)],
        "source_channel": int(source_channel),
        "target_channel": int(target_channel),
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
