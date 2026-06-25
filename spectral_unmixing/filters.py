"""
Filtering and projection helpers for spectral unmixing workflows.

Author: Fabrizio Musacchio
Date: June 2026
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy.ndimage import gaussian_filter, median_filter
from skimage.exposure import match_histograms

from .io import CANONICAL_AXIS_ORDER

SUPPORTED_FILTERS = {"median", "gaussian"}
SECOND_CHANNEL_INDEX = 1


def _normalize_filter_sequence(filters: str | Sequence[str]) -> list[str]:
    if isinstance(filters, str):
        filter_sequence = [filters]
    else:
        filter_sequence = list(filters)

    if not filter_sequence:
        raise ValueError("filters must contain at least one filter name.")

    normalized = []
    for filter_name in filter_sequence:
        normalized_name = str(filter_name).strip().lower()
        if normalized_name not in SUPPORTED_FILTERS:
            raise ValueError(
                f"Unsupported filter {filter_name!r}. Supported filters: {sorted(SUPPORTED_FILTERS)}."
            )
        normalized.append(normalized_name)
    return normalized


def _normalize_optional_filter_sequence(filters: str | Sequence[str] | None) -> list[str] | None:
    if filters is None:
        return None
    return _normalize_filter_sequence(filters)


def _ensure_tzcyx_stack(stack) -> np.ndarray:
    stack = np.asarray(stack)
    if stack.ndim == 2:
        return stack[np.newaxis, np.newaxis, np.newaxis, :, :]
    if stack.ndim == 3:
        return stack[np.newaxis, :, np.newaxis, :, :]
    if stack.ndim == 5:
        return stack
    raise ValueError(
        "Expected a stack with shape YX, ZYX, or TZCYX. "
        f"Got shape {stack.shape!r}."
    )


def _restore_original_shape(filtered_stack: np.ndarray, original_ndim: int) -> np.ndarray:
    if original_ndim == 2:
        return filtered_stack[0, 0, 0, :, :]
    if original_ndim == 3:
        return filtered_stack[0, :, 0, :, :]
    return filtered_stack


def _normalize_zrange(
    zrange: tuple[int, int] | Sequence[int] | None,
    z_count: int,
) -> tuple[int, int]:
    if zrange is None:
        return 0, z_count

    if len(zrange) != 2:
        raise ValueError("zrange must be None or a tuple/list with exactly two integers.")

    start = int(zrange[0])
    stop = int(zrange[1])

    start = max(0, min(start, z_count))
    stop = max(0, min(stop, z_count))

    if stop < start:
        start, stop = stop, start

    if start == stop:
        if start >= z_count:
            start = max(0, z_count - 1)
            stop = z_count
        else:
            stop = min(z_count, start + 1)

    return start, stop


def _normalize_time_dependent_parameter(
    value,
    *,
    time_count: int,
    name: str,
    cast,
):
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        values = list(value)
        if not values:
            raise ValueError(f"{name} must not be an empty list.")
        if len(values) == time_count:
            return [cast(item) for item in values]
        fallback = cast(values[0])
        return [fallback for _ in range(time_count)]
    scalar_value = cast(value)
    return [scalar_value for _ in range(time_count)]


def _resolve_channel2_sequence(
    primary_sequence: list,
    channel2_value,
    *,
    time_count: int,
    name: str,
    cast,
) -> list:
    if channel2_value is None:
        return list(primary_sequence)
    return _normalize_time_dependent_parameter(
        channel2_value,
        time_count=time_count,
        name=name,
        cast=cast,
    )


def _apply_filter_sequence_to_volume(
    volume_zyx: np.ndarray,
    *,
    filter_sequence: Sequence[str],
    median_size: int,
    gaussian_sigma: float,
    apply_3d: bool,
) -> np.ndarray:
    working_volume = np.asarray(volume_zyx, dtype=np.float32).copy()

    for filter_name in filter_sequence:
        if apply_3d:
            if filter_name == "median":
                working_volume = median_filter(
                    working_volume,
                    size=(median_size, median_size, median_size),
                ).astype(np.float32, copy=False)
            else:
                working_volume = gaussian_filter(
                    working_volume,
                    sigma=(gaussian_sigma, gaussian_sigma, gaussian_sigma),
                ).astype(np.float32, copy=False)
        else:
            filtered = np.empty_like(working_volume, dtype=np.float32)
            for z in range(working_volume.shape[0]):
                plane = working_volume[z, :, :]
                if filter_name == "median":
                    filtered[z, :, :] = median_filter(
                        plane,
                        size=(median_size, median_size),
                    )
                else:
                    filtered[z, :, :] = gaussian_filter(
                        plane,
                        sigma=(gaussian_sigma, gaussian_sigma),
                    )
            working_volume = filtered

    return working_volume


def _apply_filter_sequences_tzcyx(
    stack: np.ndarray,
    *,
    filter_sequence: Sequence[str],
    second_channel_filter_sequence: Sequence[str],
    median_sizes: Sequence[int],
    gaussian_sigmas: Sequence[float],
    second_channel_median_sizes: Sequence[int],
    second_channel_gaussian_sigmas: Sequence[float],
    apply_3d: bool,
) -> np.ndarray:
    filtered = np.empty_like(stack, dtype=np.float32)
    time_count, z_count, channel_count = stack.shape[:3]

    for t in range(time_count):
        for c in range(channel_count):
            volume = np.asarray(stack[t, :, c, :, :], dtype=np.float32)
            if c == SECOND_CHANNEL_INDEX and channel_count > SECOND_CHANNEL_INDEX:
                current_filter_sequence = second_channel_filter_sequence
                current_median_size = int(second_channel_median_sizes[t])
                current_gaussian_sigma = float(second_channel_gaussian_sigmas[t])
            else:
                current_filter_sequence = filter_sequence
                current_median_size = int(median_sizes[t])
                current_gaussian_sigma = float(gaussian_sigmas[t])

            if current_median_size < 1:
                raise ValueError(
                    f"median_size must be >= 1. Got {current_median_size!r} at t={t}, c={c}."
                )
            if current_gaussian_sigma <= 0:
                raise ValueError(
                    "gaussian_sigma must be > 0. "
                    f"Got {current_gaussian_sigma!r} at t={t}, c={c}."
                )

            filtered[t, :, c, :, :] = _apply_filter_sequence_to_volume(
                volume,
                filter_sequence=current_filter_sequence,
                median_size=current_median_size,
                gaussian_sigma=current_gaussian_sigma,
                apply_3d=apply_3d,
            )
    return filtered


def apply_filters(
    stack,
    filters: str | Sequence[str],
    *,
    filters_channel2: str | Sequence[str] | None = None,
    median_size: int | Sequence[int] = 3,
    gaussian_sigma: float | Sequence[float] = 1.0,
    median_size_channel2: int | Sequence[int] | None = None,
    gaussian_sigma_channel2: float | Sequence[float] | None = None,
    apply_3d: bool = False,
) -> np.ndarray:
    """
    Apply one or more filters to a microscopy stack.

    Parameters
    ----------
    stack : array-like
        Input stack in canonical ``TZCYX`` order, or a simpler ``ZYX`` / ``YX`` array.
    filters : str or sequence of str
        Either a single filter name or a sequence such as ``["median", "gaussian"]``.
        Filters are applied in the order provided.
    filters_channel2 : str or sequence of str or None, optional
        Optional filter sequence applied only to the second channel (index ``1``).
        If ``None``, the same ``filters`` sequence is used for all channels.
    median_size : int or sequence of int, optional
        Median kernel size. If a sequence with length ``T`` is provided, the value
        is applied per time point. If the sequence length does not match ``T``,
        only the first entry is used for all time points.
    gaussian_sigma : float or sequence of float, optional
        Gaussian sigma. If a sequence with length ``T`` is provided, the value
        is applied per time point. If the sequence length does not match ``T``,
        only the first entry is used for all time points.
    median_size_channel2 : int or sequence of int or None, optional
        Optional median kernel size override for the second channel (index ``1``).
        If ``None``, the values from ``median_size`` are reused.
    gaussian_sigma_channel2 : float or sequence of float or None, optional
        Optional Gaussian sigma override for the second channel (index ``1``).
        If ``None``, the values from ``gaussian_sigma`` are reused.
    apply_3d : bool, optional
        If True, apply filters in 3D over ``ZYX`` for each ``T`` and ``C`` volume.
        If False, apply them plane-wise in ``YX`` for each available ``T`` and ``Z``.

    Returns
    -------
    np.ndarray
        Filtered stack with the same shape as the input.
    """

    filter_sequence = _normalize_filter_sequence(filters)
    second_channel_filter_sequence = _normalize_optional_filter_sequence(filters_channel2)
    original_stack = np.asarray(stack)
    original_ndim = original_stack.ndim
    working_stack = _ensure_tzcyx_stack(original_stack).astype(np.float32, copy=True)
    time_count = int(working_stack.shape[0])

    median_sizes = _normalize_time_dependent_parameter(
        median_size,
        time_count=time_count,
        name="median_size",
        cast=int,
    )
    gaussian_sigmas = _normalize_time_dependent_parameter(
        gaussian_sigma,
        time_count=time_count,
        name="gaussian_sigma",
        cast=float,
    )
    second_channel_filter_sequence = (
        filter_sequence if second_channel_filter_sequence is None else second_channel_filter_sequence
    )
    second_channel_median_sizes = _resolve_channel2_sequence(
        median_sizes,
        median_size_channel2,
        time_count=time_count,
        name="median_size_channel2",
        cast=int,
    )
    second_channel_gaussian_sigmas = _resolve_channel2_sequence(
        gaussian_sigmas,
        gaussian_sigma_channel2,
        time_count=time_count,
        name="gaussian_sigma_channel2",
        cast=float,
    )

    working_stack = _apply_filter_sequences_tzcyx(
        working_stack,
        filter_sequence=filter_sequence,
        second_channel_filter_sequence=second_channel_filter_sequence,
        median_sizes=median_sizes,
        gaussian_sigmas=gaussian_sigmas,
        second_channel_median_sizes=second_channel_median_sizes,
        second_channel_gaussian_sigmas=second_channel_gaussian_sigmas,
        apply_3d=apply_3d,
    )

    return _restore_original_shape(working_stack, original_ndim)


def max_z_project(
    stack,
    *,
    zrange: tuple[int, int] | Sequence[int] | None = None,
) -> np.ndarray:
    """
    Compute a maximum-intensity projection over the Z axis while preserving ``T`` and ``C``.

    Parameters
    ----------
    stack : array-like
        Input stack in canonical ``TZCYX`` order, or a simpler ``ZYX`` / ``YX`` array.
    zrange : tuple[int, int] or None, optional
        Optional half-open Z range ``(start, stop)`` used for the projection. If
        the provided bounds fall outside the stack, they are clamped to the valid
        Z extent. If ``None``, the full Z range is used.

    Returns
    -------
    np.ndarray
        The returned stack stays in canonical ``TZCYX`` order with a singleton Z dimension.
    """

    stack_tzcyx = _ensure_tzcyx_stack(stack)
    z_start, z_stop = _normalize_zrange(zrange, stack_tzcyx.shape[1])
    projected = np.max(stack_tzcyx[:, z_start:z_stop, :, :, :], axis=1, keepdims=True)
    return projected


def match_histograms_across_time(
    stack,
    *,
    reference_t: int = 0,
) -> np.ndarray:
    """
    Match each time point to a reference time point using per-channel histogram matching.

    Parameters
    ----------
    stack : array-like
        Input stack in canonical ``TZCYX`` order.
    reference_t : int, optional
        Reference time point used for histogram matching. Default is ``0``.

    Returns
    -------
    np.ndarray
        Histogram-matched stack with the same ``TZCYX`` shape as the input.
    """

    stack_tzcyx = _ensure_tzcyx_stack(stack)
    if stack_tzcyx.ndim != 5:
        raise ValueError(
            f"Expected a {CANONICAL_AXIS_ORDER} stack. Got shape {stack_tzcyx.shape!r}."
        )
    if stack_tzcyx.shape[0] <= 1:
        raise ValueError("Histogram matching across time requires T > 1.")
    if not 0 <= int(reference_t) < stack_tzcyx.shape[0]:
        raise ValueError(
            f"reference_t must be between 0 and {stack_tzcyx.shape[0] - 1}. Got {reference_t!r}."
        )

    matched = stack_tzcyx.astype(np.float32, copy=True)
    reference_t = int(reference_t)

    for c in range(stack_tzcyx.shape[2]):
        reference_volume = np.asarray(stack_tzcyx[reference_t, :, c, :, :], dtype=np.float32)
        for t in range(stack_tzcyx.shape[0]):
            if t == reference_t:
                continue
            moving_volume = np.asarray(stack_tzcyx[t, :, c, :, :], dtype=np.float32)
            matched[t, :, c, :, :] = match_histograms(
                moving_volume,
                reference_volume,
                channel_axis=None,
            ).astype(np.float32, copy=False)

    return matched


__all__ = [
    "CANONICAL_AXIS_ORDER",
    "SUPPORTED_FILTERS",
    "apply_filters",
    "match_histograms_across_time",
    "max_z_project",
]
