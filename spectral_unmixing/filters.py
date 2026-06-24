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


def _apply_single_filter_tzcyx(
    stack: np.ndarray,
    *,
    filter_name: str,
    median_size: int,
    gaussian_sigma: float,
    apply_3d: bool,
) -> np.ndarray:
    filtered = np.empty_like(stack, dtype=np.float32)
    time_count, z_count, channel_count = stack.shape[:3]

    for t in range(time_count):
        for c in range(channel_count):
            volume = np.asarray(stack[t, :, c, :, :], dtype=np.float32)
            if apply_3d:
                if filter_name == "median":
                    filtered[t, :, c, :, :] = median_filter(
                        volume, size=(median_size, median_size, median_size)
                    )
                else:
                    filtered[t, :, c, :, :] = gaussian_filter(
                        volume, sigma=(gaussian_sigma, gaussian_sigma, gaussian_sigma)
                    )
            else:
                for z in range(z_count):
                    plane = volume[z, :, :]
                    if filter_name == "median":
                        filtered[t, z, c, :, :] = median_filter(
                            plane, size=(median_size, median_size)
                        )
                    else:
                        filtered[t, z, c, :, :] = gaussian_filter(
                            plane, sigma=(gaussian_sigma, gaussian_sigma)
                        )
    return filtered


def apply_filters(
    stack,
    filters: str | Sequence[str],
    *,
    median_size: int = 3,
    gaussian_sigma: float = 1.0,
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
    median_size : int, optional
        Kernel size for median filtering.
    gaussian_sigma : float, optional
        Sigma for Gaussian filtering.
    apply_3d : bool, optional
        If True, apply filters in 3D over ``ZYX`` for each ``T`` and ``C`` volume.
        If False, apply them plane-wise in ``YX`` for each available ``T`` and ``Z``.

    Returns
    -------
    np.ndarray
        Filtered stack with the same shape as the input.
    """

    if median_size < 1:
        raise ValueError(f"median_size must be >= 1. Got {median_size!r}.")
    if gaussian_sigma <= 0:
        raise ValueError(f"gaussian_sigma must be > 0. Got {gaussian_sigma!r}.")

    filter_sequence = _normalize_filter_sequence(filters)
    original_stack = np.asarray(stack)
    original_ndim = original_stack.ndim
    working_stack = _ensure_tzcyx_stack(original_stack).astype(np.float32, copy=True)

    for filter_name in filter_sequence:
        working_stack = _apply_single_filter_tzcyx(
            working_stack,
            filter_name=filter_name,
            median_size=median_size,
            gaussian_sigma=gaussian_sigma,
            apply_3d=apply_3d,
        )

    return _restore_original_shape(working_stack, original_ndim)


def max_z_project(stack) -> np.ndarray:
    """
    Compute a maximum-intensity projection over the Z axis while preserving ``T`` and ``C``.

    The returned stack stays in canonical ``TZCYX`` order with a singleton Z dimension.
    """

    stack_tzcyx = _ensure_tzcyx_stack(stack)
    projected = np.max(stack_tzcyx, axis=1, keepdims=True)
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
