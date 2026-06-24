"""
Registration helpers for TZCYX microscopy stacks.

Author: Fabrizio Musacchio
Date: June 2026
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy.ndimage import median_filter, shift as ndi_shift
from skimage.registration import phase_cross_correlation

from .io import CANONICAL_AXIS_ORDER

SUPPORTED_REGISTRATION_METHODS = {"phase_cross_correlation", "pystackreg"}


def _normalize_registration_method(method: str) -> str:
    normalized = str(method).strip().lower()
    if normalized not in SUPPORTED_REGISTRATION_METHODS:
        raise ValueError(
            f"Unsupported registration method {method!r}. "
            f"Supported methods: {sorted(SUPPORTED_REGISTRATION_METHODS)}."
        )
    return normalized


def _normalize_zrange(zrange: tuple[int, int] | Sequence[int] | None, z_count: int) -> tuple[int, int]:
    if zrange is None:
        return 0, z_count

    if len(zrange) != 2:
        raise ValueError("zrange must be None or a tuple/list with exactly two integers.")

    start = int(zrange[0])
    stop = int(zrange[1])
    if not 0 <= start < stop <= z_count:
        raise ValueError(
            f"zrange must satisfy 0 <= start < stop <= {z_count}. Got {(start, stop)!r}."
        )
    return start, stop


def _ensure_tzcyx_stack(stack) -> np.ndarray:
    stack = np.asarray(stack)
    if stack.ndim != 5:
        raise ValueError(
            f"Expected a {CANONICAL_AXIS_ORDER} stack with 5 dimensions. Got shape {stack.shape!r}."
        )
    return stack


def _apply_median_to_zyx(volume_zyx: np.ndarray, kernel_size: int) -> np.ndarray:
    filtered = np.empty_like(volume_zyx, dtype=np.float32)
    for z in range(volume_zyx.shape[0]):
        filtered[z, :, :] = median_filter(volume_zyx[z, :, :], size=(kernel_size, kernel_size))
    return filtered


def _build_registration_projections(
    stack: np.ndarray,
    *,
    registration_channel: int,
    zrange: tuple[int, int] | None,
    pre_median_filter: bool,
    post_median_filter: bool,
    median_kernel_size: int,
) -> np.ndarray:
    z_start, z_stop = _normalize_zrange(zrange, stack.shape[1])
    channel_stack = np.asarray(stack[:, z_start:z_stop, registration_channel, :, :], dtype=np.float32)
    working = channel_stack.copy()

    if pre_median_filter:
        for t in range(working.shape[0]):
            working[t, :, :, :] = _apply_median_to_zyx(working[t, :, :, :], median_kernel_size)

    projections = np.max(working, axis=1)

    if post_median_filter:
        for t in range(projections.shape[0]):
            projections[t, :, :] = median_filter(
                projections[t, :, :], size=(median_kernel_size, median_kernel_size)
            )

    return projections


def _phase_cross_correlation_shift(reference_projection: np.ndarray, moving_projection: np.ndarray) -> np.ndarray:
    shift_2d, _, _ = phase_cross_correlation(reference_projection, moving_projection)
    return np.asarray(shift_2d, dtype=np.float32)


def _pystackreg_shift(reference_projection: np.ndarray, moving_projection: np.ndarray) -> np.ndarray:
    from pystackreg import StackReg  # pylint: disable=import-outside-toplevel

    sr = StackReg(StackReg.TRANSLATION)
    tmat = sr.register(reference_projection.astype(np.float32), moving_projection.astype(np.float32))
    shift_yx = np.asarray([-tmat[1, 2], -tmat[0, 2]], dtype=np.float32)
    return shift_yx


def _apply_translation_to_tzyx(stack_tzyx: np.ndarray, shift_yx: np.ndarray) -> np.ndarray:
    shifted = np.empty_like(stack_tzyx, dtype=np.float32)
    for z in range(stack_tzyx.shape[0]):
        for c in range(stack_tzyx.shape[1]):
            shifted[z, c, :, :] = ndi_shift(
                np.asarray(stack_tzyx[z, c, :, :], dtype=np.float32),
                shift=tuple(float(v) for v in shift_yx),
                order=1,
                mode="constant",
                cval=0.0,
                prefilter=True,
            )
    return shifted


def _print_verbose(verbose: bool, message: str) -> None:
    if verbose:
        print(message)


def register_stack(
    stack,
    *,
    registration_channel: int,
    method: str = "phase_cross_correlation",
    zrange: tuple[int, int] | Sequence[int] | None = None,
    pre_median_filter: bool = False,
    post_median_filter: bool = False,
    median_kernel_size: int = 3,
    verbose: bool = True,
) -> np.ndarray:
    """
    Register a TZCYX stack across time using shifts estimated from Z projections.

    Parameters
    ----------
    stack : array-like
        Input stack in canonical ``TZCYX`` order.
    registration_channel : int
        Channel used to compute the time-wise registration shifts.
    method : {"phase_cross_correlation", "pystackreg"}, optional
        Backend used for shift estimation.
    zrange : tuple[int, int] or None, optional
        Optional half-open Z range ``(start, stop)`` used for the registration projection.
    pre_median_filter : bool, optional
        If True, apply a slice-wise median filter to the selected registration volume
        before max-Z projection. This affects only shift estimation, not the stack
        that is transformed.
    post_median_filter : bool, optional
        If True, apply a 2D median filter to each projection after max-Z projection.
        This affects only shift estimation, not the stack that is transformed.
    median_kernel_size : int, optional
        Median filter kernel size used by the optional pre/post filters.
    verbose : bool, optional
        If True, print the estimated shifts line-wise for each time point.

    Returns
    -------
    np.ndarray
        Registered stack with the same ``TZCYX`` shape as the input.
    """

    stack = _ensure_tzcyx_stack(stack).astype(np.float32, copy=True)
    method = _normalize_registration_method(method)

    if stack.shape[0] <= 1:
        raise ValueError("Registration requires T > 1.")
    if not 0 <= int(registration_channel) < stack.shape[2]:
        raise ValueError(
            f"registration_channel must be between 0 and {stack.shape[2] - 1}. "
            f"Got {registration_channel!r}."
        )
    if median_kernel_size < 1:
        raise ValueError(
            f"median_kernel_size must be >= 1. Got {median_kernel_size!r}."
        )

    projections = _build_registration_projections(
        stack,
        registration_channel=int(registration_channel),
        zrange=zrange,
        pre_median_filter=pre_median_filter,
        post_median_filter=post_median_filter,
        median_kernel_size=int(median_kernel_size),
    )
    reference_projection = projections[0, :, :]
    registered = stack.copy()

    _print_verbose(
        verbose,
        (
            f"Registering stack with method='{method}', registration_channel="
            f"{int(registration_channel)}, reference_t=0"
        ),
    )
    _print_verbose(verbose, "t=0 shift_y=0.000 shift_x=0.000")

    for t in range(1, stack.shape[0]):
        moving_projection = projections[t, :, :]
        if method == "phase_cross_correlation":
            shift_yx = _phase_cross_correlation_shift(reference_projection, moving_projection)
        else:
            shift_yx = _pystackreg_shift(reference_projection, moving_projection)

        _print_verbose(
            verbose,
            f"t={t} shift_y={float(shift_yx[0]):.3f} shift_x={float(shift_yx[1]):.3f}",
        )

        registered[t, :, :, :, :] = _apply_translation_to_tzyx(
            stack[t, :, :, :, :], shift_yx
        )

    return registered


__all__ = [
    "SUPPORTED_REGISTRATION_METHODS",
    "register_stack",
]
