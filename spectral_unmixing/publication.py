"""
Publication-render helpers for microscopy stacks.

Author: Fabrizio Musacchio
Date: June 2026
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy.ndimage import gaussian_filter, median_filter
from skimage.exposure import rescale_intensity
from skimage.filters import unsharp_mask
from skimage.morphology import disk, white_tophat
from skimage.restoration import denoise_bilateral

from .io import CANONICAL_AXIS_ORDER

SUPPORTED_BACKGROUND_METHODS = {"none", "gaussian", "white_tophat"}
SUPPORTED_DENOISE_METHODS = {"none", "median", "bilateral"}


def _ensure_tzcyx_stack(stack) -> np.ndarray:
    stack = np.asarray(stack)
    if stack.ndim != 5:
        raise ValueError(
            f"Expected a {CANONICAL_AXIS_ORDER} stack with 5 dimensions. Got shape {stack.shape!r}."
        )
    return stack.astype(np.float32, copy=True)


def _normalize_choice(value: str, *, supported: set[str], name: str) -> str:
    normalized = str(value).strip().lower()
    if normalized not in supported:
        raise ValueError(
            f"Unsupported {name} {value!r}. Supported values: {sorted(supported)}."
        )
    return normalized


def _normalize_channel_parameter(
    value,
    *,
    channel_count: int,
    name: str,
    cast,
) -> list:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        values = list(value)
        if len(values) != channel_count:
            raise ValueError(
                f"{name} must be a scalar or a sequence with one value per channel. "
                f"Expected length {channel_count}, got {len(values)}."
            )
        return [cast(item) for item in values]
    return [cast(value) for _ in range(channel_count)]


def _normalize_for_filtering(image_yx: np.ndarray) -> tuple[np.ndarray, float]:
    image_yx = np.asarray(image_yx, dtype=np.float32)
    image_max = float(np.max(image_yx))
    if image_max <= 0:
        return np.zeros_like(image_yx, dtype=np.float32), 1.0
    return image_yx / image_max, image_max


def subtract_background(
    stack,
    *,
    method: str = "gaussian",
    gaussian_sigma: float | Sequence[float] = 8.0,
    white_tophat_radius: int | Sequence[int] = 9,
) -> np.ndarray:
    """
    Subtract diffuse background from a TZCYX stack channel-wise.

    Parameters
    ----------
    stack : array-like
        Input stack in canonical ``TZCYX`` order.
    method : {"none", "gaussian", "white_tophat"}, optional
        Background subtraction method.
    gaussian_sigma : float or sequence of float, optional
        Per-channel sigma for ``method="gaussian"``.
    white_tophat_radius : int or sequence of int, optional
        Per-channel structuring-element radius for ``method="white_tophat"``.

    Returns
    -------
    np.ndarray
        Background-subtracted stack with the same shape as the input.
    """

    stack_tzcyx = _ensure_tzcyx_stack(stack)
    method = _normalize_choice(
        method,
        supported=SUPPORTED_BACKGROUND_METHODS,
        name="background method",
    )

    if method == "none":
        return stack_tzcyx

    channel_count = stack_tzcyx.shape[2]
    gaussian_sigmas = _normalize_channel_parameter(
        gaussian_sigma,
        channel_count=channel_count,
        name="gaussian_sigma",
        cast=float,
    )
    white_tophat_radii = _normalize_channel_parameter(
        white_tophat_radius,
        channel_count=channel_count,
        name="white_tophat_radius",
        cast=int,
    )

    background_subtracted = np.empty_like(stack_tzcyx, dtype=np.float32)

    for t in range(stack_tzcyx.shape[0]):
        for z in range(stack_tzcyx.shape[1]):
            for c in range(channel_count):
                plane = np.asarray(stack_tzcyx[t, z, c, :, :], dtype=np.float32)
                if method == "gaussian":
                    sigma = gaussian_sigmas[c]
                    if sigma <= 0:
                        raise ValueError(f"gaussian_sigma must be > 0. Got {sigma!r}.")
                    background = gaussian_filter(plane, sigma=(sigma, sigma))
                    corrected = plane - background
                else:
                    radius = white_tophat_radii[c]
                    if radius < 1:
                        raise ValueError(
                            f"white_tophat_radius must be >= 1. Got {radius!r}."
                        )
                    corrected = white_tophat(plane, footprint=disk(radius))

                background_subtracted[t, z, c, :, :] = np.maximum(corrected, 0.0)

    return background_subtracted


def render_for_publication(
    stack,
    *,
    background_method: str = "gaussian",
    gaussian_sigma: float | Sequence[float] = (8.0, 6.0),
    white_tophat_radius: int | Sequence[int] = (9, 7),
    denoise_method: str = "bilateral",
    bilateral_sigma_color: float | Sequence[float] = (0.12, 0.10),
    bilateral_sigma_spatial: float | Sequence[float] = (2.0, 2.0),
    median_size: int | Sequence[int] = 3,
    apply_unsharp_mask: bool = True,
    unsharp_radius: float | Sequence[float] = (1.0, 0.8),
    unsharp_amount: float | Sequence[float] = (1.2, 0.9),
    lower_percentile: float | Sequence[float] = (1.0, 1.0),
    upper_percentile: float | Sequence[float] = (99.85, 99.7),
    gamma: float | Sequence[float] = (0.85, 0.95),
) -> np.ndarray:
    """
    Build a figure-ready render from a TZCYX microscopy stack.

    Notes
    -----
    This function is meant for visualization and figure rendering, not for
    quantitative analysis. It applies channel-wise background subtraction,
    optional denoising, percentile-based contrast scaling, optional unsharp
    masking, and gamma correction.
    """

    stack_tzcyx = _ensure_tzcyx_stack(stack)
    background_method = _normalize_choice(
        background_method,
        supported=SUPPORTED_BACKGROUND_METHODS,
        name="background method",
    )
    denoise_method = _normalize_choice(
        denoise_method,
        supported=SUPPORTED_DENOISE_METHODS,
        name="denoise method",
    )

    channel_count = stack_tzcyx.shape[2]
    bilateral_sigma_colors = _normalize_channel_parameter(
        bilateral_sigma_color,
        channel_count=channel_count,
        name="bilateral_sigma_color",
        cast=float,
    )
    bilateral_sigma_spatials = _normalize_channel_parameter(
        bilateral_sigma_spatial,
        channel_count=channel_count,
        name="bilateral_sigma_spatial",
        cast=float,
    )
    median_sizes = _normalize_channel_parameter(
        median_size,
        channel_count=channel_count,
        name="median_size",
        cast=int,
    )
    unsharp_radii = _normalize_channel_parameter(
        unsharp_radius,
        channel_count=channel_count,
        name="unsharp_radius",
        cast=float,
    )
    unsharp_amounts = _normalize_channel_parameter(
        unsharp_amount,
        channel_count=channel_count,
        name="unsharp_amount",
        cast=float,
    )
    lower_percentiles = _normalize_channel_parameter(
        lower_percentile,
        channel_count=channel_count,
        name="lower_percentile",
        cast=float,
    )
    upper_percentiles = _normalize_channel_parameter(
        upper_percentile,
        channel_count=channel_count,
        name="upper_percentile",
        cast=float,
    )
    gammas = _normalize_channel_parameter(
        gamma,
        channel_count=channel_count,
        name="gamma",
        cast=float,
    )

    working_stack = subtract_background(
        stack_tzcyx,
        method=background_method,
        gaussian_sigma=gaussian_sigma,
        white_tophat_radius=white_tophat_radius,
    )
    rendered = np.empty_like(working_stack, dtype=np.float32)

    for c in range(channel_count):
        lower = lower_percentiles[c]
        upper = upper_percentiles[c]
        if not 0.0 <= lower < upper <= 100.0:
            raise ValueError(
                f"Percentiles must satisfy 0 <= lower < upper <= 100. "
                f"Got lower={lower!r}, upper={upper!r} for channel {c}."
            )

        channel_data = np.asarray(working_stack[:, :, c, :, :], dtype=np.float32)
        p_low = float(np.percentile(channel_data, lower))
        p_high = float(np.percentile(channel_data, upper))
        if p_high <= p_low:
            p_high = p_low + 1e-6

        for t in range(working_stack.shape[0]):
            for z in range(working_stack.shape[1]):
                plane = np.asarray(working_stack[t, z, c, :, :], dtype=np.float32)

                if denoise_method == "median":
                    size = median_sizes[c]
                    if size < 1:
                        raise ValueError(f"median_size must be >= 1. Got {size!r}.")
                    plane = median_filter(plane, size=(size, size))
                elif denoise_method == "bilateral":
                    sigma_color = bilateral_sigma_colors[c]
                    sigma_spatial = bilateral_sigma_spatials[c]
                    if sigma_color <= 0 or sigma_spatial <= 0:
                        raise ValueError(
                            "bilateral_sigma_color and bilateral_sigma_spatial must be > 0."
                        )
                    normalized_plane, plane_scale = _normalize_for_filtering(plane)
                    plane = denoise_bilateral(
                        normalized_plane,
                        sigma_color=sigma_color,
                        sigma_spatial=sigma_spatial,
                        channel_axis=None,
                    ).astype(np.float32, copy=False) * plane_scale

                plane = rescale_intensity(
                    plane,
                    in_range=(p_low, p_high),
                    out_range=(0.0, 1.0),
                ).astype(np.float32, copy=False)

                if apply_unsharp_mask:
                    radius = unsharp_radii[c]
                    amount = unsharp_amounts[c]
                    if radius <= 0 or amount < 0:
                        raise ValueError(
                            "unsharp_radius must be > 0 and unsharp_amount must be >= 0."
                        )
                    plane = unsharp_mask(
                        plane,
                        radius=radius,
                        amount=amount,
                        preserve_range=True,
                    ).astype(np.float32, copy=False)

                gamma_value = gammas[c]
                if gamma_value <= 0:
                    raise ValueError(f"gamma must be > 0. Got {gamma_value!r}.")
                plane = np.clip(plane, 0.0, 1.0)
                plane = np.power(plane, gamma_value, dtype=np.float32)
                rendered[t, z, c, :, :] = np.clip(plane, 0.0, 1.0)

    return rendered


__all__ = [
    "SUPPORTED_BACKGROUND_METHODS",
    "SUPPORTED_DENOISE_METHODS",
    "render_for_publication",
    "subtract_background",
]
