"""
Generate a synthetic TZCYX microscopy stack with controlled two-channel bleed-through.

The script creates a 5D stack with shape ``(T, Z, C, Y, X) = (9, 20, 2, 128, 128)``
containing two time-varying 3D Gaussian spheres:

- channel 0: source sphere
- channel 1: true target sphere plus bleed-through from channel 0

The measured target channel is constructed as

.. math::

    I_1 = T_1 + \alpha I_0

where ``alpha`` is a configurable bleed-through coefficient.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import sys
# PATH SETUP:
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from spectral_unmixing.io import CANONICAL_AXIS_ORDER, write_stack_with_omio
# %% INPUT AND OUTPUT PATHS
OUTPUT_DIR = PROJECT_ROOT / "example_data" / "synthetic_data"
OUTPUT_PATH = OUTPUT_DIR / "synthetic_bleedthrough_T9_Z20_C2.tif"
REPORT_PATH = OUTPUT_PATH.with_suffix(OUTPUT_PATH.suffix + ".json")
# %% SYNTHETIC STACK PARAMETERS
SIZE_T = 9
SIZE_Z = 20
SIZE_C = 2
SIZE_Y = 128
SIZE_X = 128

ALPHA = 0.28
RANDOM_SEED = 42
BACKGROUND_LEVEL = 12.0
NOISE_SIGMA = 1.5

# %% SYNTHETIC STACK GENERATION FUNCTIONS
def gaussian_sphere_3d(
    z_grid: np.ndarray,
    y_grid: np.ndarray,
    x_grid: np.ndarray,
    *,
    center_z: float,
    center_y: float,
    center_x: float,
    sigma_z: float,
    sigma_y: float,
    sigma_x: float,
    amplitude: float,
) -> np.ndarray:
    """
    Return one anisotropic 3D Gaussian sphere volume.

    Parameters
    ----------
    z_grid, y_grid, x_grid : numpy.ndarray
        Broadcast-compatible coordinate grids.
    center_z, center_y, center_x : float
        Sphere center coordinates.
    sigma_z, sigma_y, sigma_x : float
        Gaussian standard deviations along each axis.
    amplitude : float
        Peak amplitude of the sphere.

    Returns
    -------
    numpy.ndarray
        ``float32`` volume with the same shape as the coordinate grids.
    """

    squared_distance = (
        ((z_grid - center_z) / sigma_z) ** 2
        + ((y_grid - center_y) / sigma_y) ** 2
        + ((x_grid - center_x) / sigma_x) ** 2
    )
    return (amplitude * np.exp(-0.5 * squared_distance)).astype(np.float32, copy=False)

def build_synthetic_stack(
    *,
    size_t: int = SIZE_T,
    size_z: int = SIZE_Z,
    size_y: int = SIZE_Y,
    size_x: int = SIZE_X,
    alpha: float = ALPHA,
    background_level: float = BACKGROUND_LEVEL,
    noise_sigma: float = NOISE_SIGMA,
    random_seed: int = RANDOM_SEED,
) -> tuple[np.ndarray, dict[str, object]]:
    """
    Build the synthetic measured stack and a metadata report.

    Parameters
    ----------
    size_t, size_z, size_y, size_x : int, optional
        Output stack sizes.
    alpha : float, optional
        Bleed-through coefficient from channel 0 into channel 1.
    background_level : float, optional
        Constant background intensity added to both measured channels.
    noise_sigma : float, optional
        Standard deviation of additive Gaussian noise.
    random_seed : int, optional
        Seed for reproducible noise generation.

    Returns
    -------
    tuple
        ``(stack, report)`` where ``stack`` is the measured ``TZCYX`` array and
        ``report`` stores the generative parameters.
    """

    rng = np.random.default_rng(int(random_seed))

    z_coords = np.arange(size_z, dtype=np.float32)[:, None, None]
    y_coords = np.arange(size_y, dtype=np.float32)[None, :, None]
    x_coords = np.arange(size_x, dtype=np.float32)[None, None, :]

    stack = np.zeros((size_t, size_z, SIZE_C, size_y, size_x), dtype=np.float32)
    true_source = np.zeros((size_t, size_z, size_y, size_x), dtype=np.float32)
    true_target = np.zeros((size_t, size_z, size_y, size_x), dtype=np.float32)

    for t in range(size_t):
        phase = 2.0 * np.pi * float(t) / float(max(size_t, 1))

        source_center_z = 9.5 + 1.5 * np.sin(phase)
        source_center_y = 44.0 + 8.0 * np.sin(phase)
        source_center_x = 42.0 + 10.0 * np.cos(phase)
        source_amplitude = 120.0 + 10.0 * np.cos(phase)

        target_center_z = 11.0 + 1.0 * np.cos(phase + 0.4)
        target_center_y = 82.0 + 6.0 * np.cos(phase + 0.7)
        target_center_x = 86.0 + 7.0 * np.sin(phase + 0.5)
        target_amplitude = 92.0 + 8.0 * np.sin(phase + 0.3)

        source_volume = gaussian_sphere_3d(
            z_coords,
            y_coords,
            x_coords,
            center_z=source_center_z,
            center_y=source_center_y,
            center_x=source_center_x,
            sigma_z=2.8,
            sigma_y=8.0,
            sigma_x=8.0,
            amplitude=source_amplitude,
        )
        target_volume = gaussian_sphere_3d(
            z_coords,
            y_coords,
            x_coords,
            center_z=target_center_z,
            center_y=target_center_y,
            center_x=target_center_x,
            sigma_z=3.2,
            sigma_y=9.0,
            sigma_x=9.0,
            amplitude=target_amplitude,
        )

        measured_source = source_volume + background_level
        measured_target = target_volume + alpha * source_volume + background_level

        if noise_sigma > 0.0:
            measured_source = measured_source + rng.normal(
                loc=0.0,
                scale=noise_sigma,
                size=measured_source.shape,
            ).astype(np.float32)
            measured_target = measured_target + rng.normal(
                loc=0.0,
                scale=noise_sigma,
                size=measured_target.shape,
            ).astype(np.float32)

        stack[t, :, 0, :, :] = np.clip(measured_source, 0.0, None)
        stack[t, :, 1, :, :] = np.clip(measured_target, 0.0, None)
        true_source[t] = source_volume
        true_target[t] = target_volume

    report = {
        "description": "Synthetic two-channel TZCYX stack with Gaussian 3D spheres and channel-0 bleed-through into channel 1.",
        "axis_order": CANONICAL_AXIS_ORDER,
        "shape": [int(v) for v in stack.shape],
        "size_t": int(size_t),
        "size_z": int(size_z),
        "size_c": int(SIZE_C),
        "size_y": int(size_y),
        "size_x": int(size_x),
        "alpha": float(alpha),
        "background_level": float(background_level),
        "noise_sigma": float(noise_sigma),
        "random_seed": int(random_seed),
        "channel_0_mean": float(np.mean(stack[:, :, 0, :, :])),
        "channel_1_mean": float(np.mean(stack[:, :, 1, :, :])),
        "true_source_max": float(np.max(true_source)),
        "true_target_max": float(np.max(true_target)),
    }
    return stack.astype(np.float32, copy=False), report

def build_omio_metadata(shape: tuple[int, int, int, int, int]) -> dict[str, object]:
    """
    Create a minimal OMIO metadata mapping for a synthetic ``TZCYX`` stack.

    Parameters
    ----------
    shape : tuple of int
        Stack shape in canonical ``TZCYX`` order.

    Returns
    -------
    dict
        Minimal metadata mapping compatible with the package OMIO writer helper.
    """

    return {
        "axes": CANONICAL_AXIS_ORDER,
        "shape": tuple(int(v) for v in shape),
        "SizeT": int(shape[0]),
        "SizeZ": int(shape[1]),
        "SizeC": int(shape[2]),
        "SizeY": int(shape[3]),
        "SizeX": int(shape[4]),
        "PhysicalSizeX": 0.2,
        "PhysicalSizeY": 0.2,
        "PhysicalSizeZ": 1.0,
        "PhysicalSizeXUnit": "um",
        "PhysicalSizeYUnit": "um",
        "PhysicalSizeZUnit": "um",
        "TimeIncrement": 1.0,
        "TimeIncrementUnit": "min",
        "Name": "synthetic_bleedthrough_T9_Z20_C2",
    }

# %% MAIN
def main() -> None:
    """Generate and save the synthetic stack plus a JSON sidecar report."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    stack, report = build_synthetic_stack()
    metadata = build_omio_metadata(tuple(int(v) for v in stack.shape))

    actual_output_path = write_stack_with_omio(OUTPUT_PATH, stack, metadata)
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"Wrote synthetic stack to: {actual_output_path}")
    print(f"Wrote synthetic report to: {REPORT_PATH}")
    print(json.dumps(report, indent=2))

# %% MAIN GUARD
if __name__ == "__main__":
    main()
# %% END
