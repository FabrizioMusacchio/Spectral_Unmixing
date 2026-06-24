"""
Napari viewer helpers for spectral unmixing results.

Author: Fabrizio Musacchio
Date: June 2026
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .io import _configure_omio_runtime_environment, load_stack_with_omio

_VIEWER = None


def import_napari():
    """Import napari after configuring a writable runtime environment."""

    _configure_omio_runtime_environment()
    import napari  # pylint: disable=import-outside-toplevel

    return napari


def _find_layer(viewer, layer_name: str):
    for layer in viewer.layers:
        if layer.name == layer_name:
            return layer
    return None


def _get_or_create_viewer(title: str = "Spectral Unmixing Results"):
    global _VIEWER

    napari = import_napari()
    if _VIEWER is not None:
        try:
            _ = len(_VIEWER.layers)
            return _VIEWER
        except Exception:
            _VIEWER = None

    current_viewer = None
    try:
        current_viewer = napari.current_viewer()
    except Exception:
        current_viewer = None

    if current_viewer is not None:
        _VIEWER = current_viewer
        return _VIEWER

    _VIEWER = napari.Viewer(title=title)
    return _VIEWER


def _metadata_scale_from_tzcyx(metadata: dict[str, Any]) -> tuple[float, float, float, float]:
    return (
        float(metadata.get("TimeIncrement", 1.0) or 1.0),
        float(metadata.get("PhysicalSizeZ", 1.0) or 1.0),
        float(metadata.get("PhysicalSizeY", 1.0) or 1.0),
        float(metadata.get("PhysicalSizeX", 1.0) or 1.0),
    )


def _upsert_image_layer(
    viewer,
    layer_name: str,
    data: np.ndarray,
    *,
    scale: tuple[float, float, float, float],
    colormap: str,
    blending: str = "additive",
    opacity: float = 0.8,
) -> None:
    layer = _find_layer(viewer, layer_name)
    contrast_limits = (
        float(np.min(data)),
        float(np.max(data)) if float(np.max(data)) > float(np.min(data)) else float(np.min(data)) + 1.0,
    )

    if layer is None:
        viewer.add_image(
            data,
            name=layer_name,
            scale=scale,
            colormap=colormap,
            blending=blending,
            opacity=opacity,
            contrast_limits=contrast_limits,
        )
        return

    layer.data = data
    layer.scale = scale
    layer.colormap = colormap
    layer.blending = blending
    layer.opacity = opacity
    layer.contrast_limits = contrast_limits
    layer.visible = True


def show_unmixed_channels_in_napari(
    output_path: str | Path,
    *,
    source_channel: int = 0,
    target_channel: int = 1,
    layer_prefix: str = "Unmixed",
    source_colormap: str = "green",
    target_colormap: str = "magenta",
):
    """
    Show source and corrected target channel from an unmixed stack in a shared napari viewer.

    Repeated calls reuse the same napari viewer and update layers with matching names
    instead of opening a new viewer.
    """

    output_path = Path(output_path)
    stack, metadata = load_stack_with_omio(output_path)
    viewer = _get_or_create_viewer()
    scale = _metadata_scale_from_tzcyx(metadata)

    source_data = np.asarray(stack[:, :, source_channel, :, :], dtype=np.float32)
    target_data = np.asarray(stack[:, :, target_channel, :, :], dtype=np.float32)

    source_layer_name = f"{layer_prefix} | source C{source_channel}"
    target_layer_name = f"{layer_prefix} | target C{target_channel}"

    _upsert_image_layer(
        viewer,
        source_layer_name,
        source_data,
        scale=scale,
        colormap=source_colormap,
    )
    _upsert_image_layer(
        viewer,
        target_layer_name,
        target_data,
        scale=scale,
        colormap=target_colormap,
    )
    return viewer
