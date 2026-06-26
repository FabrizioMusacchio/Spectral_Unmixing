"""
Napari viewer helpers for spectral unmixing results.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .io import _configure_omio_runtime_environment, load_stack_with_omio
# %% CONSTANTS
_VIEWER = None
DEFAULT_COLORMAPS = [
    "cyan",
    "magenta",
    "yellow",
    "green",
    "red",
    "blue",
    "bop orange",
    "bop purple",
    "gray",
    "gray_r",
]

# %% INTERNAL HELPERS
def import_napari():
    """
    Import napari after configuring a writable runtime environment.

    Returns
    -------
    module
        Imported :mod:`napari` module.
    """

    _configure_omio_runtime_environment()
    import napari  # pylint: disable=import-outside-toplevel

    return napari


def _find_layer(viewer, layer_name: str):
    """Return the napari layer with a matching name, or ``None`` if absent."""

    for layer in viewer.layers:
        if layer.name == layer_name:
            return layer
    return None


def _viewer_is_open(viewer, napari_module) -> bool:
    """
    Return ``True`` when a viewer instance still appears to be open and reusable.

    The check is intentionally defensive because napari state differs slightly
    between real GUI sessions and mocked test environments.
    """

    if viewer is None:
        return False
    if bool(getattr(viewer, "closed", False)):
        return False

    current_viewer = None
    try:
        current_viewer = napari_module.current_viewer()
    except Exception:
        current_viewer = None

    if current_viewer is viewer:
        return True
    if current_viewer is not None and current_viewer is not viewer:
        return False

    qt_window = getattr(getattr(viewer, "window", None), "_qt_window", None)
    if qt_window is not None and hasattr(qt_window, "isVisible"):
        try:
            return bool(qt_window.isVisible())
        except Exception:
            return False

    # If napari no longer reports a current viewer and no GUI visibility state
    # is available, assume that the stored viewer should not be reused.
    return False


def _get_or_create_viewer(title: str = "Spectral Unmixing Results"):
    """Reuse an existing napari viewer when possible, otherwise create one."""

    global _VIEWER

    napari = import_napari()

    current_viewer = None
    try:
        current_viewer = napari.current_viewer()
    except Exception:
        current_viewer = None

    if _viewer_is_open(current_viewer, napari):
        _VIEWER = current_viewer
        return _VIEWER

    if _viewer_is_open(_VIEWER, napari):
        return _VIEWER

    _VIEWER = napari.Viewer(title=title)
    return _VIEWER


def _metadata_scale_from_tzcyx(metadata: dict[str, Any]) -> tuple[float, float, float, float]:
    """Extract napari axis scaling for ``T``, ``Z``, ``Y``, and ``X`` from OMIO metadata."""

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
    """Create or update one image layer in the shared napari viewer."""

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


def show_all_channels_in_napari(
    stack_path: str | Path,
    *,
    layer_prefix: str = "Channels",
    colormaps: list[str] | None = None,
    title: str = "Spectral Unmixing Results",
):
    """
    Show every channel of a canonical ``TZCYX`` stack as a separate napari layer.

    Parameters
    ----------
    stack_path : str or Path
        Path to a microscopy stack readable by OMIO.
    layer_prefix : str, optional
        Prefix used when naming napari layers.
    colormaps : list of str or None, optional
        Colormaps used for successive channels. If more channels than colormaps
        are present, the list is reused cyclically.
    title : str, optional
        Window title used when a new napari viewer is created.

    Returns
    -------
    napari.Viewer
        Shared napari viewer containing one layer per channel.
    """

    stack_path = Path(stack_path)
    stack, metadata = load_stack_with_omio(stack_path)
    viewer = _get_or_create_viewer(title=title)
    scale = _metadata_scale_from_tzcyx(metadata)
    colormaps = DEFAULT_COLORMAPS if colormaps is None else list(colormaps)

    for channel_index in range(stack.shape[2]):
        channel_data = np.asarray(stack[:, :, channel_index, :, :], dtype=np.float32)
        _upsert_image_layer(
            viewer,
            f"{layer_prefix} | C{channel_index}",
            channel_data,
            scale=scale,
            colormap=colormaps[channel_index % len(colormaps)],
        )

    return viewer


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

    Parameters
    ----------
    output_path : str or Path
        Path to an unmixed microscopy stack readable by OMIO.
    source_channel : int, optional
        Source channel index to display.
    target_channel : int, optional
        Corrected target channel index to display.
    layer_prefix : str, optional
        Prefix used when naming napari layers.
    source_colormap : str, optional
        Colormap assigned to the source-channel layer.
    target_colormap : str, optional
        Colormap assigned to the target-channel layer.

    Returns
    -------
    napari.Viewer
        Shared napari viewer containing the requested layers.

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
# %% END
