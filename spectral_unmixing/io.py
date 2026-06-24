"""
OMIO-based I/O helpers for spectral unmixing workflows.

Author: Fabrizio Musacchio
Date: June 2026
"""

from __future__ import annotations

import copy
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

CANONICAL_AXIS_ORDER = "TZCYX"


def _configure_omio_runtime_environment() -> Path:
    """
    Configure a writable cache location before importing OMIO.

    OMIO 0.2.2 currently imports napari at module import time. In restricted
    environments this can fail unless a writable cache directory is available and
    numba JIT caching is disabled. The settings applied here are conservative and
    make headless batch usage reproducible.
    """

    cache_root = Path(
        os.environ.get(
            "XDG_CACHE_HOME",
            Path(tempfile.gettempdir()) / "spectral_unmixing_cache",
        )
    )
    cache_root.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))
    os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
    return cache_root


def import_omio():
    """Import OMIO after configuring a writable runtime environment."""

    _configure_omio_runtime_environment()
    import omio as om  # pylint: disable=import-outside-toplevel

    return om


def validate_tzcxy_stack(stack: np.ndarray, metadata: dict[str, Any]) -> None:
    """Validate that the stack is a 5D array with metadata axis order ``TZCYX``."""

    axes = metadata.get("axes")
    if axes is None:
        raise ValueError("OMIO metadata do not contain an 'axes' entry.")
    if axes != CANONICAL_AXIS_ORDER:
        raise ValueError(
            f"Expected OMIO metadata axes {CANONICAL_AXIS_ORDER!r}, got {axes!r}."
        )

    if np.ndim(stack) != 5:
        raise ValueError(
            f"Expected a 5D stack in {CANONICAL_AXIS_ORDER} order, got shape "
            f"{np.shape(stack)!r}."
        )

    if len(axes) != np.ndim(stack):
        raise ValueError(
            f"Metadata axes length {len(axes)} does not match stack.ndim {np.ndim(stack)}."
        )


def load_stack_with_omio(input_path: str | Path) -> tuple[np.ndarray, dict[str, Any]]:
    """Read a microscopy stack with OMIO and validate canonical axis order."""

    om = import_omio()
    image, metadata = om.imread(str(input_path), verbose=False)
    stack = np.asarray(image)
    validate_tzcxy_stack(stack, metadata)
    return stack, metadata


def update_metadata_shape(metadata: dict[str, Any], shape: tuple[int, ...]) -> dict[str, Any]:
    """Return a deep-copied metadata mapping updated for a canonical ``TZCYX`` shape."""

    if len(shape) != 5:
        raise ValueError(f"Expected a 5D shape in TZCYX order, got {shape!r}.")

    updated = copy.deepcopy(metadata)
    updated["axes"] = CANONICAL_AXIS_ORDER
    updated["shape"] = tuple(int(v) for v in shape)
    updated["SizeT"] = int(shape[0])
    updated["SizeZ"] = int(shape[1])
    updated["SizeC"] = int(shape[2])
    updated["SizeY"] = int(shape[3])
    updated["SizeX"] = int(shape[4])
    return updated


def write_stack_with_omio(
    output_path: str | Path,
    stack: np.ndarray,
    metadata: dict[str, Any],
) -> Path:
    """
    Write a corrected stack with OMIO and rename the result to the requested path.
    """

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    om = import_omio()
    metadata_to_write = update_metadata_shape(metadata, tuple(np.shape(stack)))
    metadata_to_write["original_filename"] = output_path.name

    written_paths = om.imwrite(
        str(output_path),
        np.asarray(stack),
        metadata_to_write,
        overwrite=True,
        return_fnames=True,
        verbose=False,
    )

    actual_path = Path(written_paths[0])
    if actual_path.resolve() != output_path.resolve():
        if output_path.exists():
            output_path.unlink()
        shutil.move(str(actual_path), str(output_path))

    return output_path
