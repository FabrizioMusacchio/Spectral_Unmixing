"""
Interactive VS Code user script for PICASSO-style blind unmixing of a 5-color simulation.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

"""Import the helpers used throughout the multi-channel PICASSO tutorial.

This script demonstrates a full blind-unmixing workflow for the file
``example_data/PICASSO_examples/5-color unmixing simulation.tif``.

The imported helpers cover:

- ``unmix_picasso(...)`` for iterative multi-channel blind unmixing.
- ``report_path_from_output_path(...)`` for reading the JSON sidecar report.
- ``load_stack_with_omio(...)`` and ``write_stack_with_omio(...)`` for a small
  preparatory reshaping step before unmixing.
- a local napari helper in this script for visualizing all five channels as
  separate layers in one shared viewer.
"""

import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from spectral_unmixing import report_path_from_output_path, unmix_picasso
from spectral_unmixing.io import load_stack_with_omio, write_stack_with_omio
from spectral_unmixing.viewer import import_napari
# %% INPUT AND OUTPUT PATHS
"""Define the example input stack and all output targets used below.

In fact, you just need to set ``INPUT_PATH`` to your own data and the rest will be 
automatically generated in a subfolder of the input file's parent directory.
"""

INPUT_PATH = (PROJECT_ROOT / "example_data" / "PICASSO_examples" / "5-color unmixing simulation.tif")
GROUND_TRUTH_PATH = ( PROJECT_ROOT / "example_data" / "PICASSO_examples" / "5-color unmixing simulation Ground-truth.tif")
OUTPUT_DIR = INPUT_PATH.parent / "unmixed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PREPARED_INPUT_PATH = OUTPUT_DIR / "5-color unmixing simulation_T_as_C.tif"
PREPARED_GROUND_TRUTH_PATH = OUTPUT_DIR / "5-color unmixing simulation Ground-truth_T_as_C.tif"
OUTPUT_PICASSO_MATLAB_N = OUTPUT_DIR / "5-color unmixing simulation_picasso_matlab_n_reference_t0.tif"
OUTPUT_PICASSO_SOURCE_SINK = OUTPUT_DIR / "5-color unmixing simulation_picasso_source_sink_reference_t0.tif"
# %% VIEWER HELPERS
"""Define small local helpers for channel visualization in one shared napari viewer.

Why a local helper is useful here:

- the existing package helper ``show_unmixed_channels_in_napari(...)`` is built
  for two-channel source/target display
- the present PICASSO tutorial needs to show all five channels as separate
  layers
- repeated execution should reuse the same napari viewer and update layers
  instead of opening a new window every time
"""

_VIEWER = None
DEFAULT_COLORMAPS = ["cyan", "yellow", "magenta", "green", "red"]

def _get_or_create_viewer(title: str = "PICASSO 5-Color Unmixing") -> object:
    """Reuse an existing napari viewer when possible, otherwise create one."""

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

def _metadata_scale_from_tzcyx(metadata: dict) -> tuple[float, float, float, float]:
    """Extract napari axis scaling for ``T``, ``Z``, ``Y``, and ``X`` from metadata."""

    return (
        float(metadata.get("TimeIncrement", 1.0) or 1.0),
        float(metadata.get("PhysicalSizeZ", 1.0) or 1.0),
        float(metadata.get("PhysicalSizeY", 1.0) or 1.0),
        float(metadata.get("PhysicalSizeX", 1.0) or 1.0),
    )

def show_all_channels_in_napari(
    stack_path: str | Path,
    *,
    layer_prefix: str,
    colormaps: list[str] | None = None,
) -> object:
    """Open a canonical ``TZCYX`` stack and show every channel as its own napari layer."""

    stack, metadata = load_stack_with_omio(stack_path)
    viewer = _get_or_create_viewer()
    scale = _metadata_scale_from_tzcyx(metadata)
    colormaps = DEFAULT_COLORMAPS if colormaps is None else list(colormaps)

    for c in range(stack.shape[2]):
        channel_data = np.asarray(stack[:, :, c, :, :], dtype=np.float32)
        layer_name = f"{layer_prefix} | C{c}"
        colormap = colormaps[c % len(colormaps)]
        contrast_limits = (
            float(np.min(channel_data)),
            float(np.max(channel_data))
            if float(np.max(channel_data)) > float(np.min(channel_data))
            else float(np.min(channel_data)) + 1.0,
        )

        existing_layer = None
        for layer in viewer.layers:
            if layer.name == layer_name:
                existing_layer = layer
                break

        if existing_layer is None:
            viewer.add_image(
                channel_data,
                name=layer_name,
                scale=scale,
                colormap=colormap,
                blending="additive",
                opacity=0.8,
                contrast_limits=contrast_limits,
            )
        else:
            existing_layer.data = channel_data
            existing_layer.scale = scale
            existing_layer.colormap = colormap
            existing_layer.blending = "additive"
            existing_layer.opacity = 0.8
            existing_layer.contrast_limits = contrast_limits
            existing_layer.visible = True

    return viewer

# %% PREPARE INPUT STACKS
"""Convert the PICASSO example from time-encoded pages to channel-encoded TZCYX data.

Important detail:

- OMIO reads the 5-color simulation as ``T=5, Z=1, C=1``.
- For ``unmix_picasso(...)`` we instead need the five measured color images on
  the channel axis, i.e. ``T=1, Z=1, C=5``.

This cell therefore:

1. loads the measured simulation stack
2. moves the five pages from ``T`` to ``C``
3. writes the converted stack back to disk
4. repeats the same conversion for the provided ground-truth stack

What can be adjusted:

- If another dataset already arrives as ``C>1``, you can skip this conversion
  logic or adapt it accordingly.
"""

def convert_time_encoded_stack_to_channel_stack(
    input_path: str | Path,
    output_path: str | Path,
) -> Path:
    """Convert a ``TZCYX`` stack with ``T>1`` and ``C=1`` into ``T=1`` and ``C=T``."""

    stack, metadata = load_stack_with_omio(input_path)
    if stack.shape[2] != 1:
        raise ValueError(
            f"Expected a single-channel stack for T-to-C conversion. Got shape {stack.shape!r}."
        )
    converted = np.moveaxis(stack, 0, 2)
    return write_stack_with_omio(output_path, converted, metadata)

prepared_input = convert_time_encoded_stack_to_channel_stack(
    INPUT_PATH, PREPARED_INPUT_PATH)
prepared_ground_truth = convert_time_encoded_stack_to_channel_stack(
    GROUND_TRUTH_PATH, PREPARED_GROUND_TRUTH_PATH)
print(f"Prepared measured input: {prepared_input}")
print(f"Prepared ground truth: {prepared_ground_truth}")
# %% INSPECT PREPARED STACKS IN NAPARI
"""Open the converted measured and ground-truth stacks in napari.

This is a sanity-check cell:

- the measured stack should now read as ``T=1, Z=1, C=5``
- the ground truth should have the same layout
- each channel is shown as its own color layer so you can verify that the
  conversion from time pages to channels worked as intended
"""

show_all_channels_in_napari(prepared_input, layer_prefix="Measured 5-color simulation")
show_all_channels_in_napari(prepared_ground_truth, layer_prefix="Ground truth 5-color simulation")
# %% PICASSO MATLAB-N EXAMPLE
"""Run the explicit N-channel generalization of the MATLAB PICASSO workflow.

Method summary:

- ``implementation="matlab_n"`` generalizes the original MATLAB 3-channel
  PICASSO routine to an arbitrary number of channels
- each iteration estimates pairwise subtraction coefficients and applies the
  same MATLAB-style incremental update logic used by the original 3-channel
  code
- ``alpha_mode="reference_t"`` means that the update sequence is estimated once
  from the chosen reference time point and then applied to the full stack

What can be adjusted:

- ``channels``:
  Explicitly choose a subset of channels if you do not want to unmix all five
  at once.
- ``implementation``:
  Keep ``"matlab_n"`` here for the five-channel simulation. The default
  ``"matlab_3c"`` is intentionally stricter and only works with exactly three
  selected channels.
- ``alpha_max``:
  Not used by the MATLAB-like implementation, but retained in the shared API
  and JSON report.
- ``mi_bins``:
  Also not used by the MATLAB-like implementation; the MATLAB port uses ``qN``
  instead for histogram quantization.
- ``max_iter``:
  Number of MATLAB-style unmixing iterations.
- ``step_size``:
  Strength of each incremental matrix update. Larger values unmix more
  aggressively but can become unstable.
- ``qN``:
  Quantization parameter for the MATLAB mutual-information calculation.
- ``pixel_bin_size``:
  Two-dimensional binning factor applied before mutual-information estimation.
- ``alpha_clip``:
  Hard bound applied to each pairwise coefficient before update.
"""

picasso_matlab_n_output = unmix_picasso(
    input_path=prepared_input,
    output_path=OUTPUT_PICASSO_MATLAB_N,
    channels=[0, 1, 2, 3, 4],
    implementation="matlab_n", # "matlab_3c" or "matlab_n" or "source_sink_n"
    alpha_mode="reference_t",
    alpha_reference_t=0,
    background_percentile=1.0,
    mi_bins=64,
    alpha_max=1.0,
    max_iter=50,
    tolerance=1e-4,
    max_alpha_voxels=250_000,
    step_size=0.2,
    qn=100,
    pixel_bin_size=16,
    alpha_clip=0.5,
    random_state=42,
    clip_negative=True,
    output_dtype="float32",
    verbose=True)
print(picasso_matlab_n_output)
print(report_path_from_output_path(picasso_matlab_n_output).read_text(encoding="utf-8"))
show_all_channels_in_napari(
    picasso_matlab_n_output,
    layer_prefix="PICASSO MATLAB-N unmixed 5-color simulation")
# %% PICASSO SOURCE-SINK-N EXAMPLE
"""Run the source-sink N-channel variant inspired by the napari PICASSO plugin.

Method summary:

- ``implementation="source_sink_n"`` treats every selected sink channel as a
  channel that can receive modeled spillover from one or more source channels
- the relation graph can be specified either manually in
  ``source_sink_matrix`` or more readably via ``sink_channels`` and
  ``neutral_channels``
- each allowed source-to-sink coefficient is estimated by minimizing mutual
  information between the source image and the corrected sink image

What can be adjusted:

- ``sink_channels``:
  Actual channel indices that should be corrected as sinks.
- ``neutral_channels``:
  Actual channel indices that should stay neutral, meaning they are neither
  corrected as sinks nor used as sources.
- ``source_sink_matrix``:
  Optional lower-level alternative if you want to specify the full relation
  matrix manually.
- ``alpha_max``:
  Upper bound for the source-to-sink coefficients.
- ``mi_bins``:
  Histogram bin count used by the mutual-information objective.
- ``max_alpha_voxels``:
  Optional voxel cap for faster coefficient estimation on larger stacks.

This mode is not the original MATLAB PICASSO algorithm. It is a more direct
source-sink formulation that is often easier to reason about when explicit
cross-talk relations are known or suspected.

For this 5-color example the ground-truth bleed-through directions are usually
not known a priori. A practical starting point is therefore:

- choose the channels you currently want to clean as ``sink_channels``
- mark channels that should stay untouched as ``neutral_channels``
- let all remaining non-neutral selected channels act as possible sources

The example below starts with the broadest non-neutral model: every selected
channel is treated as a potential sink and a potential source. For a more
targeted run, simply restrict ``sink_channels`` or add channels to
``neutral_channels``.
"""

sink_channels = [0, 1, 2, 3, 4]
neutral_channels = []

# Example for a more selective run:
#
# sink_channels = [1, 3]
# neutral_channels = [4]

picasso_source_sink_output = unmix_picasso(
    input_path=prepared_input,
    output_path=OUTPUT_PICASSO_SOURCE_SINK,
    channels=[0, 1, 2, 3, 4],
    implementation="source_sink_n",
    alpha_mode="reference_t",
    alpha_reference_t=0,
    sink_channels=sink_channels,
    neutral_channels=neutral_channels,
    background_percentile=1.0,
    mi_bins=64,
    alpha_max=1.0,
    max_iter=50,
    tolerance=1e-4,
    max_alpha_voxels=250_000,
    random_state=42,
    clip_negative=True,
    output_dtype="float32",
    verbose=True)
print(picasso_source_sink_output)
print(report_path_from_output_path(picasso_source_sink_output).read_text(encoding="utf-8"))
show_all_channels_in_napari(picasso_source_sink_output, layer_prefix="PICASSO source-sink-N unmixed 5-color simulation")
# %% END
