"""
Interactive VS Code user script for PICASSO-style blind unmixing of a 5-color simulation.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

from spectral_unmixing import report_path_from_output_path, unmix_picasso
from spectral_unmixing.viewer import show_all_channels_in_napari
# %% INPUT AND OUTPUT PATHS
"""Define the example input stack and all output targets used below.

In fact, you just need to set ``INPUT_PATH`` to your own data and the rest will be
automatically generated in a subfolder of the input file's parent directory.
"""
# define the input path to the example dataset:
INPUT_PATH = (PROJECT_ROOT / "example_data" / "PICASSO_examples" / "5_color_unmixing_simulation.tif")
OUTPUT_DIR = INPUT_PATH.parent / "unmixed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# inspect the stack in Napari:
show_all_channels_in_napari(INPUT_PATH, layer_prefix="Measured 5-color simulation")
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
# define the output path for the PICASSO MATLAB-N unmixing result:
OUTPUT_PICASSO_MATLAB_N = OUTPUT_DIR / "5_color_unmixing_simulation_picasso_matlab_n_reference_t0.tif"
picasso_matlab_n_output = unmix_picasso(
    input_path=INPUT_PATH,
    output_path=OUTPUT_PICASSO_MATLAB_N,
    channels=[0, 1, 2, 3, 4],
    # method="picasso",  # default
    implementation="matlab_n",  # "matlab_3c" or "matlab_n" or "source_sink_n"
    alpha_mode="reference_t",
    alpha_reference_t=0,
    background_percentile=1.0,
    # preprocess_alpha_inputs=True,  # recorded for compatibility
    mi_bins=64,
    alpha_max=1.0,
    max_iter=50,
    tolerance=1e-4,
    max_alpha_voxels=250_000,
    step_size=0.2,
    qn=100,
    pixel_bin_size=16,
    alpha_clip=0.5,
    # negativity_threshold=0.0009,
    # clip_every_n_iterations=50,
    random_state=42,
    clip_negative=True,
    output_dtype="float32",
    verbose=True,
)
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

# define the output path for the PICASSO source-sink-N unmixing result:
OUTPUT_PICASSO_SOURCE_SINK = OUTPUT_DIR / "5_color_unmixing_simulation_picasso_source_sink_reference_t0.tif"

sink_channels = [0, 1, 2, 3, 4]
neutral_channels = []

# Example for a more selective run:
#
# sink_channels = [1, 3]
# neutral_channels = [4]

picasso_source_sink_output = unmix_picasso(
    input_path=INPUT_PATH,
    output_path=OUTPUT_PICASSO_SOURCE_SINK,
    channels=[0, 1, 2, 3, 4],
    # method="picasso",  # default
    implementation="source_sink_n",
    alpha_mode="reference_t",
    alpha_reference_t=0,
    sink_channels=sink_channels,
    neutral_channels=neutral_channels,
    # source_sink_matrix=[
    #     [1, -1, 0, 0, 0],
    #     [0, 1, 0, -1, 0],
    #     [0, 0, 1, 0, 0],
    #     [0, -1, 0, 1, 0],
    #     [0, 0, 0, 0, 1],
    # ],
    background_percentile=1.0,
    # preprocess_alpha_inputs=True,  # recorded for compatibility
    mi_bins=64,
    alpha_max=1.0,
    max_iter=50,
    tolerance=1e-4,
    max_alpha_voxels=250_000,
    random_state=42,
    clip_negative=True,
    output_dtype="float32",
    verbose=True,
)
print(picasso_source_sink_output)
print(report_path_from_output_path(picasso_source_sink_output).read_text(encoding="utf-8"))
show_all_channels_in_napari(picasso_source_sink_output, layer_prefix="PICASSO source-sink-N unmixed 5-color simulation")
# %% END
