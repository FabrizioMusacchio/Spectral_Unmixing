"""
Interactive VS Code user script for fine-tuned filtering and registration of an
unmixed TZCYX stack.

This variant demonstrates the extended filtering API with:
- a dedicated filter sequence for the second channel (index 1),
- optional per-time-point filter strengths,
- standard Z-drift correction, time registration, histogram matching, and
  max-Z projection.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

from spectral_unmixing.filters import apply_filters, match_histograms_across_time, max_z_project
from spectral_unmixing.io import load_stack_with_omio, write_stack_with_omio
from spectral_unmixing.registration import correct_intra_stack_z_drift, register_stack

import omio as om
# %% INPUT AND OUTPUT PATHS
INPUT_PATH = (PROJECT_ROOT / "example_data" / "Gockel_Nieves_Rivera_2026" / "Gockel_Nieves_Rivera_2026_5D_stack.tif")
INPUT_NAME = INPUT_PATH.stem

OUTPUT_DIR = INPUT_PATH.parent / "registered_fine_filtered"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_PATH = OUTPUT_DIR / f"{INPUT_NAME}_registered_fine_filtered.tif"
# %% USER SETTINGS
REGISTRATION_CHANNEL = 0
INTRA_STACK_METHOD = "pystackreg"
TEMPORAL_REGISTRATION_METHOD = "pystackreg"
NEIGHBOR_WINDOW_SIZE = 3
PROJECTION_ZRANGE = (10, 19)
# %% FILTER SETTINGS
# Base filter settings applied to all channels unless overridden for channel 2
# (the second channel, index 1).
FILTERS = ["median", "gaussian"]
FILTERS_CHANNEL2 = ["median"]

# Per-time-point filtering: when list length matches T, each time point gets its
# own value. Otherwise only the first entry is used for all time points.
MEDIAN_SIZE = [3, 3, 3, 3, 3, 3, 3, 3, 3]
GAUSSIAN_SIGMA = [0.8, 0.8, 0.9, 1.0, 1.0, 0.9, 0.8, 0.8, 0.8]

# Optional second-channel-specific overrides.
MEDIAN_SIZE_CHANNEL2 = [5, 5, 5, 5, 5, 5, 5, 5, 5]
GAUSSIAN_SIGMA_CHANNEL2 = None

POST_FILTERS = ["median", "gaussian"]
POST_FILTERS_CHANNEL2 = ["median"]
POST_MEDIAN_SIZE = 3
POST_GAUSSIAN_SIGMA = 1.2
POST_MEDIAN_SIZE_CHANNEL2 = 5
POST_GAUSSIAN_SIGMA_CHANNEL2 = None
# %% LOAD STACK WITH OMIO
stack, metadata = load_stack_with_omio(INPUT_PATH)
print(f"Loaded stack: {stack.shape}, axes={metadata.get('axes')}")
om.open_in_napari(stack, metadata, "Unmixed raw |")
# %% CORRECT INTRA-STACK Z-DRIFT
z_corrected_stack = correct_intra_stack_z_drift(
    stack,
    registration_channel=REGISTRATION_CHANNEL,
    method=INTRA_STACK_METHOD,
    reference_mode="neighbor",
    neighbor_window_size=NEIGHBOR_WINDOW_SIZE,
    pre_median_filter=True,
    post_median_filter=False,
    median_kernel_size=3,
    verbose=True)
print(f"Z-corrected stack: {z_corrected_stack.shape}")
z_corrected_metadata = om.update_metadata_from_image(metadata, z_corrected_stack)
om.open_in_napari(z_corrected_stack, z_corrected_metadata, "Z-corrected |")
# %% REGISTER STACK ACROSS TIME
registered_stack = register_stack(
    z_corrected_stack,
    registration_channel=REGISTRATION_CHANNEL,
    method=TEMPORAL_REGISTRATION_METHOD,
    zrange=PROJECTION_ZRANGE,
    pre_median_filter=True,
    post_median_filter=True,
    median_kernel_size=5,
    verbose=True,)
print(f"Registered stack: {registered_stack.shape}")
registered_metadata = om.update_metadata_from_image(metadata, registered_stack)
om.open_in_napari(registered_stack, registered_metadata, "Registered |")
# %% MATCH HISTOGRAMS ACROSS TIME
matched_stack = match_histograms_across_time(registered_stack, reference_t=0)
print(f"Histogram matched stack: {matched_stack.shape}")
matched_metadata = om.update_metadata_from_image(metadata, matched_stack)
om.open_in_napari(matched_stack, matched_metadata, "Registered + hist matched |")
# %% FILTER REGISTERED STACK
filtered_stack = apply_filters(
    matched_stack,
    filters=FILTERS,
    filters_channel2=FILTERS_CHANNEL2,
    median_size=MEDIAN_SIZE,
    gaussian_sigma=GAUSSIAN_SIGMA,
    median_size_channel2=MEDIAN_SIZE_CHANNEL2,
    gaussian_sigma_channel2=GAUSSIAN_SIGMA_CHANNEL2,
    apply_3d=False)
print(f"Filtered stack: {filtered_stack.shape}")
filtered_metadata = om.update_metadata_from_image(metadata, filtered_stack)
om.open_in_napari(filtered_stack, filtered_metadata, "Fine filtered |")
# %% MAX-Z-PROJECT
projected_stack = max_z_project(filtered_stack, zrange=PROJECTION_ZRANGE)
print(f"Projected stack: {projected_stack.shape}")
projected_metadata = om.update_metadata_from_image(metadata, projected_stack)
om.open_in_napari(projected_stack, projected_metadata, "Projected |")
# %% FILTER PROJECTED STACK AGAIN
filtered_projected_stack = apply_filters(
    projected_stack,
    filters=POST_FILTERS,
    filters_channel2=POST_FILTERS_CHANNEL2,
    median_size=POST_MEDIAN_SIZE,
    gaussian_sigma=POST_GAUSSIAN_SIGMA,
    median_size_channel2=POST_MEDIAN_SIZE_CHANNEL2,
    gaussian_sigma_channel2=POST_GAUSSIAN_SIGMA_CHANNEL2,
    apply_3d=False)
print(f"Filtered projected stack: {filtered_projected_stack.shape}")
filtered_projected_metadata = om.update_metadata_from_image(metadata, filtered_projected_stack)
om.open_in_napari(filtered_projected_stack, filtered_projected_metadata, "Fine filtered projected |")
# %% SAVE FILTERED PROJECTED STACK WITH OMIO
saved_output = write_stack_with_omio( OUTPUT_PATH, filtered_projected_stack, metadata,)
print(saved_output)
# %% END
