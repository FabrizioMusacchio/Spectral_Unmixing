"""
Interactive VS Code user script for creating a publication-ready render from an
unmixed TZCYX microscopy stack.

Pipeline overview
-----------------
1. Load the unmixed stack with OMIO.
2. Correct intra-stack Z-slice drift.
3. Register the stack across time.
4. Optionally match time-wise histograms to t=0.
5. Compute a max-Z projection over a user-defined Z range.
6. Build a figure-ready render with background subtraction, edge-preserving
   denoising, channel-wise contrast scaling, and mild sharpening.
7. Open the intermediate and final results in napari and save the final render.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

import sys
from pathlib import Path

# PATH SETUP:
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from spectral_unmixing.filters import match_histograms_across_time, max_z_project, apply_filters
from spectral_unmixing.io import load_stack_with_omio, write_stack_with_omio
from spectral_unmixing.publication import render_for_publication
from spectral_unmixing.registration import correct_intra_stack_z_drift, register_stack

import omio as om

import matplotlib.pyplot as plt
# %% INPUT AND OUTPUT PATHS
INPUT_PATH = (PROJECT_ROOT / "example_data" / "MicroSynDep_private" / "unmixed"
    / "ID14135_TP0_d2_unmixed_fixed_alpha.tif")
OUTPUT_DIR = INPUT_PATH.parent / "publication_render"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_PATH = OUTPUT_DIR / "ID14135_TP0_d2_unmixed_fixed_alpha_publication_render.tif"
# %% USER SETTINGS
REGISTRATION_CHANNEL = 0
INTRA_STACK_METHOD = "pystackreg"  # "phase_cross_correlation" or "pystackreg"
TEMPORAL_REGISTRATION_METHOD = "pystackreg"  # "phase_cross_correlation" or "pystackreg"
NEIGHBOR_WINDOW_SIZE = 3  # 3 -> z-1, z, z+1; 5 -> z-2 ... z+2
PROJECTION_ZRANGE =  (7, 13)  # None or e.g. (6, 14) to use only a focused subset of slices
MATCH_HISTOGRAMS = True
# %% PUBLICATION RENDER SETTINGS
# Channel 0 is typically dendrite/spines and channel 1 microglia.
BACKGROUND_METHOD = "gaussian"  # "gaussian", "white_tophat", or "none"
BACKGROUND_GAUSSIAN_SIGMA = (12.0, 8.0)
BACKGROUND_WHITE_TOPHAT_RADIUS = (9, 7)

DENOISE_METHOD = "bilateral"  # "bilateral", "median", or "none"
BILATERAL_SIGMA_COLOR = (0.06, 0.05)
BILATERAL_SIGMA_SPATIAL = (1.5, 1.5)

MEDIAN_SIZE = (3, 3)

APPLY_UNSHARP_MASK = True
UNSHARP_RADIUS = (0.8, 0.6)
UNSHARP_AMOUNT = (0.45, 0.20)

LOWER_PERCENTILE = (2.0, 2.0)
UPPER_PERCENTILE = (99.6, 99.3)
GAMMA = (0.78, 0.92)

PRE_PROJECTION_MEDIAN_FILTER = False
POST_PROJECTION_MEDIAN_FILTER = False

SHOW_INTERMEDIATE_RESULTS_IN_NAPARI = False
# %% LOAD STACK WITH OMIO
stack, metadata = load_stack_with_omio(INPUT_PATH)
print(f"Loaded stack: {stack.shape}, axes={metadata.get('axes')}")
if SHOW_INTERMEDIATE_RESULTS_IN_NAPARI:
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
    verbose=True,)
print(f"Z-corrected stack: {z_corrected_stack.shape}")
z_corrected_metadata = om.update_metadata_from_image(metadata, z_corrected_stack)
if SHOW_INTERMEDIATE_RESULTS_IN_NAPARI:
    om.open_in_napari(z_corrected_stack, z_corrected_metadata, "Z-corrected |")
# %% REGISTER STACK ACROSS TIME
if z_corrected_stack.shape[0] > 1:
    registered_stack = register_stack(
        z_corrected_stack,
        registration_channel=REGISTRATION_CHANNEL,
        method=TEMPORAL_REGISTRATION_METHOD,
        zrange=PROJECTION_ZRANGE,
        pre_median_filter=True,
        post_median_filter=True,
        median_kernel_size=5,
        verbose=True,)
else:
    registered_stack = z_corrected_stack.copy()
    print("Skipped temporal registration because T <= 1.")

print(f"Registered stack: {registered_stack.shape}")
registered_metadata = om.update_metadata_from_image(metadata, registered_stack)
if SHOW_INTERMEDIATE_RESULTS_IN_NAPARI:
    om.open_in_napari(registered_stack, registered_metadata, "Registered |")
# %% MATCH HISTOGRAMS ACROSS TIME
if MATCH_HISTOGRAMS and registered_stack.shape[0] > 1:
    matched_stack = match_histograms_across_time(registered_stack, reference_t=0)
    print(f"Histogram-matched stack: {matched_stack.shape}")
else:
    matched_stack = registered_stack.copy()
    print("Skipped histogram matching because T <= 1 or MATCH_HISTOGRAMS is False.")

matched_metadata = om.update_metadata_from_image(metadata, matched_stack)
if SHOW_INTERMEDIATE_RESULTS_IN_NAPARI:
    om.open_in_napari(matched_stack, matched_metadata, "Registered + hist matched |")
# %% FILTER REGISTERED STACK
if PRE_PROJECTION_MEDIAN_FILTER:
    filtered_stack = apply_filters(
        matched_stack,
        filters=["median"],
        median_size=3,
        gaussian_sigma=1.0,
        apply_3d=False)
    print(f"Filtered stack: {filtered_stack.shape}")
    filtered_metadata = om.update_metadata_from_image(metadata, filtered_stack)
    if SHOW_INTERMEDIATE_RESULTS_IN_NAPARI:
        om.open_in_napari(filtered_stack, filtered_metadata, "Filtered |")
else:
    filtered_stack = matched_stack.copy()
    filtered_metadata = om.update_metadata_from_image(metadata, filtered_stack)
# %% MAX-Z-PROJECT
projected_stack = max_z_project(filtered_stack, zrange=PROJECTION_ZRANGE)
print(f"Projected stack: {projected_stack.shape}")
projected_metadata = om.update_metadata_from_image(metadata, projected_stack)
if SHOW_INTERMEDIATE_RESULTS_IN_NAPARI:
    om.open_in_napari(projected_stack, projected_metadata, "Projected raw |")
# %% FILTER PROJECTED STACK AGAIN
if POST_PROJECTION_MEDIAN_FILTER:
    filtered_projected_stack = apply_filters(
        projected_stack,
        filters=["median"],
        median_size=3,
        gaussian_sigma=0.5,
        apply_3d=False)
    print(f"Filtered projected stack: {filtered_projected_stack.shape}")
    filtered_projected_metadata = om.update_metadata_from_image(metadata, filtered_projected_stack)
    if SHOW_INTERMEDIATE_RESULTS_IN_NAPARI:
        om.open_in_napari(filtered_projected_stack, filtered_projected_metadata, "Filtered Projected |")
else:
    filtered_projected_stack = projected_stack.copy()
    filtered_projected_metadata = om.update_metadata_from_image(metadata, filtered_projected_stack)
# %% BUILD PUBLICATION RENDER
publication_render_stack = render_for_publication(
    filtered_projected_stack,
    background_method=          BACKGROUND_METHOD,
    gaussian_sigma=             BACKGROUND_GAUSSIAN_SIGMA,
    white_tophat_radius=        BACKGROUND_WHITE_TOPHAT_RADIUS,
    denoise_method=             DENOISE_METHOD,
    bilateral_sigma_color=      BILATERAL_SIGMA_COLOR,
    bilateral_sigma_spatial=    BILATERAL_SIGMA_SPATIAL,
    median_size=                MEDIAN_SIZE,
    apply_unsharp_mask=         APPLY_UNSHARP_MASK,
    unsharp_radius=             UNSHARP_RADIUS,
    unsharp_amount=             UNSHARP_AMOUNT,
    lower_percentile=           LOWER_PERCENTILE,
    upper_percentile=           UPPER_PERCENTILE,
    gamma=                      GAMMA)
print(f"Publication render stack: {publication_render_stack.shape}")
publication_render_metadata = om.update_metadata_from_image(metadata, publication_render_stack)
if SHOW_INTERMEDIATE_RESULTS_IN_NAPARI:
    om.open_in_napari(publication_render_stack, publication_render_metadata,"Publication render |",)
# %% SAVE PUBLICATION RENDER WITH OMIO
saved_output = write_stack_with_omio(OUTPUT_PATH, publication_render_stack, publication_render_metadata)
print(saved_output)

# also save first T-slice as a single image for figure panels as png using matplotlib:
first_t_slice = publication_render_stack[0]
plt.imshow(first_t_slice[0,0,:,:], cmap="gray", vmin=0, vmax=1)
plt.axis("off")
png_output_path = OUTPUT_DIR / "ID14135_TP0_d2_unmixed_fixed_alpha_publication_render_first_t_slice_channel_0.png"
plt.savefig(png_output_path, dpi=300, bbox_inches="tight", pad_inches=0)
print(f"Saved first T-slice as PNG: {png_output_path}")
plt.close()
plt.imshow(first_t_slice[0,1,:,:], cmap="gray", vmin=0, vmax=1)
plt.axis("off")
png_output_path = OUTPUT_DIR / "ID14135_TP0_d2_unmixed_fixed_alpha_publication_render_first_t_slice_channel_1.png"
plt.savefig(png_output_path, dpi=300, bbox_inches="tight", pad_inches=0)
print(f"Saved first T-slice as PNG: {png_output_path}")
plt.close()
# %% END
