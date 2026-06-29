"""
Interactive VS Code user script for registering an unmixed TZCYX stack across time.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

from spectral_unmixing.io import load_stack_with_omio, write_stack_with_omio
from spectral_unmixing.registration import correct_intra_stack_z_drift, register_stack
from spectral_unmixing.filters import (
    apply_filters,
    match_histograms_across_time,
    max_z_project)

import omio as om
# %% INPUT AND OUTPUT PATHS
INPUT_PATH = (PROJECT_ROOT / "example_data" / "Gockel_Nieves_Rivera_2026" / "Gockel_Nieves_Rivera_2026_5D_stack.tif")
INPUT_NAME = INPUT_PATH.stem

OUTPUT_DIR = INPUT_PATH.parent / "registered"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_PATH = OUTPUT_DIR / f"{INPUT_NAME}_registered.tif"
# %% LOAD STACK WITH OMIO
stack, metadata = load_stack_with_omio(INPUT_PATH)
print(f"Loaded stack: {stack.shape}, axes={metadata.get('axes')}")
# %% INSPECT STACK IN NAPARI
om.open_in_napari(stack, metadata, "Unregistered |")
# %% CORRECT INTRA-STACK Z-DRIFT
z_corrected_stack = correct_intra_stack_z_drift(
    stack,
    registration_channel=0,  # can also be 1 if channel 1 is the more stable structure
    method="pystackreg",  # "phase_cross_correlation" "pystackreg" 
    reference_mode="neighbor",  # or "full_projection"
    neighbor_window_size=3,  # 3 -> z-1, z, z+1; 5 -> z-2 ... z+2
    pre_median_filter=True,
    post_median_filter=False,
    median_kernel_size=3,
    verbose=True)

print(f"Z-corrected stack: {z_corrected_stack.shape}")
z_corrected_metadata = om.update_metadata_from_image(metadata, z_corrected_stack)
om.open_in_napari(stack, metadata, "Unregistered |")
om.open_in_napari(z_corrected_stack, z_corrected_metadata, "Z-corrected |")
# %% REGISTER STACK ACROSS TIME
registered_stack = register_stack(
    z_corrected_stack,
    registration_channel=0,
    method="pystackreg", # phase_cross_correlation or pystackreg
    zrange=None,
    pre_median_filter=True,
    post_median_filter=True,
    median_kernel_size=5)
print(f"Registered stack: {registered_stack.shape}")
registered_metadata = om.update_metadata_from_image(metadata, registered_stack)
om.open_in_napari(registered_stack, registered_metadata, "Registered |")
# %% HISTOGRAM MATCH ACROSS TIME
# Recommendation: do this after registration and before Z projection so that
# geometry is already aligned, but the full 3D time stacks are still available.
matched_registered_stack = match_histograms_across_time(registered_stack, reference_t=0)
print(f"Histogram matched stack: {matched_registered_stack.shape}")
matched_registered_metadata = om.update_metadata_from_image(metadata, matched_registered_stack)
om.open_in_napari(matched_registered_stack, matched_registered_metadata, "Registered + hist matched |")
# %% FILTER REGISTERED STACK
filtered_stack = apply_filters(
    matched_registered_stack,
    filters=["median", "gaussian"],
    median_size=3,
    gaussian_sigma=1.0,
    apply_3d=False)
print(f"Filtered stack: {filtered_stack.shape}")
filtered_metadata = om.update_metadata_from_image(metadata, filtered_stack)
om.open_in_napari(filtered_stack, filtered_metadata, "Filtered |")
# %% MAX-Z-PROJECT
zrange=(0,10) # None or (start_z, end_z) to specify a range of z-slices to project. 
            # If None, the full z-range is projected.
projected_stack = max_z_project(filtered_stack, zrange=zrange)
print(f"Projected stack: {projected_stack.shape}")
projected_metadata = om.update_metadata_from_image(metadata, projected_stack)
# temp_projected_path = OUTPUT_DIR / "ID14135_TP0_d2_unmixed_fixed_alpha_registered_histmatched_projected_tmp.tif"
# temp_projected_saved = write_stack_with_omio(temp_projected_path, projected_stack, metadata)
om.open_in_napari(projected_stack, projected_metadata, "Projected |")
# %% FILTER PROJECTED STACK AGAIN
filtered_projected_stack = apply_filters(
    projected_stack,
    filters=["median", "gaussian"],
    median_size=3,
    gaussian_sigma=1.5,
    apply_3d=False)
print(f"Filtered projected stack: {filtered_projected_stack.shape}")
filtered_projected_metadata = om.update_metadata_from_image(metadata, filtered_projected_stack)
om.open_in_napari(filtered_projected_stack, filtered_projected_metadata, "Filtered Projected |")
# %% SAVE FILTERED PROJECTED STACK WITH OMIO
filtered_projected_path = OUTPUT_DIR / f"{INPUT_NAME}_registered_histmatched_filtered_projected_z{zrange[0]}_to_{zrange[1]}.tif"
filtered_projected_saved = write_stack_with_omio(filtered_projected_path, filtered_projected_stack, metadata)
# %% END
