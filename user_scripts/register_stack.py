"""
Interactive VS Code user script for registering an unmixed TZCYX stack across time.

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

from spectral_unmixing.io import load_stack_with_omio, write_stack_with_omio
from spectral_unmixing.registration import register_stack
from spectral_unmixing.viewer import show_unmixed_channels_in_napari
from spectral_unmixing.filters import apply_filters, max_z_project

import omio as om
# %% INPUT AND OUTPUT PATHS
INPUT_PATH = (PROJECT_ROOT / "example_data" / "MicroSynDep_private" / "unmixed"
    / "ID14135_TP0_d2_unmixed_fixed_alpha.tif")
OUTPUT_DIR = INPUT_PATH.parent / "registered"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_PATH = OUTPUT_DIR / "ID14135_TP0_d2_unmixed_fixed_alpha_registered.tif"
# %% LOAD STACK WITH OMIO
stack, metadata = load_stack_with_omio(INPUT_PATH)
print(f"Loaded stack: {stack.shape}, axes={metadata.get('axes')}")
# %% REGISTER STACK
registered_stack = register_stack(
    stack,
    registration_channel=0,
    method="pystackreg", # phase_cross_correlation or pystackreg
    zrange=None,
    pre_median_filter=True,
    post_median_filter=True,
    median_kernel_size=5)
print(f"Registered stack: {registered_stack.shape}")

projected_metadata = om.update_metadata_from_image(metadata, registered_stack)

om.open_in_napari(registered_stack, projected_metadata, "Registered |")

om.open_in_napari(stack, metadata, "Unregistered |")

# projected_unreg_stack = max_z_project(stack)
# unreg_metadata = om.update_metadata_from_image(metadata, projected_unreg_stack)
# om.open_in_napari(projected_unreg_stack, unreg_metadata, "Unregistered |")
# %% FILTER REGISTERED STACK
filtered_stack = apply_filters(
    registered_stack,
    filters=["median", "gaussian"],
    median_size=3,
    gaussian_sigma=1.0,
    apply_3d=False,)
print(f"Filtered stack: {filtered_stack.shape}")

om.open_in_napari(filtered_stack, projected_metadata, "Filtered |")
# %% MAX-Z-PROJECT
projected_stack = max_z_project(filtered_stack)
print(f"Projected stack: {projected_stack.shape}")

projected_metadata = om.update_metadata_from_image(metadata, projected_stack)

om.open_in_napari(projected_stack, projected_metadata, "Projected |")
# %% FILTER PROJECTED STACK AGAIN
filtered_projected_stack = apply_filters(
    projected_stack,
    filters=["median", "gaussian"],
    median_size=3,
    gaussian_sigma=1.5,
    apply_3d=False,)
print(f"Filtered projected stack: {filtered_projected_stack.shape}")

om.open_in_napari(filtered_projected_stack, projected_metadata, "Filtered Projected |")
# %% SAVE REGISTERED STACK WITH OMIO
saved_output = write_stack_with_omio(OUTPUT_PATH, registered_stack, projected_metadata)
print(saved_output)
# %% END