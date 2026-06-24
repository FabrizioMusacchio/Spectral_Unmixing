"""
Interactive VS Code user script for spectral bleed-through correction.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from spectral_unmixing import (
    report_path_from_output_path,
    show_unmixed_channels_in_napari,
    unmix)
# %% INPUT AND OUTPUT PATHS
INPUT_PATH = PROJECT_ROOT / "example_data" / "MicroSynDep_private" / "ID14135_TP0_d2.tif"
OUTPUT_DIR = INPUT_PATH.parent / "unmixed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FIXED = OUTPUT_DIR / "ID14135_TP0_d2_unmixed_fixed_alpha.tif"
OUTPUT_REFERENCE = OUTPUT_DIR / "ID14135_TP0_d2_unmixed_reference_t0.tif"
OUTPUT_PER_T = OUTPUT_DIR / "ID14135_TP0_d2_unmixed_per_t.tif"
# %% FIXED ALPHA EXAMPLE
fixed_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_FIXED,
    alpha=0.62,
    alpha_mode="fixed")

show_unmixed_channels_in_napari(
    fixed_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Fixed alpha",
    source_colormap="cyan",
    target_colormap="yellow")
# %% REFERENCE-TIME-POINT ALPHA EXAMPLE
reference_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE,
    alpha_mode="reference_t",
    alpha_reference_t=0,
    signal_percentile=99.0,
    background_percentile=1.0,
)
print(reference_output)
print(report_path_from_output_path(reference_output).read_text(encoding="utf-8"))
show_unmixed_channels_in_napari(
    reference_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Reference t0",
)
# %% PER-TIME-POINT ALPHA EXAMPLE
per_t_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_PER_T,
    alpha_mode="per_t",
    signal_percentile=99.0,
    background_percentile=1.0,
)
print(per_t_output)
print(report_path_from_output_path(per_t_output).read_text(encoding="utf-8"))
show_unmixed_channels_in_napari(
    per_t_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Per t",
)
