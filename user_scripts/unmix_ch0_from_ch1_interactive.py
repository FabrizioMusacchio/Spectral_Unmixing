"""
Interactive VS Code user script for spectral bleed-through correction.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

import json
import sys
from pathlib import Path

# PATH SETUP:
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from spectral_unmixing import unmix_ch0_from_ch1
# %% INPUT AND OUTPUT PATHS
INPUT_PATH = (PROJECT_ROOT / "example_data" / "MicroSynDep_private" / "ID14135_TP0_d2.tif")
OUTPUT_DIR = INPUT_PATH.parent / "unmixed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FIXED        = OUTPUT_DIR / "ID14135_TP0_d2_unmixed_fixed_alpha.tif"
OUTPUT_REFERENCE    = OUTPUT_DIR / "ID14135_TP0_d2_unmixed_reference_t0.tif"
OUTPUT_PER_T        = OUTPUT_DIR / "ID14135_TP0_d2_unmixed_per_t.tif"
# %% FIXED ALPHA EXAMPLE
fixed_report = unmix_ch0_from_ch1(
    input_path=INPUT_PATH,
    output_path=OUTPUT_FIXED,
    alpha=0.72,
    alpha_mode="fixed",)
print(json.dumps(fixed_report, indent=2))
# %% REFERENCE-TIME-POINT ALPHA EXAMPLE
reference_report = unmix_ch0_from_ch1(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE,
    alpha_mode="reference_t",
    alpha_reference_t=0,
    signal_percentile=99.0,
    background_percentile=1.0,)
print(json.dumps(reference_report, indent=2))
# %% PER-TIME-POINT ALPHA EXAMPLE
per_t_report = unmix_ch0_from_ch1(
    input_path=INPUT_PATH,
    output_path=OUTPUT_PER_T,
    alpha_mode="per_t",
    signal_percentile=99.0,
    background_percentile=1.0,)
print(json.dumps(per_t_report, indent=2))
# %% END