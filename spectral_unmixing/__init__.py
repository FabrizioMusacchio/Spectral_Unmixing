"""
Spectral unmixing utilities for microscopy stacks.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from .estimation import (
    estimate_alpha_from_volume,
    estimate_picasso_unmixing_matrix_from_volume,
    make_alpha_mask,
    mutual_information_1d,
    prepare_source_target_for_alpha,
)
from .filters import apply_filters, match_histograms_across_time, max_z_project
from .registration import correct_intra_stack_z_drift, register_stack
from .unmixing import report_path_from_output_path, unmix, unmix_ch0_from_ch1, unmix_picasso
from .viewer import show_unmixed_channels_in_napari

__all__ = [
    "apply_filters",
    "correct_intra_stack_z_drift",
    "estimate_alpha_from_volume",
    "estimate_picasso_unmixing_matrix_from_volume",
    "make_alpha_mask",
    "match_histograms_across_time",
    "max_z_project",
    "mutual_information_1d",
    "prepare_source_target_for_alpha",
    "register_stack",
    "report_path_from_output_path",
    "show_unmixed_channels_in_napari",
    "unmix",
    "unmix_ch0_from_ch1",
    "unmix_picasso",
]
# %% END
