"""
Spectral unmixing utilities for microscopy stacks.

Author: Fabrizio Musacchio
Date: June 2026
"""

from .estimation import estimate_alpha_from_volume
from .filters import apply_filters, max_z_project
from .registration import register_stack
from .unmixing import report_path_from_output_path, unmix, unmix_ch0_from_ch1
from .viewer import show_unmixed_channels_in_napari

__all__ = [
    "apply_filters",
    "estimate_alpha_from_volume",
    "max_z_project",
    "register_stack",
    "report_path_from_output_path",
    "show_unmixed_channels_in_napari",
    "unmix",
    "unmix_ch0_from_ch1",
]
