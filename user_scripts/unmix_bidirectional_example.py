"""
Interactive VS Code user script for bidirectional spectral bleed-through correction.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

"""Import the helpers used throughout the bidirectional unmixing tutorial.

This cell adds the repository root to ``sys.path`` so the local
``spectral_unmixing`` package can be imported without installation.
The imported functions cover:

- ``unmix(...)`` for the actual bidirectional spectral bleed-through correction.
- ``report_path_from_output_path(...)`` for loading the JSON sidecar report.
- ``show_unmixed_channels_in_napari(...)`` for reusing one napari viewer and
  updating its layers after each run.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from spectral_unmixing import (
    report_path_from_output_path,
    show_unmixed_channels_in_napari,
    unmix,
)

# %% INPUT AND OUTPUT PATHS
"""Define the example input stack and all output paths used below.

What can be adjusted here:

- ``INPUT_PATH``:
  Path to the raw two-channel TIFF stack to be unmixed.
- ``OUTPUT_DIR``:
  Subfolder that will receive all bidirectional unmixing outputs and JSON
  reports.
- ``OUTPUT_*``:
  One output file per alpha-estimation strategy.

Effect of changing these settings:

- A different ``INPUT_PATH`` switches to another dataset.
- A different ``OUTPUT_DIR`` changes where the results are stored.
- Renaming an ``OUTPUT_*`` path changes only the saved filename, not the
  processing logic.
"""

INPUT_PATH = (
    PROJECT_ROOT
    / "example_data"
    / "PICASSO_examples"
    / "Quantitative analysis of unmixing Before unmixing.tif"
)
OUTPUT_DIR = INPUT_PATH.parent / "unmixed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FIXED = OUTPUT_DIR / "bidirectional_unmixed_fixed_alpha.tif"
OUTPUT_REFERENCE = OUTPUT_DIR / "bidirectional_unmixed_reference_t0_mean_ratio.tif"
OUTPUT_REFERENCE_LINEAR_FIT = OUTPUT_DIR / "bidirectional_unmixed_reference_t0_linear_fit.tif"
OUTPUT_REFERENCE_CORR_MIN = OUTPUT_DIR / "bidirectional_unmixed_reference_t0_corr_min.tif"
OUTPUT_REFERENCE_MI_MIN = OUTPUT_DIR / "bidirectional_unmixed_reference_t0_mi_min.tif"
OUTPUT_PER_T = OUTPUT_DIR / "bidirectional_unmixed_per_t_mean_ratio.tif"

# %% FIXED BIDIRECTIONAL ALPHA EXAMPLE
"""Run bidirectional unmixing with manually chosen fixed coefficients.

Method summary:

- ``bidirectional=True`` activates a two-direction model in which channel 0 may
  bleed into channel 1 and channel 1 may also bleed back into channel 0.
- ``alpha_mode="fixed"`` means that neither direction is estimated from the
  data.
- The forward coefficient is supplied via ``alpha``.
- The reverse coefficient can optionally be supplied via ``alpha_reverse``.
  If ``alpha_reverse`` remains ``None``, the pipeline reuses ``alpha`` for the
  reverse direction.

What can be adjusted:

- ``alpha``:
  Forward bleed-through from ``source_channel`` into ``target_channel``.
- ``alpha_reverse``:
  Reverse bleed-through from ``target_channel`` back into ``source_channel``.
- ``source_channel`` and ``target_channel``:
  Define the forward direction. The reverse direction is inferred
  automatically.

Why this is useful:

- Best choice when both directional coefficients were measured from proper
  single-label controls acquired with the same imaging settings.
"""

fixed_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_FIXED,
    bidirectional=True,
    alpha=0.60,
    alpha_reverse=0.50,
    alpha_mode="fixed",
    method="manual",
    # alpha_reverse=0.08,
    # method_reverse="manual",
    # signal_percentile_reverse=99.0,
    # background_percentile_reverse=1.0,
    # target_low_percentile_reverse=None,
    # alpha_max_reverse=1.0,
    # mi_bins_reverse=64,
)
print(fixed_output)
print(report_path_from_output_path(fixed_output).read_text(encoding="utf-8"))
show_unmixed_channels_in_napari(
    fixed_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Bidirectional fixed",
    source_colormap="cyan",
    target_colormap="yellow")
# %% REFERENCE-TIME-POINT MEAN-RATIO EXAMPLE
"""Estimate one forward and one reverse alpha from the same reference time point.

Method summary:

- ``alpha_mode="reference_t"`` estimates one coefficient per direction from a
  chosen reference time point and then applies both coefficients to the whole
  stack.
- ``method="mean_ratio"`` is used for the forward direction.
- If ``method_reverse`` is left commented out, the same ``mean_ratio`` method
  is reused automatically for the reverse direction.

What can be adjusted:

- ``alpha_reference_t``:
  Reference time point used for both directions.
- ``signal_percentile`` and ``target_low_percentile``:
  Control the forward-direction estimation mask.
- The commented ``*_reverse`` parameters allow a different reverse-direction
  mask if needed.

Effect of these settings:

- This is the easiest automatic bidirectional mode to start with, because all
  reverse-direction settings inherit the forward values unless overridden.
"""

reference_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE,
    bidirectional=True,
    alpha_mode="reference_t",
    method="mean_ratio",
    alpha_reference_t=0,
    signal_percentile=50.0,
    background_percentile=1.0,
    # alpha_reverse=None,
    # method_reverse="mean_ratio",
    #signal_percentile_reverse=99.0,
    # background_percentile_reverse=0.5,
    # target_low_percentile_reverse=95.0,
    # alpha_max_reverse=1.0,
    # mi_bins_reverse=64,
)
print(reference_output)
print(report_path_from_output_path(reference_output).read_text(encoding="utf-8"))
show_unmixed_channels_in_napari(
    reference_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Bidirectional reference mean_ratio")
# %% REFERENCE-TIME-POINT LINEAR-FIT EXAMPLE
"""Estimate forward and reverse coefficients with masked least-squares fitting.

Method summary:

- ``method="linear_fit"`` is applied to the forward direction.
- The reverse direction inherits the same method unless ``method_reverse`` is
  set explicitly.
- Both coefficients are estimated from the same reference time point but may
  use different masking settings through the commented reverse parameters.

What can be adjusted:

- ``signal_percentile``:
  Forward bright-source mask threshold.
- ``signal_percentile_reverse``:
  Optional reverse-direction bright-source mask threshold.
- ``background_percentile`` and ``background_percentile_reverse``:
  Optional direction-specific preprocessing percentiles.
"""

reference_linear_fit_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE_LINEAR_FIT,
    bidirectional=True,
    alpha_mode="reference_t",
    method="linear_fit",
    alpha_reference_t=0,
    signal_percentile=50.0,
    background_percentile=1.0,
    # alpha_reverse=None,
    # method_reverse="linear_fit",
    # signal_percentile_reverse=98.0,
    # background_percentile_reverse=1.0,
    # target_low_percentile_reverse=None,
    # alpha_max_reverse=1.0,
    # mi_bins_reverse=64,
)
print(reference_linear_fit_output)
print(report_path_from_output_path(reference_linear_fit_output).read_text(encoding="utf-8"))
show_unmixed_channels_in_napari(
    reference_linear_fit_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Bidirectional reference linear_fit")
# %% REFERENCE-TIME-POINT CORR-MIN EXAMPLE
"""Estimate forward and reverse coefficients by minimizing residual correlation.

Method summary:

- ``method="corr_min"`` searches for a coefficient that minimizes the
  correlation between the source channel and the corrected target channel.
- In bidirectional mode this optimization is done once for the forward
  direction and once for the reverse direction.

What can be adjusted:

- ``alpha_max``:
  Forward optimization bound.
- ``alpha_max_reverse``:
  Optional reverse optimization bound if the reverse contamination is expected
  to be stronger or weaker.
- ``method_reverse``:
  Can be changed independently if one direction should use another estimator.
"""

reference_corr_min_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE_CORR_MIN,
    bidirectional=True,
    alpha_mode="reference_t",
    method="corr_min",
    alpha_reference_t=0,
    signal_percentile=50.0,
    background_percentile=1.0,
    alpha_max=1.0,
    # alpha_reverse=None,
    # method_reverse="corr_min",
    # signal_percentile_reverse=99.0,
    # background_percentile_reverse=1.0,
    # target_low_percentile_reverse=None,
    # alpha_max_reverse=0.5,
    # mi_bins_reverse=64,
)
print(reference_corr_min_output)
print(report_path_from_output_path(reference_corr_min_output).read_text(encoding="utf-8"))
show_unmixed_channels_in_napari(
    reference_corr_min_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Bidirectional reference corr_min",)
# %% REFERENCE-TIME-POINT MI-MIN EXAMPLE
"""Estimate forward and reverse coefficients by minimizing mutual information.

Method summary:

- ``method="mi_min"`` uses a PICASSO-like two-channel criterion for each
  direction independently.
- ``mi_bins`` controls the histogram-based mutual-information estimate for the
  forward direction.
- The reverse direction can inherit that setup or use its own ``mi_bins`` and
  mask settings via the commented reverse arguments.

What can be adjusted:

- ``mi_bins`` and ``mi_bins_reverse``:
  Histogram resolution used by the mutual-information estimate.
- ``alpha_max`` and ``alpha_max_reverse``:
  Search bounds for the forward and reverse optimizations.
"""

reference_mi_min_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE_MI_MIN,
    bidirectional=True,
    alpha_mode="reference_t",
    method="mi_min",
    alpha_reference_t=0,
    signal_percentile=50.0,
    background_percentile=1.0,
    alpha_max=1.0,
    mi_bins=64,
    # alpha_reverse=None,
    # method_reverse="mi_min",
    # signal_percentile_reverse=99.0,
    # background_percentile_reverse=1.0,
    # target_low_percentile_reverse=None,
    # alpha_max_reverse=1.0,
    # mi_bins_reverse=32,
)
print(reference_mi_min_output)
print(report_path_from_output_path(reference_mi_min_output).read_text(encoding="utf-8"))
show_unmixed_channels_in_napari(
    reference_mi_min_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Bidirectional reference mi_min")
# %% END
