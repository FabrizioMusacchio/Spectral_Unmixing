"""
Interactive VS Code user script for bidirectional spectral bleed-through correction.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

from spectral_unmixing import (
    report_path_from_output_path,
    show_unmixed_channels_in_napari,
    unmix)
# %% INPUT AND OUTPUT PATHS
"""Define the example input stack and all output targets used below.

In fact, you just need to set ``INPUT_PATH`` to your own data and the rest will be 
automatically generated in a subfolder of the input file's parent directory.
"""
# define the input path to the example dataset:
INPUT_PATH = (PROJECT_ROOT / "example_data" / "PICASSO_examples" / "bidirectional_example.tif")
OUTPUT_DIR = INPUT_PATH.parent / "unmixed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
# %% FIXED BIDIRECTIONAL ALPHA EXAMPLE
"""Run bidirectional unmixing with manually chosen fixed coefficients.

Method summary:

- ``bidirectional=True`` activates a two-direction model in which channel 0 may
  bleed into channel 1 and channel 1 may also bleed back into channel 0.
- ``alpha_mode="fixed"`` means that neither direction is estimated from the
  data.
- ``method="manual"`` documents that both directions use user-provided
  coefficients rather than automatic estimation.
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
- ``clip_negative``:
  If left at its default ``True``, negative corrected values in both channels
  are clipped to zero after inversion of the two-channel mixing model.
- ``output_dtype``:
  Controls the saved data type of the corrected stack. The default
  ``"float32"`` is usually the safest choice for preserving the corrected
  values.

Why this is useful:

- Best choice when both directional coefficients were measured from proper
  single-label controls acquired with the same imaging settings.
"""
# define the output path for the fixed bidirectional unmixing result:
OUTPUT_FIXED = OUTPUT_DIR / "bidirectional_unmixed_fixed_alpha.tif"

fixed_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_FIXED,
    bidirectional=True,
    # source_channel=0,  # default: 0
    # target_channel=1,  # default: 1
    alpha=0.60,
    alpha_reverse=0.50,
    alpha_mode="fixed",
    method="manual",
    # clip_negative=True,  # default: True
    # output_dtype="float32",  # default: "float32"
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
# %% MEAN-RATIO EXAMPLE
"""Estimate one forward and one reverse alpha from the same reference time point.

Method summary:

- ``method="mean_ratio"`` is used for the forward direction.
- If ``method_reverse`` is left commented out, the same ``mean_ratio`` method
  is reused automatically for the reverse direction.
- If ``alpha_mode`` is left unset, as shown here, the pipeline defaults to
  ``reference_t`` with ``alpha_reference_t=0`` because no manual ``alpha`` is
  supplied.

What can be adjusted:

- ``source_channel`` and ``target_channel``:
  Define the forward direction. The reverse direction is inferred
  automatically from that pairing.
- ``alpha_mode="reference_t"`` can be set in case of, e.g., multi-time point
  stacks. This mode estimates one coefficient per direction from a chosen
  reference time point and then applies both coefficients to the whole stack.
  You can also switch to ``alpha_mode="per_t"`` when a separate coefficient
  should be estimated for each time point.
- ``alpha_reference_t``:
  Reference time point used for both directions. Only relevant for
  ``alpha_mode="reference_t"``; default: ``0``.
- ``signal_percentile``:
  Controls how selective the forward bright-source mask is. Higher values keep
  only brighter source voxels; lower values include more voxels.
- ``background_percentile``:
  Controls the percentile-based background estimate subtracted before forward
  alpha estimation.
- ``target_low_percentile``:
  Optionally restricts forward alpha estimation to voxels with comparatively
  low target intensity, which can help suppress genuine reverse-direction
  structure in the mask.
- ``preprocess_alpha_inputs``:
  Enables or disables percentile-based background subtraction and clipping for
  alpha estimation.
- ``max_alpha_voxels`` and ``random_state``:
  Control optional subsampling of masked voxels before estimating the
  coefficient. This can matter for very large volumes.
- ``alpha_reverse``:
  If set, replaces automatic reverse-direction estimation by a fixed manual
  reverse coefficient while the forward direction still uses ``mean_ratio``.
- ``method_reverse``:
  Allows the reverse direction to use another estimation method instead of
  inheriting ``mean_ratio``.
- ``signal_percentile_reverse``, ``background_percentile_reverse``, and
  ``target_low_percentile_reverse``:
  Reverse-direction counterparts of the forward mask and preprocessing
  settings. They are only needed if the reverse direction should use different
  selection rules.

Effect of these settings:

- This is the easiest automatic bidirectional mode to start with, because all
  reverse-direction settings inherit the forward values unless overridden.
"""
# define the output path for the bidirectional reference-time-point mean-ratio result:
OUTPUT_REFERENCE = OUTPUT_DIR / "bidirectional_unmixed_reference_t0_mean_ratio.tif"

reference_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE,
    bidirectional=True,
    # source_channel=0,  # default: 0
    # target_channel=1,  # default: 1
    # alpha_mode="reference_t",
    # alpha_reference_t=0,
    method="mean_ratio",
    signal_percentile=50.0,
    background_percentile=1.0,
    # target_low_percentile=95.0,
    # preprocess_alpha_inputs=True,  # default: True
    # max_alpha_voxels=500_000,  # default
    # random_state=0,  # default
    # alpha_reverse=None,
    # method_reverse="mean_ratio",
    # signal_percentile_reverse=99.0,
    # background_percentile_reverse=0.5,
    # target_low_percentile_reverse=95.0,
)
print(reference_output)
print(report_path_from_output_path(reference_output).read_text(encoding="utf-8"))
show_unmixed_channels_in_napari(
    reference_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Bidirectional reference mean_ratio")
# %% LINEAR-FIT EXAMPLE
"""Estimate forward and reverse coefficients with masked least-squares fitting.

Method summary:

- ``method="linear_fit"`` is applied to the forward direction.
- The reverse direction inherits the same method unless ``method_reverse`` is
  set explicitly.
- Both coefficients are estimated from the same reference time point but may
  use different masking settings through the commented reverse parameters.
- As in the previous example, leaving ``alpha_mode`` unset makes the pipeline
  default to ``reference_t`` with ``alpha_reference_t=0``.

What can be adjusted:

- ``source_channel`` and ``target_channel``:
  Define the forward correction direction.
- ``alpha_mode`` and ``alpha_reference_t``:
  Control whether one fit is estimated from one reference time point or one fit
  per time point in time-lapse data.
- ``signal_percentile``:
  Forward bright-source mask threshold.
- ``target_low_percentile``:
  Optional additional forward mask constraint that keeps only target-dim voxels.
- ``background_percentile``:
  Forward background percentile used during alpha-estimation preprocessing.
- ``preprocess_alpha_inputs``:
  Switches the percentile-based preprocessing for alpha estimation on or off.
- ``max_alpha_voxels`` and ``random_state``:
  Control optional voxel subsampling before fitting.
- ``alpha_reverse``:
  Lets you replace reverse-direction estimation by a manually chosen reverse
  coefficient.
- ``method_reverse``:
  Lets the reverse direction use another method instead of inheriting
  ``linear_fit``.
- ``signal_percentile_reverse``:
  Optional reverse-direction bright-source mask threshold.
- ``background_percentile_reverse`` and ``target_low_percentile_reverse``:
  Optional reverse-direction counterparts of the forward preprocessing and mask
  settings.
"""
# define the output path for the bidirectional reference-time-point linear-fit result:
OUTPUT_REFERENCE_LINEAR_FIT = OUTPUT_DIR / "bidirectional_unmixed_reference_t0_linear_fit.tif"

reference_linear_fit_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE_LINEAR_FIT,
    bidirectional=True,
    # source_channel=0,  # default: 0
    # target_channel=1,  # default: 1
    method="linear_fit",
    #alpha_mode="reference_t",
    #alpha_reference_t=0,
    signal_percentile=50.0,
    background_percentile=1.0,
    # target_low_percentile=95.0,
    # preprocess_alpha_inputs=True,  # default: True
    # max_alpha_voxels=500_000,  # default
    # random_state=0,  # default
    # alpha_reverse=None,
    # method_reverse="linear_fit",
    # signal_percentile_reverse=98.0,
    # background_percentile_reverse=1.0,
    # target_low_percentile_reverse=None,
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
- With ``alpha_mode`` omitted, the pipeline again defaults to
  ``reference_t`` with ``alpha_reference_t=0`` unless you explicitly switch to
  ``per_t``.

What can be adjusted:

- ``source_channel`` and ``target_channel``:
  Define the forward correction direction.
- ``alpha_mode`` and ``alpha_reference_t``:
  Choose between one reference-time-point estimate and one estimate per time
  point.
- ``signal_percentile``, ``background_percentile``, and
  ``target_low_percentile``:
  Define the forward estimation mask and preprocessing.
- ``preprocess_alpha_inputs``:
  Enables or disables percentile-based preprocessing before the correlation
  minimization.
- ``alpha_max``:
  Forward optimization bound.
- ``max_alpha_voxels`` and ``random_state``:
  Control optional voxel subsampling before optimization.
- ``alpha_reverse``:
  Lets you replace reverse-direction estimation by a fixed reverse
  coefficient.
- ``alpha_max_reverse``:
  Optional reverse optimization bound if the reverse contamination is expected
  to be stronger or weaker.
- ``method_reverse``:
  Can be changed independently if the reverse direction should use another
  estimator.
- ``signal_percentile_reverse``, ``background_percentile_reverse``, and
  ``target_low_percentile_reverse``:
  Reverse-direction counterparts of the forward mask and preprocessing
  settings.
"""


# define the output path for the bidirectional reference-time-point corr-min result:
OUTPUT_REFERENCE_CORR_MIN = OUTPUT_DIR / "bidirectional_unmixed_reference_t0_corr_min.tif"

reference_corr_min_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE_CORR_MIN,
    bidirectional=True,
    # source_channel=0,  # default: 0
    # target_channel=1,  # default: 1
    method="corr_min",
    # alpha_mode="reference_t",
    # alpha_reference_t=0,
    signal_percentile=50.0,
    background_percentile=1.0,
    alpha_max=1.0,
    # target_low_percentile=95.0,
    # preprocess_alpha_inputs=True,  # default: True
    # max_alpha_voxels=500_000,  # default
    # random_state=0,  # default
    # alpha_reverse=None,
    # method_reverse="corr_min",
    # signal_percentile_reverse=99.0,
    # background_percentile_reverse=1.0,
    # target_low_percentile_reverse=None,
    # alpha_max_reverse=0.5,
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
- With ``alpha_mode`` omitted, the pipeline defaults to ``reference_t`` with
  ``alpha_reference_t=0`` unless ``per_t`` is explicitly requested.

What can be adjusted:

- ``source_channel`` and ``target_channel``:
  Define the forward correction direction.
- ``alpha_mode`` and ``alpha_reference_t``:
  Choose between one mutual-information-based estimate from a reference time
  point and one estimate per time point.
- ``signal_percentile``, ``background_percentile``, and
  ``target_low_percentile``:
  Define the forward estimation mask and preprocessing.
- ``preprocess_alpha_inputs``:
  Enables or disables percentile-based preprocessing before MI minimization.
- ``mi_bins`` and ``mi_bins_reverse``:
  Histogram resolution used by the mutual-information estimate.
- ``alpha_max`` and ``alpha_max_reverse``:
  Search bounds for the forward and reverse optimizations.
- ``max_alpha_voxels`` and ``random_state``:
  Control optional voxel subsampling before optimization.
- ``alpha_reverse``:
  Lets you replace reverse-direction estimation by a fixed reverse
  coefficient.
- ``method_reverse``:
  Can be changed independently if the reverse direction should use another
  estimator.
- ``signal_percentile_reverse``, ``background_percentile_reverse``, and
  ``target_low_percentile_reverse``:
  Reverse-direction counterparts of the forward mask and preprocessing
  settings.
"""
# define the output path for the bidirectional reference-time-point mi-min result:
OUTPUT_REFERENCE_MI_MIN = OUTPUT_DIR / "bidirectional_unmixed_reference_t0_mi_min.tif"

reference_mi_min_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE_MI_MIN,
    bidirectional=True,
    # source_channel=0,  # default: 0
    # target_channel=1,  # default: 1
    method="mi_min",
    # alpha_mode="reference_t",
    # alpha_reference_t=0,
    signal_percentile=50.0,
    background_percentile=1.0,
    alpha_max=1.0,
    mi_bins=64,
    # target_low_percentile=95.0,
    # preprocess_alpha_inputs=True,  # default: True
    # max_alpha_voxels=500_000,  # default
    # random_state=0,  # default
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
