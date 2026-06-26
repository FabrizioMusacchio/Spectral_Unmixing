"""
Interactive VS Code user script for spectral bleed-through correction.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

"""Import the helpers used throughout the interactive unmixing tutorial.

This cell wires the repository root into ``sys.path`` so the local
``spectral_unmixing`` package can be imported without installation.
The imported functions cover:

- ``unmix(...)`` for the actual spectral bleed-through correction.
- ``report_path_from_output_path(...)`` for loading the JSON sidecar report.
- ``show_unmixed_channels_in_napari(...)`` for reusing one napari viewer and
  updating its layers after each run.
"""

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
INPUT_PATH = PROJECT_ROOT / "example_data" / "synthetic_data" / "synthetic_bleedthrough_T9_Z20_C2.tif"
INPUT_NAME = INPUT_PATH.stem
OUTPUT_DIR = INPUT_PATH.parent / "unmixed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
# %% FIXED ALPHA EXAMPLE
"""Run unmixing with a manually chosen fixed bleed-through coefficient.

Method summary:

- ``alpha_mode="fixed"`` means no alpha is estimated from the data.
- ``method="manual"`` documents that the user provides ``alpha`` directly.
- The same alpha is applied to every time point and every z-slice.

What can be adjusted:

- ``alpha``:
  Strength of the source-channel subtraction from the target channel.
  Larger values remove more source contribution from the target channel.
- ``source_channel`` and ``target_channel``:
  Select which channel is treated as the bleeding source and which one is
  corrected.
- ``clip_negative`` inside ``unmix(...)``:
  If enabled, negative corrected values are clipped to zero.

When this is useful:

- Best choice when alpha was measured independently from a proper control.
- Most reproducible and scientifically preferred if acquisition settings are
- stable across experiments.
"""

OUTPUT_FIXED = OUTPUT_DIR / f"{INPUT_NAME}_unmixed_fixed_alpha.tif"

fixed_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_FIXED,
    #source_channel=0,  # default: 0
    #target_channel=1,  # default: 1
    alpha=0.62,
    alpha_mode="fixed",
    method="manual")

show_unmixed_channels_in_napari(
    fixed_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Fixed alpha",
    source_colormap="cyan",
    target_colormap="yellow")
# %% REFERENCE-TIME-POINT ALPHA EXAMPLE (mean-ratio)
"""Estimate one alpha from a reference time point using the ``mean_ratio`` rule.

Method summary:

- ``alpha_mode="reference_t"`` estimates one scalar alpha from a selected time
  point and reuses it for the full stack.
- ``method="mean_ratio"`` computes alpha as the mean target intensity divided
  by the mean source intensity within bright source voxels.
- All z-slices of the chosen reference time point contribute to that estimate.

What can be adjusted:

- ``alpha_reference_t``:
  Chooses which time point is used for alpha estimation.
- ``signal_percentile``:
  Defines how bright source voxels must be to enter the estimation mask.
  Higher values focus more strongly on the brightest source signal.
- ``background_percentile``:
  Controls the low-percentile background estimate subtracted before alpha
  estimation.
- ``target_low_percentile``:
  Optional extra restriction to prefer voxels with low target intensity, which
  can reduce bias from true biological target signal.

Effect of these settings:

- A higher ``signal_percentile`` usually makes alpha estimation more selective
  but also reduces the number of voxels used.
- A different ``alpha_reference_t`` matters when bleed-through or biology
  changes over time.
"""

OUTPUT_REFERENCE = OUTPUT_DIR / f"{INPUT_NAME}_unmixed_reference_t0_mean_ratio.tif"

reference_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE,
    #source_channel=0,  # default: 0
    #target_channel=1,  # default: 1
    alpha_mode="reference_t",
    method="mean_ratio",
    alpha_reference_t=0,
    signal_percentile=50.5,
    target_low_percentile=96.0,
    background_percentile=0.5,
    preprocess_alpha_inputs=False,
    clip_negative=True)
print(reference_output)
print(report_path_from_output_path(reference_output).read_text(encoding="utf-8"))
show_unmixed_channels_in_napari(
    reference_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Reference t0")
# %% REFERENCE-TIME-POINT LINEAR-FIT EXAMPLE
"""Estimate one alpha from a reference time point via masked least squares.

Method summary:

- ``method="linear_fit"`` fits the model ``target ≈ alpha * source`` inside
  the selected voxel mask.
- No intercept is fitted; background is handled by the percentile-based
  preprocessing used for alpha estimation.
- The resulting alpha is then applied to all time points and z-slices.

What can be adjusted:

- ``signal_percentile``:
  Controls which bright source voxels define the fitting mask.
- ``background_percentile``:
  Influences the background-subtracted data used before fitting.
- ``alpha_reference_t``:
  Selects the time point from which the fit is derived.

Effect of these settings:

- Compared with ``mean_ratio``, ``linear_fit`` is often a bit closer to a true
  least-squares estimate and can behave differently when masked intensities
  have broad dynamic ranges.
"""

OUTPUT_REFERENCE_LINEAR_FIT = OUTPUT_DIR / f"{INPUT_NAME}_unmixed_reference_t0_linear_fit.tif"

reference_linear_fit_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE_LINEAR_FIT,
    alpha_mode="reference_t",
    method="linear_fit",
    alpha_reference_t=0,
    signal_percentile=99.0,
    background_percentile=1.0)
print(reference_linear_fit_output)
print(report_path_from_output_path(reference_linear_fit_output).read_text(encoding="utf-8"))
show_unmixed_channels_in_napari(
    reference_linear_fit_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Reference linear_fit")
# %% REFERENCE-TIME-POINT CORR-MIN EXAMPLE
"""Estimate alpha by minimizing residual correlation after correction.

Method summary:

- ``method="corr_min"`` searches for the alpha that minimizes the Pearson
  correlation between the source channel and the corrected target channel.
- Intuition: after ideal bleed-through removal, source structure should be less
  visible inside the corrected target channel.

What can be adjusted:

- ``alpha_max``:
  Upper search bound for alpha. Increase it if stronger bleed-through is
  plausible; keep it conservative to avoid overly aggressive subtraction.
- ``signal_percentile`` and ``background_percentile``:
  Still control the source mask and preprocessing used for the optimization.

Effect of these settings:

- ``corr_min`` can be more aggressive than ``mean_ratio`` because it explicitly
  tries to remove statistical dependence.
- If source and target channels are biologically correlated, this method may
  subtract true target signal together with bleed-through.
"""

OUTPUT_REFERENCE_CORR_MIN = OUTPUT_DIR / f"{INPUT_NAME}_unmixed_reference_t0_corr_min.tif"

reference_corr_min_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE_CORR_MIN,
    #source_channel=0,  # default: 0
    #target_channel=1,  # default: 1
    alpha_mode="reference_t",
    method="corr_min",
    alpha_reference_t=0,
    signal_percentile=95.0,
    background_percentile=0.5,
    alpha_max=1.0,
    preprocess_alpha_inputs=True)
print(reference_corr_min_output)
print(report_path_from_output_path(reference_corr_min_output).read_text(encoding="utf-8"))
show_unmixed_channels_in_napari(
    reference_corr_min_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Reference corr_min")
# %% REFERENCE-TIME-POINT MI-MIN EXAMPLE
"""Estimate alpha by minimizing mutual information in a PICASSO-like way.

Method summary:

- ``method="mi_min"`` chooses alpha so that the mutual information between the
  source channel and the corrected target channel becomes as small as possible.
- This is a two-channel PICASSO-inspired criterion, not the full multi-channel
  blind-unmixing algorithm.

What can be adjusted:

- ``mi_bins``:
  Number of histogram bins used for the mutual-information estimate.
  More bins can capture finer structure but may become noisier.
- ``alpha_max``:
  Upper bound of the optimization range for alpha.
- ``max_alpha_voxels`` and ``random_state``:
  Control optional subsampling when very many voxels are available.
- ``signal_percentile`` and ``background_percentile``:
  Define the source mask and preprocessing used for the estimation.

Effect of these settings:

- ``mi_min`` can outperform simpler methods when residual nonlinear dependence
  is still visible after correction.
- It is also the slowest of the scalar alpha estimators and can be sensitive to
  histogram settings.
"""

OUTPUT_REFERENCE_MI_MIN = OUTPUT_DIR / f"{INPUT_NAME}_unmixed_reference_t0_mi_min.tif"

reference_mi_min_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_REFERENCE_MI_MIN,
    #source_channel=0,  # default: 0
    #target_channel=1,  # default: 1
    alpha_mode="reference_t",
    method="mi_min",
    alpha_reference_t=0,
    signal_percentile=50.0,
    background_percentile=1.0,
    preprocess_alpha_inputs=False,
    alpha_max=1.0,
    mi_bins=64)
print(reference_mi_min_output)
print(report_path_from_output_path(reference_mi_min_output).read_text(encoding="utf-8"))
show_unmixed_channels_in_napari(
    reference_mi_min_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Reference mi_min",)
# %% PER-TIME-POINT ALPHA EXAMPLE
"""Estimate one alpha per time point and correct each time point separately.

Method summary:

- ``alpha_mode="per_t"`` derives a separate alpha for each time point.
- Here the estimator is ``method="mean_ratio"``, but the same pattern can be
  combined with the other automatic alpha-estimation methods as well.
- Each estimated alpha uses all z-slices belonging to that time point.

What can be adjusted:

- ``method``:
  Swap in ``linear_fit``, ``corr_min`` or ``mi_min`` if needed.
- ``signal_percentile`` and ``background_percentile``:
  Influence the per-time-point alpha estimation exactly as in the
  reference-time-point examples.
- ``target_low_percentile``:
  Optional extra mask restriction when true target signal contaminates the
  estimation heavily.

Effect of these settings:

- Useful when illumination or gain changes over time and one global alpha would
  be too rigid.
- More flexible, but it can introduce time-dependent artifacts if biology
  changes in a way that biases the alpha estimate.
"""

OUTPUT_PER_T = OUTPUT_DIR / f"{INPUT_NAME}_unmixed_per_t_mean_ratio.tif"

per_t_output = unmix(
    input_path=INPUT_PATH,
    output_path=OUTPUT_PER_T,
    #source_channel=0,  # default: 0
    #target_channel=1,  # default: 1
    alpha_mode="per_t",
    method="mean_ratio",
    signal_percentile=99.0,
    background_percentile=1.0)
print(per_t_output)
print(report_path_from_output_path(per_t_output).read_text(encoding="utf-8"))
show_unmixed_channels_in_napari(
    per_t_output,
    source_channel=0,
    target_channel=1,
    layer_prefix="Per t")
# %% END
