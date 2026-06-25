# Spectral unmixing for microscopy image stacks

`spectral-unmixing` is a small Python package focused on spectral bleed-through correction in microscopy TIFF stacks read with OMIO. The package assumes OMIO's canonical axis order `TZCYX`, so downstream code can safely address time, z, channel, and spatial axes in a predictable way.

The main goal of the project is reproducible spectral unmixing. Additional modules for filtering, registration, and projection are included as optional helpers for further image processing, but they are intentionally secondary to the unmixing workflow.

## Installation

```bash
pip install -e .
```

Core dependencies include:

- `numpy`
- `omio-microscopy`
- `scipy`
- `scikit-image`
- `pystackreg`

## Spectral unmixing model
The implemented correction assumes that one channel contributes linearly to another channel:

$$
I_{\text{source}} \approx S
$$

$$
I_{\text{target}} \approx T + \alpha \, S
$$

Here,

- $I_{\text{source}}$ is the measured intensity in the source channel
- $I_{\text{target}}$ is the measured intensity in the target channel
- $S$ is the true source-channel signal
- $T$ is the true target-channel signal
- $\alpha$ is the bleed-through coefficient from source into target

Under this model, the source channel contaminates the target channel linearly.
The actual unmixing step therefore subtracts the estimated source contribution
from the measured target signal:

$$
I_{\text{target, corrected}}^{\ast}
= I_{\text{target}} - \alpha \, I_{\text{source}}
$$

This is the core linear spectral unmixing equation.

In practice, the subtraction can produce negative values in pixels or voxels
where the estimated bleed-through contribution is slightly larger than the
measured target intensity. Since negative intensities are not physically
meaningful for the final corrected image, the pipeline can optionally apply a
second, purely post-processing step:

$$
I_{\text{target, corrected}}
= \max\!\left(
I_{\text{target}} - \alpha \, I_{\text{source}},
0
\right)
$$

So these are not two different models. They are two consecutive steps:

1. linear unmixing by subtraction
2. optional clipping of negative values to zero

Only the chosen target channel is corrected. The source channel is left
unchanged.

## Core unmixing function
The main entry point is `spectral_unmixing.unmix(...)`.

```python
from spectral_unmixing import unmix

output_path = unmix(
    input_path="input.tif",
    output_path="unmixed/input_unmixed_reference_t0.tif",
    alpha_mode="reference_t",
    alpha_reference_t=0,
    source_channel=0,
    target_channel=1,
    signal_percentile=99.0,
)
```

`unmix(...)` performs the following steps:

- reads the input stack with OMIO
- validates that the image is in canonical `TZCYX` order
- obtains `alpha` either as a fixed value or by estimation from the data
- applies the correction to the target channel over all `T` and `Z`
- clips negative corrected values to zero if requested
- writes the corrected TIFF stack
- writes a JSON sidecar report next to the output TIFF for reproducibility

The function returns the path to the written TIFF file.

### Alpha modes
Three `alpha_mode` values are available:

- `fixed`
  Use a user-provided scalar `alpha` for the full stack.
- `reference_t`
  Estimate one scalar `alpha` from a chosen reference time point, using all
  z-slices at that time point.
- `per_t`
  Estimate one `alpha` value per time point, again using all z-slices for each
  time point.

The helper `spectral_unmixing.estimate_alpha_from_volume(...)` implements the
actual estimation on matching source and target volumes. It:

- converts inputs to `float32`
- subtracts a low-percentile background estimate from both channels
- clips negative values to zero
- creates a source-signal mask from bright source voxels
- estimates `alpha` from the ratio of masked mean target and source intensities

### Output and reproducibility
Each unmixing run writes a JSON sidecar report next to the output TIFF, for
example:

```text
input_unmixed_reference_t0.tif.json
```

This report stores the main processing settings such as alpha mode, estimated
alpha values, source and target channels, axis order, and output dtype.

Terminal progress output is enabled by default and can be disabled with
`verbose=False`.

### Scientific Note
A fixed `alpha` measured from a proper single-label control recording is scientifically preferable.

Estimating `alpha` from the mixed experimental stack is available as a pragmatic first-pass workflow, but it can be biased when source and target biology overlap spatially.

`alpha_mode="reference_t"` assumes that the bleed-through factor is stable across time.

`alpha_mode="per_t"` can compensate for slow intensity changes, but may also introduce time-dependent artifacts when biology changes over time.

## Add-ons: Filtering, registration, and projection
In addition to spectral unmixing, the package also includes optional helper
functions for:

- filtering: `apply_filters(...)`
- time-wise histogram matching: `match_histograms_across_time(...)`
- max-z projection: `max_z_project(...)`
- intra-stack z-drift correction: `correct_intra_stack_z_drift(...)`
- time registration: `register_stack(...)`

These helper modules are meant to support follow-up image processing after unmixing. They are not the primary focus of the project, and their full documentation will be expanded later in Read the Docs.

For now, see the tutorial-style user scripts:

- [user_scripts/unmix_ch0_from_ch1_interactive.py](/Users/husker/Science/Python/Projekte/Spectral%20Unmixing/user_scripts/unmix_ch0_from_ch1_interactive.py)
- [user_scripts/filter_and_project_stack.py](/Users/husker/Science/Python/Projekte/Spectral%20Unmixing/user_scripts/filter_and_project_stack.py)
- [user_scripts/filter_and_register_stack.py](/Users/husker/Science/Python/Projekte/Spectral%20Unmixing/user_scripts/filter_and_register_stack.py)
- [user_scripts/fine_filter_and_register_stack.py](/Users/husker/Science/Python/Projekte/Spectral%20Unmixing/user_scripts/fine_filter_and_register_stack.py)
