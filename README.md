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

`alpha_mode` answers the question:

> From which part of the dataset should the coefficient be obtained?

It does **not** determine how the coefficient is computed numerically.

### Alpha estimation methods
`method` controls how `alpha` is estimated once the relevant source and target
volumes have been chosen by `alpha_mode`.

Available methods are:

- `manual`
  Use a user-provided `alpha`. This is only meaningful for
  `alpha_mode="fixed"`.
- `mean_ratio`
  Estimate `alpha` as the ratio of mean target and source intensities inside a
  bright-source mask.
- `linear_fit`
  Estimate `alpha` by masked least-squares fitting without intercept.
- `corr_min`
  Estimate `alpha` by minimizing the correlation between the source channel and
  the corrected target channel.
- `mi_min`
  Estimate `alpha` by minimizing the mutual information between the source
  channel and the corrected target channel.

So the logic is:

- `alpha_mode` decides **where** alpha is estimated from
- `method` decides **how** alpha is estimated

The default estimation method is `mean_ratio`:

```python
method="mean_ratio"
```

#### Shared alpha-estimation preprocessing
The helper `spectral_unmixing.prepare_source_target_for_alpha(...)` implements a
shared optional preprocessing step for alpha estimation.

If `preprocess_alpha_inputs=True`, it:

- converts inputs to `float32`
- subtracts a low-percentile background estimate from both channels
- clips negative values to zero

Mathematically, if the raw source and target volumes are denoted by
\(X_{\mathrm{raw}}\) and \(Y_{\mathrm{raw}}\), the preprocessing step computes
background estimates

$$
b_X = \operatorname{percentile}(X_{\mathrm{raw}}, p_{\mathrm{bg}})
$$

and

$$
b_Y = \operatorname{percentile}(Y_{\mathrm{raw}}, p_{\mathrm{bg}})
$$

and then forms

$$
X = \max(X_{\mathrm{raw}} - b_X, 0)
$$

$$
Y = \max(Y_{\mathrm{raw}} - b_Y, 0)
$$

where \(p_{\mathrm{bg}}\) is the chosen `background_percentile`.

This preprocessing is used only for **estimating** `alpha`. The final image
correction is still applied to the measured working array inside the unmixing
pipeline.

#### Shared alpha mask
The helper `spectral_unmixing.make_alpha_mask(...)` creates the voxel mask used
for alpha estimation.

By default, the mask is defined from bright source voxels:

$$
\mathcal{M}
=
\left\{
i \;\middle|\; X_i \ge
\operatorname{percentile}(X, p_{\mathrm{sig}})
\right\}
$$

where \(p_{\mathrm{sig}}\) is `signal_percentile`.

Optionally, the mask can be restricted further to voxels with comparatively low
target intensity:

$$
\mathcal{M}
=
\mathcal{M}
\cap
\left\{
i \;\middle|\; Y_i \le
\operatorname{percentile}(Y, p_{\mathrm{target,low}})
\right\}
$$

where \(p_{\mathrm{target,low}}\) is `target_low_percentile`.

This can be useful when one wants to estimate bleed-through primarily from
voxels with strong source signal but as little genuine target signal as
possible.

If this stricter mask becomes too small, the implementation falls back to the
source-only mask when possible and records that behavior in the JSON report.

#### Method: `manual`
For `method="manual"`, no coefficient is estimated from the data.

The user supplies

$$
\alpha \ge 0
$$

directly, and the correction uses that fixed value.

This is the scientifically preferred mode when `alpha` has been determined from
a suitable single-label control measurement acquired with the same imaging
settings.

#### Method: `mean_ratio`
This is the original default behavior of the package.

After preprocessing and masking, let

$$
x_i = X_i \quad \text{for } i \in \mathcal{M}
$$

and

$$
y_i = Y_i \quad \text{for } i \in \mathcal{M}.
$$

Then the estimate is

$$
\hat{\alpha}_{\mathrm{mean\_ratio}}
=
\frac{\frac{1}{|\mathcal{M}|}\sum_{i \in \mathcal{M}} y_i}
{\frac{1}{|\mathcal{M}|}\sum_{i \in \mathcal{M}} x_i}.
$$

This estimator is simple and often stable, but it is not identical to a
least-squares fit.

#### Method: `linear_fit`
This method performs masked least-squares fitting **without intercept**:

$$
y_i \approx \alpha x_i
\qquad \text{for } i \in \mathcal{M}.
$$

The resulting estimator is

$$
\hat{\alpha}_{\mathrm{linear\_fit}}
=
\frac{\sum_{i \in \mathcal{M}} x_i y_i}
{\sum_{i \in \mathcal{M}} x_i^2}.
$$

No intercept is fitted, because the optional background subtraction is already
handled during preprocessing and because an intercept would blur the
interpretation of \(\alpha\) as a bleed-through coefficient.

#### Method: `corr_min`
This method chooses \(\alpha\) so that the corrected target channel becomes as
uncorrelated as possible with the source channel:

$$
Y^{(\alpha)} = Y - \alpha X.
$$

The estimate is obtained by solving

$$
\hat{\alpha}_{\mathrm{corr\_min}}
=
\arg\min_{0 \le \alpha \le \alpha_{\max}}
\operatorname{corr}(X, Y^{(\alpha)})^2.
$$

In practice, the implementation uses Pearson correlation and bounded scalar
optimization on the interval \([0, \alpha_{\max}]\).

This can be more aggressive than `mean_ratio` or `linear_fit`, especially when
the true biology in source and target channels is itself correlated.

#### Method: `mi_min`
This method follows the two-channel version of the PICASSO idea: choose
\(\alpha\) such that the statistical dependence between source and corrected
target becomes minimal.

Again define

$$
Y^{(\alpha)} = Y - \alpha X.
$$

Then the estimate is obtained from

$$
\hat{\alpha}_{\mathrm{mi\_min}}
=
\arg\min_{0 \le \alpha \le \alpha_{\max}}
\operatorname{MI}(X, Y^{(\alpha)}),
$$

where \(\operatorname{MI}\) denotes mutual information. The current
implementation uses a histogram-based mutual-information estimate with a user
controlled number of bins `mi_bins`.

The two-channel `mi_min` method is inspired by the PICASSO criterion, but it is
**not** the full multi-channel PICASSO algorithm.

### Multi-channel PICASSO-like blind unmixing
The package additionally provides a separate function:

```python
from spectral_unmixing import unmix_picasso
```

This implements an iterative multi-channel blind-unmixing workflow under the
assumption that:

- the number of measured channels equals the number of fluorophores
- the mixture is approximately linear
- one wants to reduce pairwise statistical dependence between reconstructed
  channels

The underlying linear model is

$$
I = M F,
$$

where

- \(I\) is the vector of measured channels
- \(F\) is the vector of latent fluorophore signals
- \(M\) is an unknown mixing matrix

The iterative PICASSO-like implementation starts from the measured channels and
repeatedly updates channel pairs using

$$
F_j \leftarrow F_j - a_{ij} F_i
$$

with

$$
a_{ij}
=
\arg\min_{0 \le a \le \alpha_{\max}}
\operatorname{MI}(F_i, F_j - a F_i).
$$

These pairwise updates induce an estimated unmixing matrix \(U\), which is then
applied to the selected channels of the full stack.

This is a PICASSO-like iterative blind-unmixing criterion. It is not a
deep-learning method and it is conceptually separate from the simpler
two-channel `unmix(...)` workflow.

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
