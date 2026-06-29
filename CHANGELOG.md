## Spectral Unmixing changelog

See here for a detailed list of changes made in each release of *Spectral Unmixing*.
Please, also refer to the Repository [Releases page](https://github.com/FabrizioMusacchio/spectral-unmixing/releases).

Each release is also archived on Zenodo for long-term preservation and citation purposes: <https://doi.org/10.5281/zenodo.20933784>



---

## 🚀 spectral-unmixing v0.0.4

June 29, 2026

PICASSO source-sink refinement release for **spectral-unmixing**.

This release provides:

### PICASSO source-sink-N improvements

- Refined the `source_sink_n` blind-unmixing implementation to follow the PICASSO napari plugin's source-sink formulation more closely while remaining free of `torch` and neural-MINE dependencies.
- Changed the default source-sink optimization strategy from greedy pairwise fitting to joint optimization of all modeled sources contributing to a given sink.
- Added optional optimization of one small background-offset parameter per modeled source-sink relation.
- Added multi-start optimization for the source-sink workflow to improve robustness against poor local minima.
- Preserved the older greedy strategy as an explicit fallback via configuration.

### API and reproducibility

- Extended `unmix_picasso(...)` with new `source_sink_n` options:
  - `source_sink_optimize_background`
  - `source_sink_max_background`
  - `source_sink_n_restarts`
  - `source_sink_joint_optimization`
- Extended the PICASSO JSON sidecar reports to record the learned source-sink background parameters and the new optimization settings.

### Scripts and documentation

- Updated the public 2-color, 3-color, and 5-color PICASSO example scripts so their `source_sink_n` example cells expose the new optimization options directly.
- Updated the README, overview page, and PICASSO tutorials to describe the refined `source_sink_n` formulation mathematically and operationally.
- Clarified in the documentation that the current `source_sink_n` implementation is histogram-MI-based and napari-plugin-inspired, but not a neural or MINE-based reimplementation of the plugin.

---

## 🚀 spectral-unmixing v0.0.3

June 28, 2026

Documentation, dataset, and API refinement release for **spectral-unmixing**.

This release provides:

### Unmixing API refinement

- Changed the default `unmix(...)` behavior from an implicit `alpha_mode="fixed"` to `alpha_mode=None`.
- Added automatic alpha-mode resolution for `unmix(...)`:
  - if `alpha` is provided, the effective mode becomes `fixed`
  - if `alpha` is omitted and `method="manual"`, the pipeline now raises a clear error
  - if `alpha` is omitted and `method!="manual"`, the effective mode becomes `reference_t` with `alpha_reference_t=0`
- Preserved explicit user intent:
  - explicitly setting `alpha_mode="fixed"` still requires a user-provided `alpha`
  - explicitly setting `alpha_mode="reference_t"` or `alpha_mode="per_t"` still behaves exactly as requested
- Extended the JSON sidecar report with the requested and effective alpha-mode information to make this default resolution reproducible.

### Documentation

- Added a full Sphinx / Read the Docs documentation tree for the package.
- Expanded the overview and functionality pages with more complete mathematical and workflow descriptions.
- Added a broad tutorial collection covering:
  - standard two-channel unmixing
  - full `TZCYX` synthetic-stack unmixing
  - bidirectional unmixing
  - PICASSO-family blind unmixing for 2-channel, 3-channel, and 5-channel examples
  - filtering, registration, histogram matching, and projection helper workflows
- Improved tutorial style and consistency across the interactive scripts and the corresponding documentation pages.
- Added dataset provenance and license documentation for the public tutorial datasets.

### Public tutorial datasets

- Published a dedicated public example-dataset release on Zenodo for reproducing the tutorials:
  - Musacchio, F. (2026). *Example datasets for the spectral-unmixing pipeline* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.20984021
- Added transparent dataset-level README files for:
  - `Gockel_Nieves_Rivera_2026`
  - `PICASSO_examples`
  - `synthetic_data`

### Packaging and CI

- Added `docs/requirements.txt` so that Read the Docs builds can install the documentation dependencies correctly.
- Added and corrected a GitHub Actions test workflow for automated tests and Codecov integration.

### Testing and docs synchronization

- Updated unit tests to cover the new `alpha_mode=None` default behavior.
- Updated the documentation to describe the new default behavior consistently for both `T=1` and `T>1` stacks.

---

## 🚀 spectral-unmixing v0.0.2

June 26, 2026

PyPI and archival release of **spectral-unmixing**.

This release contains no code changes relative to `v0.0.1`.

It was created to:

- publish the package on PyPI via `pip install spectral-unmixing`
- create a follow-up GitHub release suitable for Zenodo archiving and citation
- align the public release metadata with the first official package publication

---

## 🚀 spectral-unmixing v0.0.1

June 26, 2026

First public main release of **spectral-unmixing**.

This initial release provides:

### Core unmixing

- A reusable `unmix(...)` pipeline for spectral bleed-through correction of microscopy TIFF stacks read with OMIO.
- Canonical `TZCYX` stack validation and processing throughout the main workflow.
- One-direction linear unmixing with correction of a user-selected target channel from a user-selected source channel.
- Optional bidirectional two-channel unmixing via inversion of the corresponding 2x2 linear mixing model.
- JSON sidecar reports written next to output TIFF files for reproducibility of parameters and estimated coefficients.
- Verbose terminal reporting during processing runs.

### Alpha estimation

- Fixed-alpha mode for user-supplied bleed-through coefficients.
- `reference_t` alpha mode to estimate one coefficient from a selected reference time point using all z-slices.
- `per_t` alpha mode to estimate one coefficient per time point.
- Multiple alpha-estimation methods:
  - `manual`
  - `mean_ratio`
  - `linear_fit`
  - `corr_min`
  - `mi_min`
- Shared optional preprocessing for alpha estimation, including percentile-based background subtraction and clipping.
- Optional reverse-direction parameter set for bidirectional unmixing, including inherited defaults when reverse values are omitted.

### PICASSO-family blind unmixing

- A separate `unmix_picasso(...)` workflow for multi-channel blind unmixing.
- `matlab_3c` implementation as a Python port of the original MATLAB 3-channel PICASSO workflow.
- `matlab_n` implementation as an explicit N-channel generalization of the MATLAB 3-channel workflow.
- `source_sink_n` implementation for direct source-sink multi-channel unmixing inspired by the napari PICASSO plugin.
- Support for explicit `source_sink_matrix` definitions.
- Support for higher-level `sink_channels` and `neutral_channels` configuration to auto-build source-sink relations without writing matrices manually.

### Filtering, projection, and registration add-ons

- `apply_filters(...)` with median and Gaussian filtering for `TZCYX` stacks.
- Support for applying multiple filters sequentially.
- Support for per-time-point filter strengths and optional alternate filtering settings for channel 2.
- `match_histograms_across_time(...)` for time-wise histogram matching against the `t=0` reference stack.
- `max_z_project(...)` with optional `zrange` selection and safe out-of-bounds clamping.
- `register_stack(...)` for time registration using either `pystackreg` or `phase_cross_correlation`.
- `correct_intra_stack_z_drift(...)` for optional within-stack z-drift correction.

### I/O and visualization

- OMIO-based TIFF reading and writing helpers.
- Output-path reconciliation for OMIO-written TIFF files, including `.ome.tif` edge cases.
- Shared napari viewer helpers for visualizing unmixed results without opening duplicate viewers on repeated calls.
- Configurable source and target colormaps in napari display helpers.

### Examples and tutorials

- Interactive example scripts for:
  - standard two-channel unmixing
  - bidirectional unmixing
  - PICASSO 2-channel, 3-channel, and 5-channel examples
  - filtering, projection, and registration workflows
- Synthetic example generation for full `TZCYX` test data with controlled bleed-through.

### Packaging, docs, and testing

- GPL-3.0 licensing.
- A package layout prepared for open-source distribution and future PyPI publication.
- Public README documentation focused on the core unmixing workflow, mathematical model, and available estimation methods.
- Publishable docstrings across the pipeline intended to integrate cleanly with future Sphinx / Read the Docs API documentation.
- Unit tests covering estimation, unmixing, filtering, registration, viewer reuse, and OMIO output handling.
