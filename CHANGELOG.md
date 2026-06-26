## Spectral Unmixing changelog

See here for a detailed list of changes made in each release of *Spectral Unmixing*.
Please, also refer to the Repository [Releases page](https://github.com/FabrizioMusacchio/spectral-unmixing/releases).

Each release is also archived on Zenodo for long-term preservation and citation purposes: *Will be provided soon.*



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
