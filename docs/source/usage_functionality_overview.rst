Functionality overview
======================

This page summarizes the currently available unmixing and filtering functions
at a conceptual level before the tutorial pages show them in script form.


Spectral unmixing workflows
---------------------------

Two main unmixing entry points are available:

- ``spectral_unmixing.unmix(...)`` for directed linear spectral unmixing
- ``spectral_unmixing.unmix_picasso(...)`` for PICASSO-family blind unmixing


``unmix(...)``: directed linear unmixing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the main workflow for the common case in which one channel bleeds into
another channel and a direct correction model is desired.

Forward model
^^^^^^^^^^^^^

.. math::

   I_{\mathrm{source}} \approx S

.. math::

   I_{\mathrm{target}} \approx T + \alpha S

The corrected target channel is computed by:

.. math::

   I_{\mathrm{target, corrected}}^{\ast}
   =
   I_{\mathrm{target}} - \alpha I_{\mathrm{source}}

Optionally:

.. math::

   I_{\mathrm{target, corrected}}
   =
   \max \left(
   I_{\mathrm{target}} - \alpha I_{\mathrm{source}},
   0
   \right)

Main configuration dimensions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``alpha_mode`` controls **where** alpha is obtained from in multi-time-point stacks:
  ``fixed``, ``reference_t``, or ``per_t``.

  - ``fixed`` uses a user-provided scalar ``alpha`` for the full stack. Only relevant for single-time-point stacks and for ``method="manual"``.
  - ``reference_t`` estimates one scalar ``alpha`` from a chosen reference time point, using all z-slices at that time point. Relevant for multi-time-point stacks and for any estimation method.
  - ``per_t`` estimates one ``alpha`` value per time point, again using all z-slices for each time point. Relevant for multi-time-point stacks and for any estimation method.

- ``method`` controls **how** alpha is estimated:
  ``manual``, ``mean_ratio``, ``linear_fit``, ``corr_min``, or ``mi_min``.
- ``source_channel`` and ``target_channel`` define the directed correction.
- ``bidirectional=True`` enables simultaneous correction of both channels by
  solving a :math:`2 \times 2` linear mixing model.

Alpha-estimation methods
^^^^^^^^^^^^^^^^^^^^^^^^

The package currently supports:

- ``manual``:
  user-supplied fixed alpha
- ``mean_ratio``:
  ratio of mean intensities inside a bright-source mask
- ``linear_fit``:
  masked least-squares fit without intercept
- ``corr_min``:
  coefficient chosen to minimize residual source-target correlation after
  correction
- ``mi_min``:
  coefficient chosen to minimize mutual information after correction

Masking and preprocessing
^^^^^^^^^^^^^^^^^^^^^^^^^

Alpha estimation can use:

- low-percentile background subtraction,
- clipping of negative background-corrected values,
- a bright-source mask controlled by ``signal_percentile``,
- and an optional low-target constraint controlled by
  ``target_low_percentile``.


``unmix_picasso(...)``: PICASSO-family blind unmixing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This workflow is intended for blind unmixing scenarios in which mixing
relations are inferred from the data rather than imposed as a single manual
coefficient.

The implemented PICASSO-family workflows are motivated by the original
PICASSO paper:

   Seo, J., Sim, Y., Kim, J. et al. *PICASSO allows ultra-multiplexed
   fluorescence imaging of spatially overlapping proteins without reference
   spectra measurements*. Nature Communications 13, 2475 (2022).
   https://doi.org/10.1038/s41467-022-30168-z

Currently available implementations are:

- ``matlab_3c``:
  close Python port of the original MATLAB 3-channel PICASSO workflow
- ``matlab_n``:
  explicit N-channel generalization of the MATLAB-style iterative update scheme
- ``source_sink_n``:
  explicit source-sink multi-channel correction with either a full
  ``source_sink_matrix`` or the higher-level ``sink_channels`` /
  ``neutral_channels`` interface

PICASSO-family workflows still rely on linear channel-mixing assumptions. The
blind component lies in estimating those relations from the data.


Filtering and projection helpers
--------------------------------

The package also includes post-unmixing helpers for stack processing.

``apply_filters(...)``
~~~~~~~~~~~~~~~~~~~~~~

Supported filter types currently include:

- ``median``
- ``gaussian``
- ordered filter chains such as ``["median", "gaussian"]``

The filters are applied on canonical ``TZCYX`` stacks and work for:

- full ``T>1, Z>1`` stacks,
- ``T=1`` or ``Z=1`` edge cases,
- and per-time-point parameter lists for adaptive filtering across time.

At a conceptual level:

- the median filter replaces each pixel or voxel by the median within a local
  neighborhood
- the Gaussian filter performs weighted local smoothing with a Gaussian kernel

The implementation supports:

- 2D slice-wise operation
- optional 3D filtering via ``apply_3D=True``
- optional alternate filter settings for channel 2 through dedicated keyword
  arguments


``match_histograms_across_time(...)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This helper matches the intensity distribution of each time point to the
reference stack at ``t=0``. It is useful when brightness drifts over time and a
consistent visual or analytical scale is needed before later processing.


``max_z_project(...)``
~~~~~~~~~~~~~~~~~~~~~~

This helper computes a maximum-intensity projection over the z-axis while
preserving ``T`` and ``C``.

If :math:`V(t, z, c, y, x)` is the input stack, the projected output is:

.. math::

   P(t, c, y, x) = \max_{z \in Z_{\mathrm{sel}}} V(t, z, c, y, x)

where :math:`Z_{\mathrm{sel}}` is either the full z-range or a user-specified
``zrange``.


Registration helpers
--------------------

The package currently provides:

- ``register_stack(...)`` for time registration,
- ``correct_intra_stack_z_drift(...)`` for within-stack z-alignment problems.

``register_stack(...)`` supports:

- ``pystackreg``-based rigid registration,
- ``phase_cross_correlation`` from scikit-image,
- channel-specific registration references,
- optional z-projection selection,
- and optional median pre-/post-filtering for registration reference images.


Visualization helpers
---------------------

For napari-based inspection, the package provides:

- ``show_unmixed_channels_in_napari(...)`` for source/target two-channel views
- ``show_all_channels_in_napari(...)`` for multi-channel inspection

Both helpers reuse an open napari viewer when possible and reopen one when the
previous viewer has been closed.
