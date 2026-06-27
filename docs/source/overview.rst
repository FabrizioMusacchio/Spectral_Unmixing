Overview
========

.. figure:: _static/blead_throug_1600px.png
   :alt: Schematic illustration of spectral bleed-through
   :align: center
   :figwidth: 90%

   **Figure 1.** Schematic illustration of spectral bleed-through in
   fluorescence imaging setups. Shown are four different fluorophores (blue,
   green, yellow, red) and their respective emission spectra as a function of
   wavelength. Bleed-through occurs when the emission of one fluorophore is
   also detected in another channel, for example when green emission is
   partially recorded in the yellow channel. This can happen if detection
   windows are not cleanly separated, if the selected detection range is too
   broad, or if the emission peaks of two fluorophores lie too close together.
   Source:
   `fabriziomusacchio.com <https://www.fabriziomusacchio.com/teaching/teaching_bioimage_analysis/09_napari_bleach_correction>`_
   (license: CC BY-NC-SA 4.0)


.. figure:: _static/spectra_unmixing_example.jpg
   :alt: Example of spectral unmixing in a two-channel fluorescence image
   :align: center
   :figwidth: 100%

   **Figure 2.** Example of spectral unmixing in a two-channel 2D fluorescence
   image. Panel **a)** shows the spectrally mixed input image, in which signal
   from the cyan channel bleeds into the magenta channel. Panel **b)** shows
   the result after correction with ``spectral-unmixing`` using a fixed alpha,
   which improves the separation of the two channels and reduces false-positive
   magenta signal originating from cyan bleed-through. Source image data:
   `figshare dataset <https://figshare.com/articles/figure/PICASSO_allows_ultra-multiplexed_fluorescence_imaging_of_spatially_overlapping_proteins_without_reference_spectra_measurements/19596682/1?file=34810114>`_
   (CC BY 4.0), from Seo, J., Sim, Y., Kim, J. et al. *PICASSO allows
   ultra-multiplexed fluorescence imaging of spatially overlapping proteins
   without reference spectra measurements*. Nature Communications 13, 2475
   (2022). https://doi.org/10.1038/s41467-022-30168-z.


What the package does
---------------------

``spectral-unmixing`` focuses on spectral bleed-through correction in
multi-dimensional microscopy stacks. It supports:

- linear two-channel unmixing via ``unmix(...)``
- optional bidirectional two-channel correction via inversion of a 2x2 linear
  mixing model
- multiple alpha-estimation strategies for linear unmixing
- PICASSO-family blind unmixing via ``unmix_picasso(...)``
- optional follow-up helpers for filtering, projection, registration, and
  napari visualization


Stack model and file formats
----------------------------

The package assumes that OMIO returns stacks in canonical ``TZCYX`` order:

- ``T``: time
- ``Z``: z-plane
- ``C``: channel
- ``Y``: image height
- ``X``: image width

This means the workflows support:

- full time-lapse z-stacks,
- single-time-point 3D stacks,
- 2D multi-channel images with ``Z=1``,
- and mixed cases such as ``T>1, Z=1``.

Input is intentionally format-agnostic on the package side. Any microscopy file
format currently supported by `OMIO <https://omio.readthedocs.io/en/latest/>`_
can be used as input. Output stacks are written back through OMIO, typically as
TIFF or OME-TIFF.


Linear spectral unmixing
------------------------

The core linear model implemented by ``unmix(...)`` is:

.. math::

   I_{\mathrm{source}} \approx S

.. math::

   I_{\mathrm{target}} \approx T + \alpha S

Here, :math:`S` is the true source signal, :math:`T` is the true target signal,
and :math:`\alpha` describes how strongly the source channel contaminates the
target channel.

The corrected target channel is obtained by:

.. math::

   I_{\mathrm{target, corrected}}^{\ast}
   =
   I_{\mathrm{target}} - \alpha I_{\mathrm{source}}

Optionally, negative values are clipped:

.. math::

   I_{\mathrm{target, corrected}}
   =
   \max \left(
   I_{\mathrm{target}} - \alpha I_{\mathrm{source}},
   0
   \right)

This is a linear unmixing workflow. The package does **not** implement a
nonlinear unmixing model.


Blind unmixing in the PICASSO family
------------------------------------

In addition to direct linear correction, the package provides
``unmix_picasso(...)`` for PICASSO-family blind unmixing. In this context,
"blind" means that mixing relations are estimated from the measured data rather
than supplied as fixed reference spectra.

Three implementation paths are available:

- ``matlab_3c``:
  a close Python port of the original MATLAB 3-channel PICASSO workflow
- ``matlab_n``:
  an explicit N-channel generalization of that 3-channel workflow
- ``source_sink_n``:
  a source-sink formulation in which selected channels are corrected as sinks
  from one or more modeled source channels

These workflows still use linear channel-mixing assumptions. What differs is
how the mixing coefficients are inferred.


Optional helper modules
-----------------------

The package also includes secondary helper functionality for post-unmixing
processing:

- ``apply_filters(...)``:
  median and Gaussian filtering for canonical ``TZCYX`` stacks
- ``match_histograms_across_time(...)``:
  time-wise histogram matching against the first time point
- ``max_z_project(...)``:
  z-projection while preserving ``T`` and ``C``
- ``register_stack(...)``:
  time registration using either ``pystackreg`` or
  ``phase_cross_correlation``
- ``correct_intra_stack_z_drift(...)``:
  optional within-stack z-drift correction
- ``show_unmixed_channels_in_napari(...)`` and
  ``show_all_channels_in_napari(...)``:
  reusable napari inspection helpers

The current documentation focuses primarily on the unmixing workflows. The
filtering and registration modules are already available and are summarized in
the usage documentation, but their dedicated tutorial pages can be expanded
later.


Reproducibility
---------------

Every unmixing run writes a JSON sidecar report next to the output stack. This
report stores the parameters used for the run, including the selected mode,
channel assignments, estimated coefficients, and implementation settings. This
makes the package suitable for transparent, script-based microscopy workflows.


Where to start
--------------

If you are new to the package, a good reading order is:

1. :doc:`installation`
2. :doc:`usage_functionality_overview`
3. :doc:`usage_unmix_example`
4. the more specialized tutorial pages in :doc:`usage`
