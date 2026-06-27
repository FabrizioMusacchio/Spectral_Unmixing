Full 3D+t example
=================

This tutorial documents the interactive script
``user_scripts/unmix_full_TZCYX_synthetic_example.py``. It extends the basic
two-channel workflow to a full canonical ``TZCYX`` stack with multiple time
points and multiple z-slices.


How to use this tutorial
------------------------

The script is meant to be run as an interactive Python script in an editor that
supports cell-based execution.

The recommended workflow is:

1. open ``user_scripts/unmix_full_TZCYX_synthetic_example.py``,
2. run the cells from top to bottom,
3. adapt the parameters that matter for your own data.

The subsections below follow the same order as the script.


What this tutorial covers
-------------------------

Compared with :doc:`usage_unmix_example`, this script adds three important
practical ingredients:

- a full ``T > 1`` and ``Z > 1`` stack,
- a synthetic dataset with known bleed-through,
- and a per-time-point alpha-estimation example.

That makes it a good bridge between the simplest public example and a real
time-lapse z-stack workflow.


Imports
-------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # import the required helper functions:
   :end-before: # %% INPUT AND OUTPUT PATHS

The imported helpers are the same as in the simpler two-channel tutorial:

- ``unmix(...)`` for the actual spectral unmixing,
- ``report_path_from_output_path(...)`` for inspecting the JSON sidecar,
- ``show_unmixed_channels_in_napari(...)`` for visual inspection.


Define input and output paths
-----------------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # define the input path to the example dataset:
   :end-before: # %% FIXED ALPHA EXAMPLE

What you would usually change here:

- ``INPUT_PATH`` to your own `OMIO <https://omio.readthedocs.io/en/latest/>`_-readable stack,
- the output directory or filename pattern.


Spectral unmixing with a manually set alpha
-------------------------------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # define the output path for the fixed-alpha unmixing result:
   :end-before: # %% REFERENCE-TIME-POINT ALPHA EXAMPLE (mean-ratio)

This is the simplest baseline on the full stack. The exact same fixed alpha is
applied across all time points and all z-slices.

The most relevant parameters are:

- ``method="manual"``
- ``alpha``
- optional ``source_channel`` and ``target_channel``
- optional ``clip_negative``

This is the best starting point when you already have a coefficient from a
proper control measurement or an empirical estimate from a similar dataset. 


``mean_ratio`` on a full ``TZCYX`` stack
----------------------------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # define the output path for the reference-time-point unmixing result:
   :end-before: # %% REFERENCE-TIME-POINT LINEAR-FIT EXAMPLE

This is the simplest automatic estimator on the synthetic full stack. One alpha
is estimated from the selected reference time point, but that estimate uses all
z-slices belonging to that time point.

The parameters that matter most are:

- ``alpha_reference_t``: defines the reference time point for alpha estimation
- ``signal_percentile``: defines the  ...
- ``background_percentile``: ...
- ``target_low_percentile``: ...
- ``preprocess_alpha_inputs``: ...

This example is especially useful because it shows how those parameters behave
when the input really is a full time-lapse z-stack rather than just a 2D image.


``linear_fit`` on a full ``TZCYX`` stack
----------------------------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # define the output path for the reference-time-point linear-fit unmixing result:
   :end-before: # %% REFERENCE-TIME-POINT CORR-MIN EXAMPLE

This variant replaces the ratio-of-means estimator by a masked least-squares
fit without intercept.

Use this when you want a fit-based coefficient but still want the same general
reference-time-point workflow. In practice, the most important settings remain
the mask and preprocessing parameters:

- ``signal_percentile``
- ``background_percentile``
- optional ``target_low_percentile``
- optional ``preprocess_alpha_inputs``


``corr_min`` on a full ``TZCYX`` stack
--------------------------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # define the output path for the reference-time-point corr-min unmixing result:
   :end-before: # %% REFERENCE-TIME-POINT MI-MIN EXAMPLE

This method chooses alpha so that residual correlation between source and
corrected target becomes minimal.

The most relevant tuning parameters are:

- ``alpha_max``
- ``signal_percentile``
- ``background_percentile``
- optional ``target_low_percentile``
- optional ``max_alpha_voxels`` and ``random_state``

Because the optimization is performed on the full-stack reference volume, this
example is useful for seeing how aggressive correlation-minimization can become
on richer data.


``mi_min`` on a full ``TZCYX`` stack
------------------------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # define the output path for the reference-time-point mi-min unmixing result:
   :end-before: # %% PER-TIME-POINT ALPHA EXAMPLE

This method uses the two-channel PICASSO-like mutual-information criterion.

The most important settings are:

- ``mi_bins``
- ``alpha_max``
- ``signal_percentile``
- ``background_percentile``
- optional ``target_low_percentile``
- optional ``max_alpha_voxels`` and ``random_state``

Compared with the simpler estimators, this is often the slowest option, but it
can be useful when residual dependence remains visible after simpler
corrections.


Per-time-point alpha estimation
-------------------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # define the output path for the per-time-point unmixing result:
   :end-before: # %% END

This is the main extra method-specific feature of this tutorial. Instead of
estimating one global alpha, the script derives one alpha per time point and
applies each value to the corresponding time slice.

The key settings are:

- ``alpha_mode="per_t"``
- ``method``
- ``signal_percentile``
- ``background_percentile``
- optional ``target_low_percentile``

Use this when one global coefficient seems too rigid because brightness or gain
changes over time. At the same time, be careful: if the underlying biology also
changes over time, the alpha estimate itself can become time-dependent in a way
that is not purely optical.
