Tutorial: ``unmix_example.py``
==============================

This tutorial documents the interactive script
``user_scripts/unmix_example.py``. It is the best starting point for learning
the core ``unmix(...)`` workflow on a simple two-channel example.


What this tutorial covers
-------------------------

The script demonstrates:

- fixed-alpha correction,
- reference-time-point alpha estimation,
- several automatic alpha-estimation methods,
- and direct napari inspection after each run.


Cell 1: Imports
---------------

This cell imports the package functions used throughout the tutorial.

The most important imported helpers are:

- ``unmix(...)`` as the core workflow,
- ``report_path_from_output_path(...)`` to inspect the JSON sidecar,
- ``show_unmixed_channels_in_napari(...)`` to visualize the result.

.. literalinclude:: ../../user_scripts/unmix_example.py
   :language: python
   :start-after: # %% IMPORTS
   :end-before: # %% INPUT AND OUTPUT PATHS


Cell 2: Input and output paths
------------------------------

This cell defines the input dataset and the output filenames used by the later
examples.

What a user typically changes here:

- ``INPUT_PATH`` to point to a new OMIO-readable microscopy stack,
- the output directory or naming scheme,
- and, if desired, the source/target channel assumptions in later cells.

.. literalinclude:: ../../user_scripts/unmix_example.py
   :language: python
   :start-after: # %% INPUT AND OUTPUT PATHS
   :end-before: # %% FIXED ALPHA EXAMPLE


Cell 3: Fixed alpha
-------------------

This is the simplest and often most scientifically interpretable mode when a
coefficient has been measured from a suitable control dataset.

The important knobs are:

- ``alpha``:
  manually supplied bleed-through coefficient
- ``alpha_mode="fixed"``
- optional visualization colormaps for napari

.. literalinclude:: ../../user_scripts/unmix_example.py
   :language: python
   :start-after: # %% FIXED ALPHA EXAMPLE
   :end-before: # %% REFERENCE-TIME-POINT ALPHA EXAMPLE (mean-ratio)


Cell 4: Reference-time-point alpha with ``mean_ratio``
------------------------------------------------------

This cell estimates one alpha from a selected reference time point and applies
it to the full stack.

Important parameters:

- ``alpha_reference_t``:
  which time point is used for estimation
- ``signal_percentile``:
  how strict the bright-source mask is
- ``background_percentile``:
  low-percentile background subtraction for estimation

Increasing ``signal_percentile`` makes the mask more selective. Lowering it
uses more voxels but may mix in less source-specific regions.

.. literalinclude:: ../../user_scripts/unmix_example.py
   :language: python
   :start-after: # %% REFERENCE-TIME-POINT ALPHA EXAMPLE (mean-ratio)
   :end-before: # %% REFERENCE-TIME-POINT LINEAR-FIT EXAMPLE


Cell 5: Reference-time-point ``linear_fit``
-------------------------------------------

This variant estimates alpha via masked least squares without an intercept.

Use this when you want a fit-based coefficient rather than the simpler ratio of
means. It relies on the same reference-time-point and mask concept as the
previous cell.

.. literalinclude:: ../../user_scripts/unmix_example.py
   :language: python
   :start-after: # %% REFERENCE-TIME-POINT LINEAR-FIT EXAMPLE
   :end-before: # %% REFERENCE-TIME-POINT CORR-MIN EXAMPLE


Cell 6: Reference-time-point ``corr_min``
-----------------------------------------

This method chooses alpha such that correlation between the source channel and
the corrected target channel is minimized.

Useful parameters include:

- ``alpha_max``:
  upper bound of the optimization interval
- the same mask-related parameters as above

This method can be more aggressive than simple ratio-based estimation when
source and target channels are biologically correlated.

.. literalinclude:: ../../user_scripts/unmix_example.py
   :language: python
   :start-after: # %% REFERENCE-TIME-POINT CORR-MIN EXAMPLE
   :end-before: # %% REFERENCE-TIME-POINT MI-MIN EXAMPLE


Cell 7: Reference-time-point ``mi_min``
---------------------------------------

This method uses a PICASSO-like criterion in the two-channel case by choosing
alpha to minimize mutual information between the source channel and the
corrected target channel.

Important user-tunable parameters:

- ``mi_bins``:
  histogram resolution for mutual-information estimation
- ``alpha_max``:
  search interval bound
- mask and background settings

.. literalinclude:: ../../user_scripts/unmix_example.py
   :language: python
   :start-after: # %% REFERENCE-TIME-POINT MI-MIN EXAMPLE
   :end-before: # %% END


What to change for your own data
--------------------------------

For a new dataset, the main things to adapt are:

- ``INPUT_PATH``
- ``source_channel`` and ``target_channel``
- your chosen alpha strategy:
  fixed, reference-time-point, or another method
- mask strictness via ``signal_percentile`` and
  ``target_low_percentile``

If a reliable single-label control exists, start with the fixed-alpha mode
first. It usually gives the most interpretable baseline.
