Basic usage example
==============================

This tutorial documents the interactive script
``user_scripts/unmix_example.py``. It is the best starting point for learning
the core ``unmix(...)`` workflow on a simple two-channel example.

How to use this tutorial
------------------------

The script is designed to be run as an interactive Python script, best
with an IDE that supports cell-based execution (e.g., Spyder, VSCode, or PyCharm). 

It is organized in cells, reflecting the structure of this tutorial.
Thus, the recommended way to follow this tutorial is:

1. download and open the ``user_scripts/unmix_example.py`` script,
2. run the cells from top to bottom,
3. adjust the configuration values that are relevant for your own data.

The subsections below follow the same order as the script cells.


What this tutorial covers
-------------------------

The script demonstrates:

- fixed-alpha correction,
- reference-time-point alpha estimation,
- several automatic alpha-estimation methods,
- and direct napari inspection after each run.


Imports
-------

The first cell imports the package functions used throughout the tutorial.

The most important imported helpers are:

- ``unmix(...)`` as the core workflow,
- ``report_path_from_output_path(...)`` to inspect the JSON sidecar,
- ``show_unmixed_channels_in_napari(...)`` to visualize the result.

.. literalinclude:: ../../user_scripts/unmix_example.py
   :language: python
   :start-after: # %% IMPORTS
   :end-before: # %% INPUT AND OUTPUT PATHS

.. note::

  ``PROJECT_ROOT = Path(__file__).resolve().parents[1]`` is used to locate the example dataset 
  and output directory relative to the project root. You can change this to a fixed path if you 
  want to run the script from a different working directory. You can completely remove this line 
  if you want to use absolute paths for input and output (see the next cell).


Define input and output paths
------------------------------

Next, we need to define the input dataset and the output filenames used by the 
examples:

.. literalinclude:: ../../user_scripts/unmix_example.py
   :language: python
   :start-after: # define the input path to the example dataset:
   :end-before: # %% FIXED ALPHA EXAMPLE

What you can change here:

- ``INPUT_PATH`` to point to a new `OMIO <https://omio.readthedocs.io/en/latest/>`_-readable microscopy stack, and
- the output directory or naming scheme.

Spectral unmixing with a manually set alpha
--------------------------------------------------

The first cell in the script demonstrates the usage of a user-provided, fixed alpha value for the full stack.
This is the simplest and often most scientifically interpretable mode when a
coefficient has been measured from a suitable control dataset or has been determined from prior knowledge.

The corresponding cell consists of three main steps:

1. set an ``OUTPUT_PATH`` (here called ``OUTPUT_FIXED``) for the corrected stack,
2. call ``unmix(...)`` with the fixed alpha,
3. visualize the result in napari using ``show_unmixed_channels_in_napari(...)``.

.. literalinclude:: ../../user_scripts/unmix_example.py
   :language: python
   :start-after: # define the output path for the fixed-alpha unmixing result:
   :end-before: # %% REFERENCE-TIME-POINT ALPHA EXAMPLE (mean-ratio)

The important knobs are:

- ``method="manual"``: this tells the workflow to use a user-provided alpha value,
- ``alpha``:
  manually supplied bleed-through coefficient
- ``alpha_mode="fixed"`` (relevant only for multi-time-point stacks)
- optional ``source_channel`` and ``target_channel`` to define the directed correction; default is 0 and 1, respectively, which means that channel 1 is corrected for bleed-through from channel 0.
- optional visualization colormaps for napari


``mean_ratio`` method
---------------------

By changing the method to ``mean_ratio``, the workflow estimates alpha from a selected reference 
time point (in case of a multi-time-point stack; otherwise single-time-point workflow is applied) 
and applies it to the full stack. Alpha is computed as the mean target intensity divided
by the mean source intensity within a mask that is defined by the ``signal_percentile`` and
 ``background_percentile`` parameters:


.. literalinclude:: ../../user_scripts/unmix_example.py
   :language: python
   :start-after: # define the output path for the reference-time-point unmixing result:
   :end-before: # %% REFERENCE-TIME-POINT LINEAR-FIT EXAMPLE

Important parameters:

- ``method="mean_ratio"``: this tells the workflow to estimate alpha from the ratio of means
- ``alpha_reference_t``:
  which time point is used for estimation in case of multi-time-point stacks; otherwise ignored and single-time-point workflow is applied
- ``signal_percentile``:
  how strict the bright-source mask is
- ``background_percentile``:
  low-percentile background subtraction for estimation

Increasing ``signal_percentile`` makes the mask more selective. Lowering it
uses more voxels but may mix in less source-specific regions.


``linear_fit`` method
---------------------

This variant estimates alpha via masked least squares without an intercept.

Use this when you want a fit-based coefficient rather than the simpler ratio of
means. It relies on the same reference-time-point (if relevant) and mask concept as the
previous cell:

.. literalinclude:: ../../user_scripts/unmix_example.py
   :language: python
   :start-after: # define the output path for the reference-time-point linear-fit unmixing result:
   :end-before: # %% REFERENCE-TIME-POINT CORR-MIN EXAMPLE


``corr_min`` method
-------------------

This method chooses alpha such that correlation between the source channel and
the corrected target channel is minimized.

Useful parameters include:

- ``alpha_max``:
  upper bound of the optimization interval
- the same mask-related parameters as above

This method can be more aggressive than simple ratio-based estimation when
source and target channels are biologically correlated:

.. literalinclude:: ../../user_scripts/unmix_example.py
   :language: python
   :start-after: # define the output path for the reference-time-point corr-min unmixing result:
   :end-before: # %% REFERENCE-TIME-POINT MI-MIN EXAMPLE


``mi_min`` method
-----------------

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
   :start-after: # define the output path for the reference-time-point mi-min unmixing result:
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
