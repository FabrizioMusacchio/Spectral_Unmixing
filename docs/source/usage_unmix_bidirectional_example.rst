Bidirectional unmixing example
==============================

This tutorial documents the interactive script
``user_scripts/unmix_bidirectional_example.py``.

It covers the case in which bleed-through can occur in both directions between
two channels, so a simple one-way subtraction is no longer sufficient.


How to use this tutorial
------------------------

The script is designed for cell-based execution in an interactive Python
environment.

The recommended workflow is:

1. open ``user_scripts/unmix_bidirectional_example.py``,
2. run the cells from top to bottom,
3. adapt the forward and reverse settings to your own dataset.

The subsections below follow the same order as the script.


Core idea
---------

In bidirectional mode, our package models:

.. math::

   I_0 = S_0 + \alpha_{10} S_1

.. math::

   I_1 = S_1 + \alpha_{01} S_0

and solves the corresponding :math:`2 \times 2` linear system jointly.

This is preferable to sequential subtraction, because sequential subtraction
would depend on the order in which the two corrections are applied.


Imports
-------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # %% IMPORTS
   :end-before: # %% INPUT AND OUTPUT PATHS


Define input and output paths
-----------------------------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # define the input path to the example dataset:
   :end-before: # %% FIXED BIDIRECTIONAL ALPHA EXAMPLE

As in the other examples, the main thing you would change here is
``INPUT_PATH``. The script automatically writes results into an ``unmixed``
subfolder next to the input file.


Bidirectional unmixing with fixed coefficients
----------------------------------------------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # define the output path for the fixed bidirectional unmixing result:
   :end-before: # %% MEAN-RATIO EXAMPLE

This is the bidirectional analogue of the standard fixed-alpha workflow.

The main parameters are:

- ``bidirectional=True``:
  activates the two-direction model instead of the simpler one-way
  source-to-target subtraction.
- ``alpha``:
  forward-direction bleed-through coefficient. Larger values remove more of the
  forward contamination.
- ``alpha_reverse``:
  reverse-direction coefficient. Larger values remove more contamination in the
  reverse direction.
- optional ``source_channel`` and ``target_channel``:
  define the forward direction; the reverse direction is inferred from that
  pairing.

This is the scientifically strongest option when both directional coefficients
come from proper single-label controls.


Bidirectional ``mean_ratio``
----------------------------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # define the output path for the bidirectional reference-time-point mean-ratio result:
   :end-before: # %% LINEAR-FIT EXAMPLE

The most important settings are:

- ``signal_percentile``:
  controls how selective the forward-direction source mask is. Higher values
  keep only brighter source voxels; lower values include more voxels.
- ``background_percentile``:
  controls the percentile-based background subtraction used before estimation.
  Higher values remove more baseline; lower values preserve more low signal.
- optional ``target_low_percentile``
- the commented ``*_reverse`` parameters if the reverse direction should use a
  different mask or method

This is usually the easiest automatic bidirectional mode to start with because
all reverse-direction settings inherit the forward values unless you override
them explicitly.

.. note::

  The bidirectional workflow generally accepts a different estimation method. 
  If you want to use a different method in the reverse direction, you can set
  ``method_reverse`` to one of the supported methods. If it is left at ``None`` 
  or omitted, the reverse direction will inherit the forward method. The same inheritance
  applies to other parameters as well.

The example script shown here uses a stack with only one time point, so
``alpha_mode`` and ``alpha_reference_t`` are left commented out in the code. For 
real multi-time-point stacks, you would usually
set ``alpha_mode="reference_t"`` explicitly when one shared coefficient per
direction should be estimated from a chosen reference time point, or
``alpha_mode="per_t"`` when one forward and one reverse coefficient should be
estimated separately for each time point. In these cases, additional relevant parameters 
are:

- ``alpha_mode``: ``reference_t`` or ``per_t`` for multi-time-point stacks.
  The former estimates one alpha from the reference time point; the latter
  estimates one alpha per time point. In the present ``T=1`` example, these
  arguments are commented out because the default ``reference_t`` behavior with
  ``t=0`` is already sufficient.
- ``alpha_reference_t``:
  defines the reference time point from which both directional coefficients are
  estimated. Only relevant when ``alpha_mode="reference_t"``.


Bidirectional ``linear_fit``
----------------------------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # define the output path for the bidirectional reference-time-point linear-fit result:
   :end-before: # %% REFERENCE-TIME-POINT CORR-MIN EXAMPLE

This variant estimates forward and reverse coefficients with masked least
squares. 

The most relevant settings are:

- ``signal_percentile``
- ``background_percentile``
- optional ``target_low_percentile``
- optional ``method_reverse``
- optional reverse-direction mask parameters such as
  ``signal_percentile_reverse``

This is a good choice when you want a fit-based estimate in both directions but
still want to retain the same reference-time-point logic.


Bidirectional ``corr_min``
--------------------------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # define the output path for the bidirectional reference-time-point corr-min result:
   :end-before: # %% REFERENCE-TIME-POINT MI-MIN EXAMPLE

Here the forward and reverse coefficients are chosen by minimizing residual
correlation after correction. 


The settings most worth tuning are:

- ``alpha_max``:
  upper bound for the forward optimization. Larger values allow stronger
  subtraction; smaller values keep the fit more conservative.
- optional ``alpha_max_reverse``
- ``signal_percentile``
- ``background_percentile``
- optional ``method_reverse``

This can be more aggressive than ratio- or fit-based estimation, especially
when the two channels are themselves biologically correlated.


Bidirectional ``mi_min``
------------------------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # define the output path for the bidirectional reference-time-point mi-min result:
   :end-before: # %% END

This method uses the two-channel PICASSO-like mutual-information criterion in
both directions. 

The key settings are:

- ``mi_bins``:
  histogram resolution for the mutual-information objective. Higher values can
  capture finer structure but can become less stable; lower values are coarser
  and often more robust.
- optional ``mi_bins_reverse``
- ``alpha_max``
- optional ``alpha_max_reverse``
- mask and preprocessing settings for each direction

This is the most flexible of the scalar bidirectional estimators, but also the
slowest and potentially the most aggressive if the channels are biologically
linked.


What to tune most often
-----------------------

For real bidirectional datasets, the most important additional controls beyond
the standard one-way workflow are:

- ``alpha_reverse``
- ``method_reverse``
- ``signal_percentile_reverse``
- ``background_percentile_reverse``
- ``target_low_percentile_reverse``
- ``alpha_max_reverse``
- ``mi_bins_reverse``

Any reverse parameter left at ``None`` falls back to the corresponding forward
setting. That inheritance model is convenient for a first pass, but for real
data it is often worth checking whether the reverse direction deserves its own
masking or optimization settings.
