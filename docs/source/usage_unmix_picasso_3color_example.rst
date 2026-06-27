PICASSO 3-color example
=======================

This tutorial documents the interactive script
``user_scripts/unmix_picasso_3color_example.py``.

It is the most direct comparison point between:

- the original 3-channel MATLAB PICASSO workflow,
- the explicit N-channel generalization,
- and the source-sink formulation.

The PICASSO-family methods documented here are motivated by the original
PICASSO paper:

   Seo, J., Sim, Y., Kim, J. et al. *PICASSO allows ultra-multiplexed
   fluorescence imaging of spatially overlapping proteins without reference
   spectra measurements*. Nature Communications 13, 2475 (2022).
   https://doi.org/10.1038/s41467-022-30168-z


How to use this tutorial
------------------------

The script is designed for interactive execution.

The recommended workflow is:

1. open ``user_scripts/unmix_picasso_3color_example.py``,
2. inspect the measured channels first,
3. compare the three different blind-unmixing implementations on the same
   stack.


Imports
-------

.. literalinclude:: ../../user_scripts/unmix_picasso_3color_example.py
   :language: python
   :start-after: # %% IMPORTS
   :end-before: # %% INPUT AND OUTPUT PATHS


Define input and output paths
-----------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_3color_example.py
   :language: python
   :start-after: # define the input path to the example dataset:
   :end-before: # %% INSPECT PREPARED STACKS IN NAPARI


Inspect the measured channels
-----------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_3color_example.py
   :language: python
   :start-after: # inspect the stack in Napari:
   :end-before: # %% PICASSO MATLAB-N EXAMPLE

Before any unmixing is applied, the script opens the measured stack in napari.
This is useful for confirming channel order and for building intuition about
which channels appear to contaminate which others.


``matlab_n`` blind unmixing
---------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_3color_example.py
   :language: python
   :start-after: # define the output path for the PICASSO MATLAB-N unmixing result:
   :end-before: # %% PICASSO MATLAB-3C EXAMPLE

This method keeps the MATLAB-style iterative logic, but uses the generalized
``matlab_n`` implementation rather than the strict 3-channel port.

The settings that matter most are:

- ``channels``:
  selects which measured channels participate in the blind-unmixing run.
- ``implementation="matlab_n"``:
  chooses the generalized MATLAB-style path rather than the strict
  three-channel port.
- ``max_iter``:
  number of iterative updates. More iterations can unmix more strongly but may
  also become less stable.
- ``step_size``:
  strength of each update step. Larger values are more aggressive; smaller
  values are more conservative.
- ``qN``:
  quantization parameter for the MATLAB-style mutual-information estimate.
- ``pixel_bin_size``:
  spatial binning factor before the mutual-information calculation. Larger
  values smooth more strongly; smaller values preserve more detail.
- ``alpha_clip``:
  clipping bound for pairwise coefficients. Larger values allow stronger
  pairwise subtraction.
- optional ``negativity_threshold`` and ``clip_every_n_iterations``:
  control how intermediate negativity is monitored and how often positivity
  enforcement is applied.

This is the best choice when you want MATLAB-like behavior but also want the
same conceptual workflow to scale to larger channel counts.


``matlab_3c`` blind unmixing
----------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_3color_example.py
   :language: python
   :start-after: # define the output path for the PICASSO MATLAB-3C unmixing result:
   :end-before: # %% PICASSO SOURCE-SINK-N EXAMPLE

This is the closest available Python analogue of the original MATLAB
three-channel workflow.

Use this mode when exact three-channel comparison to the classic PICASSO logic
is more important than having a unified N-channel code path.

The main settings remain the MATLAB-style iteration parameters:

- ``max_iter``
- ``step_size``
- ``qN``
- ``pixel_bin_size``
- ``alpha_clip``


``source_sink_n`` blind unmixing
--------------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_3color_example.py
   :language: python
   :start-after: # define the output path for the PICASSO source-sink-N unmixing result:
   :end-before: # %% END

This variant uses an explicit source-sink description of the expected bleed-
through graph.

The most relevant settings are:

- ``sink_channels``:
  defines which channels should be corrected as sinks.
- ``neutral_channels``:
  defines which channels should remain untouched and not act as sinks.
- optional ``source_sink_matrix``:
  gives explicit manual control over the allowed bleed-through graph.
- ``alpha_max``:
  upper bound for source-to-sink coefficients. Larger values allow stronger
  subtraction; smaller values constrain the fit more strongly.
- ``mi_bins``:
  histogram resolution for the mutual-information objective.
- ``max_alpha_voxels``:
  optional cap on the number of voxels used for coefficient estimation. Lower
  values speed up the fit; higher values use more data.

This is often the easiest mode to reason about biologically, because the user
can describe which channels should actually be cleaned and which channels
should remain untouched or act only as sources.
