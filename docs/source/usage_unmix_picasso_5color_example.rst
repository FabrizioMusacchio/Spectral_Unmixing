PICASSO 5-color example
=======================

This tutorial documents the interactive script
``user_scripts/unmix_picasso_5color_example.py``.

It is the main public example for genuine multi-channel blind unmixing beyond
the 3-channel case.

The PICASSO-family blind-unmixing logic used here is motivated by the original
PICASSO paper:

   Seo, J., Sim, Y., Kim, J. et al. *PICASSO allows ultra-multiplexed
   fluorescence imaging of spatially overlapping proteins without reference
   spectra measurements*. Nature Communications 13, 2475 (2022).
   https://doi.org/10.1038/s41467-022-30168-z


How to use this tutorial
------------------------

The script is intended for interactive execution in a cell-based editor.

The recommended workflow is:

1. open ``user_scripts/unmix_picasso_5color_example.py``,
2. inspect the measured channels first,
3. compare the generalized MATLAB-style path with the source-sink path.


What this tutorial is good for
------------------------------

This is the best public example in the repository for understanding how the
package behaves when there are many channels and when the cross-talk graph is
not obvious in advance.


Imports
-------

.. literalinclude:: ../../user_scripts/unmix_picasso_5color_example.py
   :language: python
   :start-after: # %% IMPORTS
   :end-before: # %% INPUT AND OUTPUT PATHS


Define input and output paths
-----------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_5color_example.py
   :language: python
   :start-after: # define the input path to the example dataset:
   :end-before: # inspect the stack in Napari:


Inspect the measured channels
-----------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_5color_example.py
   :language: python
   :start-after: # inspect the stack in Napari:
   :end-before: # %% PICASSO MATLAB-N EXAMPLE

As in the smaller PICASSO tutorials, it is worth inspecting the raw measured
channels first before deciding which blind-unmixing strategy is the better fit
for the dataset.


``matlab_n`` blind unmixing on five channels
--------------------------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_5color_example.py
   :language: python
   :start-after: # define the output path for the PICASSO MATLAB-N unmixing result:
   :end-before: # %% PICASSO SOURCE-SINK-N EXAMPLE

This method applies the explicit N-channel generalization of the MATLAB-style
PICASSO iteration to all five selected channels.

The settings that matter most are:

- ``channels``:
  selects which measured channels are included in the blind-unmixing run.
- ``implementation="matlab_n"``:
  chooses the generalized MATLAB-style iteration for the selected channels.
- ``max_iter``:
  number of iteration steps. More iterations can unmix more strongly but may
  also increase instability.
- ``step_size``:
  strength of each update step. Larger values make the update more aggressive;
  smaller values make it gentler.
- ``qN``:
  quantization parameter for the MATLAB-style mutual-information estimate.
- ``pixel_bin_size``:
  spatial binning factor before the mutual-information calculation. Larger
  values smooth more strongly; smaller values preserve more detail.
- ``alpha_clip``:
  clipping bound for pairwise coefficients. Larger values allow stronger
  pairwise subtraction.
- optional ``negativity_threshold`` and ``clip_every_n_iterations``:
  control negativity handling during the iterative updates.

This is usually the best choice when you want a broad, symmetric blind-unmixing
strategy across many channels without encoding a very explicit source-sink
model yourself.


``source_sink_n`` blind unmixing on five channels
-------------------------------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_5color_example.py
   :language: python
   :start-after: # define the output path for the PICASSO source-sink-N unmixing result:
   :end-before: # %% END

This method uses the more explicit source-sink formulation.

The most important settings are:

- ``sink_channels``:
  defines which channels should actually be corrected as sinks.
- ``neutral_channels``:
  defines which channels should stay untouched and be excluded from active
  correction roles.
- optional ``source_sink_matrix``:
  gives explicit manual control over the allowed source-to-sink graph.
- ``alpha_max``:
  upper bound for source-to-sink coefficients. Larger values allow stronger
  subtraction; smaller values keep the estimate more conservative.
- ``mi_bins``:
  histogram resolution for the mutual-information objective.
- ``max_alpha_voxels``:
  optional cap on the voxel count used for coefficient estimation. Lower values
  speed up estimation; higher values use more of the image data.

For real five-channel data this mode is often attractive because it lets you
move from a broad first-pass model to a more selective graph once you have a
better idea which channels are plausible sinks and which should remain neutral.
