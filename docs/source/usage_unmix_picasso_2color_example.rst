PICASSO 2-color example
=======================

This tutorial documents the interactive script
``user_scripts/unmix_picasso_2color_example.py``.

Even for two-channel data, ``unmix_picasso(...)`` can be useful when you want
to compare a blind-unmixing workflow against the standard directed
``unmix(...)`` approach.

The PICASSO-family blind-unmixing logic used here is motivated by the original
PICASSO paper:

   Seo, J., Sim, Y., Kim, J. et al. *PICASSO allows ultra-multiplexed
   fluorescence imaging of spatially overlapping proteins without reference
   spectra measurements*. Nature Communications 13, 2475 (2022).
   https://doi.org/10.1038/s41467-022-30168-z


How to use this tutorial
------------------------

The script is intended for interactive execution in an editor with cell support.

The recommended workflow is:

1. open ``user_scripts/unmix_picasso_2color_example.py``,
2. run the sections from top to bottom,
3. compare the blind-unmixing output against the simpler directed workflow if
   you have already worked through :doc:`usage_unmix_example`.


What this tutorial covers
-------------------------

This is the smallest PICASSO-family example in the repository. It shows two
conceptually different blind-unmixing strategies:

- ``implementation="matlab_n"``, an N-channel generalization of the MATLAB
  PICASSO workflow,
- ``implementation="source_sink_n"``, a more explicit source-sink formulation.


Imports
-------

.. literalinclude:: ../../user_scripts/unmix_picasso_2color_example.py
   :language: python
   :start-after: # %% IMPORTS
   :end-before: # %% INPUT AND OUTPUT PATHS


Define input and output paths
-----------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_2color_example.py
   :language: python
   :start-after: # define the input path to the example dataset:
   :end-before: # %% INSPECT PREPARED STACKS IN NAPARI

The main thing to adapt here is ``INPUT_PATH``. The script automatically
creates output paths in an ``unmixed`` subfolder next to the input file.


Inspect the measured channels
-----------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_2color_example.py
   :language: python
   :start-after: # inspect the stack in Napari:
   :end-before: # %% PICASSO MATLAB-N EXAMPLE

This inspection step is worth keeping in your own scripts. It lets you verify
that the channels were loaded as expected before any blind-unmixing is applied.


``matlab_n`` blind unmixing
---------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_2color_example.py
   :language: python
   :start-after: # define the output path for the PICASSO MATLAB-N unmixing result:
   :end-before: # %% PICASSO SOURCE-SINK-N EXAMPLE

This method applies the explicit N-channel generalization of the MATLAB-style
PICASSO iteration.

The most influential settings are:

- ``channels``
- ``implementation="matlab_n"``
- ``max_iter``
- ``step_size``
- ``qN``
- ``pixel_bin_size``
- ``alpha_clip``
- optional ``negativity_threshold`` and ``clip_every_n_iterations``

Even though some shared API options such as ``mi_bins`` and ``alpha_max`` are
still present in the call, the MATLAB-like implementation is primarily driven
by the MATLAB-style iteration parameters listed above.


``source_sink_n`` blind unmixing
--------------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_2color_example.py
   :language: python
   :start-after: # define the output path for the PICASSO source-sink-N unmixing result:
   :end-before: # %% END

This variant uses a more explicit source-sink description of the expected
cross-talk graph.

The most relevant settings are:

- ``sink_channels``
- ``neutral_channels``
- optional ``source_sink_matrix``
- ``alpha_max``
- ``mi_bins``
- ``max_alpha_voxels``

For two-channel data this is often the easier PICASSO-family mode to reason
about, because the user can state very directly which channel should be treated
as the sink and which one is allowed to act as the source.
