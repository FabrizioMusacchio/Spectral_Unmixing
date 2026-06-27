Tutorial: ``unmix_picasso_3color_example.py``
=============================================

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


Cell 1: Imports, paths, and first inspection
--------------------------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_3color_example.py
   :language: python
   :start-after: # %% IMPORTS
   :end-before: # %% PICASSO MATLAB-N EXAMPLE

This block imports the necessary helpers, defines the input/output paths, and
opens the measured 3-channel stack in napari.


Cell 2: ``matlab_n``
--------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_3color_example.py
   :language: python
   :start-after: # %% PICASSO MATLAB-N EXAMPLE
   :end-before: # %% PICASSO MATLAB-3C EXAMPLE

This cell runs the general N-channel version of the MATLAB-style workflow on
three channels.

Use this mode when you want the MATLAB-style iterative logic but prefer to keep
the same code path across both 3-channel and larger-N examples.


Cell 3: ``matlab_3c``
---------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_3color_example.py
   :language: python
   :start-after: # %% PICASSO MATLAB-3C EXAMPLE
   :end-before: # %% PICASSO SOURCE-SINK-N EXAMPLE

This cell runs the explicit 3-channel MATLAB port.

This is the implementation to use when you want the closest available Python
analogue of the original 3-channel MATLAB code.


Cell 4: ``source_sink_n``
-------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_3color_example.py
   :language: python
   :start-after: # %% PICASSO SOURCE-SINK-N EXAMPLE
   :end-before: # %% END

This cell demonstrates source-sink configuration on a realistic 3-channel case.

The important idea here is that the user can describe the expected bleed-through
graph in a readable way:

- which channels should be corrected as sinks,
- which channels should remain neutral,
- and which channels are therefore allowed to act as sources.

This is often easier to reason about than writing the full source-sink matrix
manually.
