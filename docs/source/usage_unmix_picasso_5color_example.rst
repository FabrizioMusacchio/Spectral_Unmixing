Tutorial: ``unmix_picasso_5color_example.py``
=============================================

This tutorial documents the interactive script
``user_scripts/unmix_picasso_5color_example.py``.

It is the main public example for genuine multi-channel blind unmixing beyond
the 3-channel case.


Cell 1: Imports, paths, and inspection
--------------------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_5color_example.py
   :language: python
   :start-after: # %% IMPORTS
   :end-before: # %% PICASSO MATLAB-N EXAMPLE

This block defines the measured 5-channel example and opens it in napari with
one layer per channel.

The current file is already stored in canonical ``TZCYX`` order with the
measured images on the channel axis, so no special preparation step is needed.


Cell 2: ``matlab_n`` on five channels
-------------------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_5color_example.py
   :language: python
   :start-after: # %% PICASSO MATLAB-N EXAMPLE
   :end-before: # %% PICASSO SOURCE-SINK-N EXAMPLE

This cell demonstrates the explicit N-channel generalization of the MATLAB
PICASSO workflow on five channels.

The most influential parameters are:

- ``max_iter``
- ``step_size``
- ``qN``
- ``pixel_bin_size``
- ``alpha_clip``

These control how aggressively and how stably the iterative blind-unmixing
updates proceed.


Cell 3: ``source_sink_n`` on five channels
------------------------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_5color_example.py
   :language: python
   :start-after: # %% PICASSO SOURCE-SINK-N EXAMPLE
   :end-before: # %% END

This cell shows the broadest source-sink formulation on the full five-channel
example.

The script starts from a permissive configuration:

- all selected channels are allowed as sinks,
- no channels are declared neutral.

For real datasets, this should often be narrowed down. The most important user
decisions are therefore:

- which channels actually need correction,
- which channels should be protected as neutral,
- and whether a more explicit source-sink matrix is preferable for the
  biological problem at hand.
