Tutorial: ``unmix_picasso_2color_example.py``
=============================================

This tutorial documents the interactive script
``user_scripts/unmix_picasso_2color_example.py``.

Even for two-channel data, ``unmix_picasso(...)`` can be useful when you want
to compare a blind-unmixing workflow against the standard directed
``unmix(...)`` approach.


Cell 1: Imports and path setup
------------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_2color_example.py
   :language: python
   :start-after: # %% IMPORTS
   :end-before: # %% INSPECT PREPARED STACKS IN NAPARI

This block imports ``unmix_picasso(...)``, the JSON-report helper, and the
multi-channel napari viewer helper. It then defines the input and output paths.


Cell 2: Inspect the measured stack
----------------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_2color_example.py
   :language: python
   :start-after: # %% INSPECT PREPARED STACKS IN NAPARI
   :end-before: # %% PICASSO MATLAB-N EXAMPLE

This cell opens the measured stack with one napari layer per channel. It is a
sanity-check step and also shows the recommended inspection pattern for
multi-channel example data.


Cell 3: ``matlab_n`` blind unmixing
-----------------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_2color_example.py
   :language: python
   :start-after: # %% PICASSO MATLAB-N EXAMPLE
   :end-before: # %% PICASSO SOURCE-SINK-N EXAMPLE

This cell runs the explicit N-channel generalization of the MATLAB PICASSO
workflow on a two-channel example.

Key parameters to tune:

- ``implementation="matlab_n"``
- ``max_iter``
- ``step_size``
- ``qN``
- ``pixel_bin_size``
- ``alpha_clip``


Cell 4: ``source_sink_n`` blind unmixing
----------------------------------------

.. literalinclude:: ../../user_scripts/unmix_picasso_2color_example.py
   :language: python
   :start-after: # %% PICASSO SOURCE-SINK-N EXAMPLE
   :end-before: # %% END

This cell demonstrates the source-sink formulation.

For two-channel data, the most important conceptual decision is whether you want
to model:

- one channel as sink and the other as source, or
- both channels as non-neutral participants.

The script uses the higher-level ``sink_channels`` / ``neutral_channels``
interface rather than a manually written source-sink matrix.
