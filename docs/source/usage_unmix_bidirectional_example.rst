Tutorial: ``unmix_bidirectional_example.py``
============================================

This tutorial documents the interactive script
``user_scripts/unmix_bidirectional_example.py``.

It covers the case in which bleed-through can occur in both directions between
two channels, so a simple one-way subtraction is no longer sufficient.


Core idea
---------

In bidirectional mode, the package models:

.. math::

   I_0 = S_0 + \alpha_{10} S_1

.. math::

   I_1 = S_1 + \alpha_{01} S_0

and solves the corresponding 2x2 linear system jointly.


Cell 1: Imports
---------------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # %% IMPORTS
   :end-before: # %% INPUT AND OUTPUT PATHS


Cell 2: Input and output paths
------------------------------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # %% INPUT AND OUTPUT PATHS
   :end-before: # %% FIXED BIDIRECTIONAL ALPHA EXAMPLE

This cell defines the example dataset and the output paths for each
bidirectional variant.


Cell 3: Fixed bidirectional coefficients
----------------------------------------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # %% FIXED BIDIRECTIONAL ALPHA EXAMPLE
   :end-before: # %% REFERENCE-TIME-POINT MEAN-RATIO EXAMPLE

This cell is the most direct bidirectional analogue of the standard fixed-alpha
workflow.

Important parameters:

- ``bidirectional=True``
- ``alpha``:
  forward coefficient
- ``alpha_reverse``:
  reverse coefficient

If ``alpha_reverse`` is omitted, the package can reuse the forward alpha, but
that should only be done when you have a reason to believe the two directions
are comparable.


Cell 4: Reference-time-point ``mean_ratio``
-------------------------------------------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # %% REFERENCE-TIME-POINT MEAN-RATIO EXAMPLE
   :end-before: # %% REFERENCE-TIME-POINT LINEAR-FIT EXAMPLE

This cell estimates one forward and one reverse coefficient from the same
reference time point.


Cell 5: Reference-time-point ``linear_fit``
-------------------------------------------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # %% REFERENCE-TIME-POINT LINEAR-FIT EXAMPLE
   :end-before: # %% REFERENCE-TIME-POINT CORR-MIN EXAMPLE

This variant estimates forward and reverse coefficients with masked least
squares.


Cell 6: Reference-time-point ``corr_min``
-----------------------------------------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # %% REFERENCE-TIME-POINT CORR-MIN EXAMPLE
   :end-before: # %% REFERENCE-TIME-POINT MI-MIN EXAMPLE


Cell 7: Reference-time-point ``mi_min``
---------------------------------------

.. literalinclude:: ../../user_scripts/unmix_bidirectional_example.py
   :language: python
   :start-after: # %% REFERENCE-TIME-POINT MI-MIN EXAMPLE
   :end-before: # %% END

This cell uses the PICASSO-like two-channel mutual-information criterion in
both directions.


What to tune
------------

The most important bidirectional controls are:

- ``alpha_reverse``
- ``method_reverse``
- ``signal_percentile_reverse``
- ``background_percentile_reverse``
- ``target_low_percentile_reverse``
- ``alpha_max_reverse``
- ``mi_bins_reverse``

Any reverse parameter left at ``None`` falls back to the corresponding forward
setting.
