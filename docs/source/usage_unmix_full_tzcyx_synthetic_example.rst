Tutorial: ``unmix_full_TZCYX_synthetic_example.py``
===================================================

This tutorial documents the interactive script
``user_scripts/unmix_full_TZCYX_synthetic_example.py``. It extends the standard
two-channel unmixing tutorial to a full canonical ``TZCYX`` stack with multiple
time points and multiple z-planes.


Why this tutorial matters
-------------------------

This script is useful because it demonstrates the same unmixing logic on a
stack that is closer to a real structural-imaging dataset:

- ``T > 1``
- ``Z > 1``
- ``C = 2``

It also uses a synthetic dataset with known bleed-through, which makes it a
safe and reproducible starting point for testing.


Cell 1: Imports
---------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # %% IMPORTS
   :end-before: # %% INPUT AND OUTPUT PATHS

This cell imports the same main helpers as the simpler two-channel tutorial,
but on a synthetic ``TZCYX`` stack.


Cell 2: Input and output paths
------------------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # %% INPUT AND OUTPUT PATHS
   :end-before: # %% FIXED ALPHA EXAMPLE

This cell points to the synthetic full-stack example and defines output paths
for all following runs.


Cell 3: Fixed alpha
-------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # %% FIXED ALPHA EXAMPLE
   :end-before: # %% REFERENCE-TIME-POINT ALPHA EXAMPLE (mean-ratio)

This cell runs the simplest correction on the full ``TZCYX`` stack.

What it demonstrates:

- the unmixing logic applies over all time points and all z-planes,
- the source channel stays unchanged,
- and only the target channel is corrected.


Cell 4: Reference-time-point ``mean_ratio``
-------------------------------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # %% REFERENCE-TIME-POINT ALPHA EXAMPLE (mean-ratio)
   :end-before: # %% REFERENCE-TIME-POINT LINEAR-FIT EXAMPLE

This cell estimates one alpha from ``t=0`` using all z-planes at that time
point and then applies that same value to the entire stack.


Cell 5: Reference-time-point ``linear_fit``
-------------------------------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # %% REFERENCE-TIME-POINT LINEAR-FIT EXAMPLE
   :end-before: # %% REFERENCE-TIME-POINT CORR-MIN EXAMPLE

This cell shows the least-squares alternative on the same full stack.


Cell 6: Reference-time-point ``corr_min``
-----------------------------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # %% REFERENCE-TIME-POINT CORR-MIN EXAMPLE
   :end-before: # %% REFERENCE-TIME-POINT MI-MIN EXAMPLE

This cell demonstrates optimization-based alpha selection using correlation
minimization.


Cell 7: Reference-time-point ``mi_min``
---------------------------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # %% REFERENCE-TIME-POINT MI-MIN EXAMPLE
   :end-before: # %% PER-TIME-POINT ALPHA EXAMPLE

This cell demonstrates mutual-information minimization on a full time-lapse
z-stack.


Cell 8: Per-time-point alpha
----------------------------

.. literalinclude:: ../../user_scripts/unmix_full_TZCYX_synthetic_example.py
   :language: python
   :start-after: # %% PER-TIME-POINT ALPHA EXAMPLE
   :end-before: # %% END

This is the most important extra cell relative to the simpler tutorial. It
shows how to estimate a separate alpha for each time point while still using
all z-planes at each ``t``.

Use this when:

- intensity drifts over time,
- one global alpha appears too rigid,
- but you still want to preserve the full ``TZCYX`` structure.

Be cautious, however: if source and target biology both change over time, a
per-time-point estimate can also introduce time-dependent artifacts.
