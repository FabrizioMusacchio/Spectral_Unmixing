Usage
=====

This section documents the practical, script-based usage model of
``spectral-unmixing``.

The package is designed to be driven from project-specific Python scripts that
can be executed cell by cell in the VS Code interactive window or a similar
notebook-like environment. The reusable package functions live in
``spectral_unmixing/``, while the repository's ``user_scripts/`` folder shows
how these functions are combined for real workflows.

.. toctree::
   :maxdepth: 2

   usage_example_datasets
   usage_functionality_overview
   usage_unmix_example
   usage_unmix_full_tzcyx_synthetic_example
   usage_unmix_bidirectional_example
   usage_unmix_picasso_2color_example
   usage_unmix_picasso_3color_example
   usage_unmix_picasso_5color_example
   usage_filter_and_register_stack


How to read the tutorials
-------------------------

Each tutorial page is based on one interactive script from the
repository's ``user_scripts/`` folder. The documentation follows 
the script cell by cell and explains:

- what each cell does,
- which parameters matter most,
- what a user can modify for a new dataset,
- and how those changes affect the unmixing result.


Recommended starting point
--------------------------

If you are new to the package, start in this order:

1. :doc:`usage_unmix_example`
2. :doc:`usage_unmix_full_tzcyx_synthetic_example`
3. :doc:`usage_filter_and_register_stack`
4. :doc:`usage_unmix_bidirectional_example`
5. the PICASSO-family blind-unmixing tutorials

The standard two-channel tutorial introduces the core package logic with the
least complexity. The  ``TZCYX`` tutorial then shows the same logic on
full time-lapse z-stacks. The helper tutorial then shows how filtering,
registration, histogram matching, and projection can be chained after
unmixing. The bidirectional tutorial adds a more advanced alpha-estimation
strategy for cases where the cross-talk is not strictly unidirectional. The
PICASSO tutorials build on the same ideas but add
blind coefficient estimation and multi-channel (more than two) configuration.
