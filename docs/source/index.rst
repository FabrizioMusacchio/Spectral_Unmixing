Spectral Unmixing Documentation
===============================

.. figure:: _static/spectra_unmixing_logo.png
   :alt: Spectral Unmixing logo
   :align: center
   :figwidth: 55%


.. image:: https://badgen.net/badge/icon/GitHub%20repository?icon=github&label
   :target: https://github.com/FabrizioMusacchio/spectral-unmixing
   :alt: GitHub Repository

.. image:: https://img.shields.io/pypi/v/spectral-unmixing.svg
   :target: https://pypi.org/project/spectral-unmixing/
   :alt: PyPI version

.. image:: https://img.shields.io/badge/License-GPL%20v3-green.svg
   :target: https://github.com/FabrizioMusacchio/spectral-unmixing/blob/main/LICENSE
   :alt: GPLv3 License

.. image:: https://img.shields.io/badge/Zenodo%20Archive-10.5281%2Fzenodo.20933784-blue
   :target: https://doi.org/10.5281/zenodo.20933784
   :alt: Zenodo Archive


`spectral-unmixing <https://github.com/FabrizioMusacchio/spectral-unmixing>`_
is a Python package for spectral bleed-through correction in microscopy stacks.
It provides linear two-channel unmixing workflows, optional bidirectional
correction, and PICASSO-family blind-unmixing workflows for multi-channel data.

The package is built around microscopy stacks in canonical ``TZCYX`` order and
uses `OMIO <https://omio.readthedocs.io/en/latest/>`_ for I/O. This means that
the workflows are not restricted to TIFF files alone, but can be applied to any
microscopy format currently supported by OMIO.

In addition to the core unmixing functions, the package includes optional
helpers for filtering, projection, time registration, intra-stack z-drift
correction, and napari-based inspection.


.. toctree::
   :maxdepth: 3
   :caption: Contents

   overview
   installation
   usage
   api
   changelog
   contributing
