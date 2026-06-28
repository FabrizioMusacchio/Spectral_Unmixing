Spectral Unmixing Documentation
===============================

.. figure:: _static/spectra_unmixing_logo.png
   :alt: Spectral Unmixing logo
   :align: center
   :figwidth: 55%


.. image:: https://badgen.net/badge/icon/GitHub%20repository?icon=github&label
   :target: https://github.com/FabrizioMusacchio/Spectral_Unmixing
   :alt: GitHub Repository

.. image:: https://img.shields.io/github/v/release/FabrizioMusacchio/Spectral_Unmixing
   :alt: GitHub Release

.. image:: https://img.shields.io/pypi/v/spectral-unmixing.svg
   :target: https://pypi.org/project/spectral-unmixing/
   :alt: PyPI version

.. image:: https://img.shields.io/badge/License-GPL%20v3-green.svg
   :target: https://github.com/FabrizioMusacchio/Spectral_Unmixing/blob/main/LICENSE
   :alt: GPLv3 License

.. image:: https://github.com/FabrizioMusacchio/Spectral_Unmixing/actions/workflows/spectral_unmixing_tests.yml/badge.svg
   :alt: Tests

.. image:: https://codecov.io/gh/FabrizioMusacchio/Spectral_Unmixing/graph/badge.svg?token=207WTB265P 
   :target: https://codecov.io/gh/FabrizioMusacchio/Spectral_Unmixing
   :alt: Codecov

.. image:: https://img.shields.io/github/last-commit/FabrizioMusacchio/Spectral_Unmixing
   :target: https://github.com/FabrizioMusacchio/Spectral_Unmixing/commits/main/
   :alt: GitHub last commit

.. image:: https://img.shields.io/github/issues/FabrizioMusacchio/Spectral_Unmixing
   :target: https://github.com/FabrizioMusacchio/Spectral_Unmixing/issues
   :alt: GitHub Issues Open

.. image:: https://img.shields.io/github/issues-closed/FabrizioMusacchio/Spectral_Unmixing?color=53c92e
   :target: https://github.com/FabrizioMusacchio/Spectral_Unmixing/issues?q=is%3Aissue%20state%3Aclosed
   :alt: GitHub Issues Closed

.. image:: https://img.shields.io/github/issues-pr/FabrizioMusacchio/Spectral_Unmixing
   :target: https://github.com/FabrizioMusacchio/Spectral_Unmixing/pulls
   :alt: GitHub Issues or Pull Requests

.. image:: https://readthedocs.org/projects/spectral-unmixing/badge/?version=latest
   :target: https://spectral-unmixing.readthedocs.io/en/latest/
   :alt: Documentation Status

.. image:: https://img.shields.io/github/languages/code-size/FabrizioMusacchio/Spectral_Unmixing
   :alt: GitHub code size in bytes

.. image:: https://img.shields.io/pypi/dm/spectral-unmixing?logo=pypy&label=PiPY%20downloads&color=blue
   :target: https://pypistats.org/packages/spectral-unmixing
   :alt: PyPI Downloads

.. image:: https://static.pepy.tech/personalized-badge/spectral-unmixing?period=total&units=INTERNATIONAL_SYSTEM&left_color=GRAY&right_color=BLUE&left_text=PiPY+total+downloads
   :target: https://pepy.tech/projects/spectral-unmixing
   :alt: PyPI Total Downloads

.. image:: https://img.shields.io/badge/Example%20Datasets-10.5281%2Fzenodo.20984021-blue
   :target: https://doi.org/10.5281/zenodo.20984021
   :alt: Spectral Unmixing example datasets on Zenodo

.. image:: https://img.shields.io/badge/Zenodo%20Archive-10.5281%2Fzenodo.20933784-blue
   :target: https://doi.org/10.5281/zenodo.20933784
   :alt: Zenodo Archive



`spectral-unmixing <https://github.com/FabrizioMusacchio/Spectral_Unmixing>`_
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
