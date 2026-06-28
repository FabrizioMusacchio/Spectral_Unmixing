Installation
============

We recommend creating a dedicated Python 3.12 conda environment first:

.. code-block:: bash

   conda create -n spectral-unmixing python=3.12
   conda activate spectral-unmixing

.. note::

   We have tested ``spectral-unmixing`` with Python 3.12. Newer versions may work but 
   are not guaranteed to be compatible. Older Python versions are not supported,
   as ``spectral-unmixing`` relies on `OMIO <https://omio.readthedocs.io/en/latest/>`_ for 
   reading microscopy data, and OMIO requires Python 3.12 or newer.



Install from PyPI
-----------------

The standard installation path is:

.. code-block:: bash

   pip install spectral-unmixing


Upgrade an existing installation
--------------------------------

To upgrade an existing PyPI installation:

.. code-block:: bash

   pip install --upgrade spectral-unmixing


Developer installation
----------------------

To work on the repository locally:

.. code-block:: bash

   git clone https://github.com/FabrizioMusacchio/spectral-unmixing.git
   cd spectral-unmixing
   pip install -e .

To develop the package *and* build the documentation locally:

.. code-block:: bash

   pip install -e .[dev,docs]


Runtime dependencies
--------------------

Main runtime dependencies currently include:

- ``numpy``
- ``omio-microscopy``
- ``scipy``
- ``scikit-image``
- ``pystackreg``
- ``ipykernel``


Documentation build dependencies
--------------------------------

If you want to build the Sphinx documentation locally, install the package
with the ``docs`` extra:

.. code-block:: bash

   pip install -e .[docs]


Quick verification
------------------

After installation, you can verify that the package imports correctly:

.. code-block:: bash

   python -c "import spectral_unmixing; print(spectral_unmixing.__all__)"
