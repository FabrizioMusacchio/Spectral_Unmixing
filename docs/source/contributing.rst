Contributing
============

Contributions are welcome.

Ways to contribute
------------------

You can contribute by:

- reporting bugs,
- suggesting new unmixing or filtering features,
- improving documentation,
- adding tests,
- or contributing new public example workflows.


Development workflow
--------------------

A typical local development workflow is:

1. create a dedicated Python 3.12 environment,
2. install the package in editable mode,
3. run the relevant user scripts or unit tests,
4. make focused changes,
5. update documentation and changelog entries when appropriate.


Testing
-------

Before proposing changes, please run the unit tests:

.. code-block:: bash

   python -m unittest discover -s tests -v


Documentation contributions
---------------------------

Documentation improvements are highly encouraged. In particular, the project
uses:

- public docstrings intended for autodoc/API pages,
- tutorial-style user scripts in ``user_scripts/``,
- and Sphinx pages in ``docs/source/``.

If you add or substantially change user-facing behavior, please update both the
relevant docstring and the corresponding tutorial or overview page where
appropriate.


Issues and pull requests
------------------------

Please use the GitHub issue tracker for bug reports, questions, and feature
requests:

`https://github.com/FabrizioMusacchio/spectral-unmixing/issues <https://github.com/FabrizioMusacchio/spectral-unmixing/issues>`_
