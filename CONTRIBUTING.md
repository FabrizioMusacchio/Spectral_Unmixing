# How to Contribute

Thank you for your interest in contributing to Spectral Unmixing. This project
welcomes improvements to the code, documentation, tests, tutorials, and overall
usability.

Spectral Unmixing is a Python package for reproducible spectral bleed-through
correction in multidimensional microscopy stacks. Its core focus is directed
two-channel unmixing, optional bidirectional two-channel correction, and
PICASSO-family blind unmixing, with additional helper modules for filtering,
projection, registration, viewing, and OMIO-based image I/O.

## Before you start

Please check the GitHub issue tracker to see whether your idea, bug report, or
enhancement has already been discussed:

[https://github.com/FabrizioMusacchio/Spectral_Unmixing/issues](https://github.com/FabrizioMusacchio/Spectral_Unmixing/issues)

If a related issue exists, comment there to indicate your interest or to add
relevant technical details. If no issue exists, open a new one with a short
description of:

- what you would like to change or add;
- why it is useful in the context of spectral unmixing or microscopy workflows;
- any thoughts on implementation, edge cases, or testing.

For small fixes such as typos, broken links, or minor documentation
improvements, opening a pull request directly is fine.

## Development Environment

Spectral Unmixing requires Python 3.12 or newer. It builds on standard
scientific Python packages used in bioimage analysis, including NumPy, SciPy,
scikit-image, pystackreg, and OMIO for microscopy image I/O.

A typical development setup using `conda` looks like this:

```sh
git clone https://github.com/FabrizioMusacchio/Spectral_Unmixing.git
cd Spectral_Unmixing

conda create -n spectral-unmixing-dev -c conda-forge python=3.12
conda activate spectral-unmixing-dev

pip install -e .
```

To install optional development dependencies such as testing and coverage tools:

```sh
pip install -e ".[dev]"
```

Documentation dependencies can be installed with:

```sh
pip install -e ".[docs]"
```

## Making Changes and Opening Pull Requests

All code contributions should be submitted as pull requests against the `main`
branch of the repository.

A recommended workflow:

1. Create a new feature branch:

   ```sh
   git checkout -b feature/my-feature
   ```

2. Implement your changes. New public functions or modules should include clear
   NumPy-style docstrings explaining purpose, inputs, outputs, assumptions, and
   limitations.

3. Add or update tests where appropriate.

4. Push your branch and open a pull request with a concise title, a short
   explanation of what changed and why, and references to related issues.

Draft pull requests are welcome if you would like feedback during development.

## Commit Conventions

Clear and consistent commit messages help keep the project history readable.
Prefixes inspired by Conventional Commits are encouraged:

- `feat:` new functionality
- `fix:` bug fixes
- `docs:` documentation changes
- `refactor:` internal code restructuring without behavior changes
- `test:` adding or modifying tests
- `chore:` maintenance tasks or tooling updates

Example:

```text
fix: preserve target dtype metadata in unmixing report
```

## Testing

Spectral Unmixing uses `pytest` for automated testing. To run the full test
suite locally:

```sh
pytest
```

If you add new functionality or fix a bug, please extend the test suite where
the behavior can be tested reliably.

Good tests for this project usually use synthetic arrays rather than large
microscopy files. This is especially useful for:

- known bleed-through coefficients;
- known bidirectional mixing matrices;
- expected alpha-estimation behavior;
- expected output shapes for `T=1`, `Z=1`, and full `TZCYX` stacks;
- registration or filtering behavior on small synthetic images.

Large microscopy datasets should not be added to the repository. If a real
example is necessary to reproduce a bug, provide the smallest cropped or
anonymized file that still demonstrates the issue, preferably via a public or
temporary external link.

## Reproducibility and Scientific Contribution Guidelines

Spectral Unmixing is intended to support reproducible scientific image analysis.
Contributions should therefore respect the following principles:

- **Deterministic behavior:** Given identical inputs and parameters, results
  should be deterministic whenever possible. Randomized algorithms should expose
  a `random_state` or equivalent parameter.

- **Transparent parameter handling:** Defaults, fallback behavior, clipping,
  masking, normalization, and background handling should be explicit and
  documented.

- **Sidecar-report consistency:** User-facing processing functions that write
  output data should also preserve enough processing information for the result
  to be audited later.

- **Documented limitations:** Known limitations or unsupported cases should be
  documented rather than silently ignored.

- **Minimal scope changes:** Pull requests should focus on a well-defined
  improvement. Large conceptual changes should be discussed in an issue before
  implementation.

These guidelines help keep the package reviewable, maintainable, and suitable
for long-term archival publication.

## Spectral-Unmixing Design Constraints

Spectral Unmixing makes several explicit design choices. Contributions that
change these choices should be discussed in an issue before implementation,
because they affect reproducibility and user expectations.

- **Canonical stack model:** Core image-processing functions operate on
  canonical `TZCYX` stacks. Simpler inputs such as `T=1` or `Z=1` should remain
  supported, but ambiguous axis handling should not be hidden inside numerical
  algorithms.

- **Separation of alpha mode and alpha method:** `alpha_mode` controls where an
  alpha value is obtained from, while `method` controls how alpha is estimated.
  New alpha-estimation routines should preserve this separation.

- **Original data correction:** Background subtraction, masking, and
  preprocessing may be used for alpha estimation, but the final correction
  should remain mathematically explicit and documented.

- **Conservative blind-unmixing claims:** PICASSO-family workflows should be
  described carefully. New blind-unmixing methods should document assumptions,
  optimization criteria, ambiguity handling, and expected failure modes.

- **OMIO as I/O layer:** File-format support is delegated to OMIO. Spectral
  Unmixing should not grow its own independent microscopy file readers unless
  there is a strong, discussed reason.

## Reporting Bugs

Please report bugs via the GitHub issue tracker:

[https://github.com/FabrizioMusacchio/Spectral_Unmixing/issues](https://github.com/FabrizioMusacchio/Spectral_Unmixing/issues)

Include the following information if possible:

- Spectral Unmixing version, for example from `pip show spectral-unmixing`;
- Python version;
- operating system;
- installation method, such as PyPI install, editable install, or local checkout;
- minimal code snippet or steps to reproduce the issue;
- full traceback or warning output;
- input stack shape and expected axis order;
- alpha mode, alpha method, source channel, target channel, and other relevant
  parameters;
- if applicable, the generated JSON sidecar report;
- if applicable, a small synthetic or cropped example file illustrating the
  problem.

If the bug appears to be caused by reading or writing a specific microscopy file
format, please include the OMIO version and, if possible, a minimal example file.
Format-reader bugs may ultimately need to be reported upstream to OMIO, but an
issue in Spectral Unmixing is still useful if the problem appears inside an
unmixing workflow.

## Requests for New Methods and Workflow Extensions

Users are encouraged to request or contribute new spectral-unmixing methods,
alpha estimators, blind-unmixing variants, or pre- and post-processing helpers.

Useful feature requests include:

- a short description of the method or workflow;
- the scientific use case it addresses;
- whether it applies to directed two-channel, bidirectional, or multi-channel
  blind unmixing;
- expected input and output shapes;
- key parameters and reasonable defaults;
- references to papers, algorithms, or existing implementations;
- ideas for synthetic tests or example data.

For new alpha-estimation methods, please describe how alpha is estimated, which
mask or preprocessing assumptions are required, and how failure cases should be
reported. For new PICASSO-family or blind-unmixing methods, please describe the
optimization objective and how channel permutation or scale ambiguity should be
handled.

## Documentation and Tutorials

Documentation changes are very welcome. Good documentation contributions include:

- fixing outdated parameter descriptions;
- adding short conceptual explanations for methods;
- improving tutorial text or comments in `user_scripts`;
- adding examples for edge cases such as `T=1`, `Z=1`, or per-time-point alpha;
- clarifying how JSON sidecar reports support reproducibility.

When updating tutorials, please keep examples reproducible and avoid relying on
private datasets.

## License and Contributions

By submitting a pull request, you agree that your contributions will be released
under the project's license as specified in the repository.

If you are unsure how to begin or would like to discuss a potential
contribution, open an issue to start a conversation.
