# Gockel_Nieves_Rivera_2026
This folder contains a cropped example stack used in the helper tutorials of `spectral-unmixing`, in particular for registration, filtering, histogram matching, and projection workflows.

Please download the data from Zenodo and place it in this folder before running the tutorials or tests:

Musacchio, F. (2026). Example datasets for the `spectral-unmixing` pipeline [Data set]. Zenodo. <https://doi.org/10.5281/zenodo.20984022>

## Contents

- `Gockel_Nieves_Rivera_2026_5D_stack.tif`
  Cropped cut-out used as the main tutorial input stack.
- `registered/`
  Derived example outputs generated from the cropped stack with the helper
  scripts in this repository.

## Source
The input stack in this folder is a cropped cut-out derived from:

Gockel, N., Nieves-Rivera, N., Musacchio, F., Druart, M., Jaako, K., Fuhrmann, F., Rozkalne, R., Poll, S., Baiba, J., Fuhrmann, M., & Le Magueresse, C. (2025). *Example Datasets for Microglial Motility Analysis Using the MotilA Pipeline* [Data set]. Zenodo. <https://doi.org/10.5281/zenodo.15061566>

## License
The source dataset is licensed under:

- Creative Commons Attribution Share Alike 4.0 International
- https://creativecommons.org/licenses/by-sa/4.0/

Because the stack in this folder is a cropped derivative of that dataset, it is distributed here under the same CC BY-SA 4.0 license.

The derived files in `registered/` were generated from this cropped stack with the scripts in this repository and are likewise distributed under CC BY-SA 4.0.

## Notes

- This folder is included as a tutorial example dataset for the
  `spectral-unmixing` package.
- The files here are intended to help users reproduce the helper tutorials and
  test the registration, filtering, histogram-matching, and projection
  workflow on a real microscopy stack.
