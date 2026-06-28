# synthetic_data

This folder contains synthetic example data generated specifically for the `spectral-unmixing` tutorials and tests.

## Contents

- `synthetic_bleedthrough_T9_Z20_C2.tif`
  Synthetic two-channel `TZCYX` microscopy stack with controlled channel-0 to channel-1 bleed-through.
- `synthetic_bleedthrough_T9_Z20_C2.tif.json`
  JSON sidecar generated together with the synthetic stack.


## Source
The synthetic input stack in this folder was generated within this repository using:

- `additional_scripts/generate_synthetic_bleedthrough_stack.py`

The generator script creates a `T=9`, `Z=20`, `C=2` stack with two time-varying 3D Gaussian spheres and controlled bleed-through from channel 0 into channel 1.

## License
The synthetic dataset and the derived tutorial outputs in this folder are distributed under:

- CC BY 4.0
- https://creativecommons.org/licenses/by/4.0/

## Notes

- This folder is included as a public tutorial example dataset for the `spectral-unmixing` package.
- It is intended to let users reproduce the full `TZCYX` unmixing tutorial on data with known synthetic ground truth and controlled bleed-through.
