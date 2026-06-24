# Simple spectral-unmixing for microscopy image stacks

`spectral-unmixing` provides a small, reusable Python package for bleed-through
correction in multi-dimensional microscopy stacks that follow OMIO's canonical
axis order `TZCYX`.

## Features

- OMIO-based reading and writing of microscopy TIFF stacks
- Validation that incoming data are in canonical `TZCYX` order
- Spectral bleed-through correction from one source channel into one target channel
- Three alpha modes:
  - fixed alpha
  - alpha estimated once from a reference time point
  - alpha estimated independently for each time point
- Small command-line interface for batch-friendly use
- Cell-structured user script for VS Code Interactive Window workflows

## Installation

```bash
pip install -e .
```

This project depends on:

- `numpy`
- `omio-microscopy`

## Quick Start

```python
from spectral_unmixing import unmix

output_path = unmix(
    input_path="input.tif",
    output_path="unmixed/input_unmixed_reference_t0.tif",
    alpha_mode="reference_t",
    alpha_reference_t=0,
    signal_percentile=99.0,
)
```

Each run writes a JSON sidecar report next to the output TIFF, for example
`input_unmixed_reference_t0.tif.json`, so the exact settings remain reproducible.
Terminal progress output is enabled by default and can be disabled with
`verbose=False`.

## Scientific Note

A fixed alpha measured from a proper single-label control recording is
scientifically preferable. Estimating alpha from the mixed experimental stack is
provided as a pragmatic first-pass workflow and can be biased when biological
signals overlap spatially.
