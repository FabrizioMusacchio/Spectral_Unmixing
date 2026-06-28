Example datasets
================

The `spectral-unmixing` tutorials are backed by a small public example-dataset
collection that is available as a Zenodo dataset release:

   Musacchio, F. (2026). *Example datasets for the spectral-unmixing pipeline*
   [Data set]. Zenodo. https://doi.org/10.5281/zenodo.20984021

The goal of this collection is to let users replay the tutorials with the same
inputs that are used throughout the documentation and the interactive scripts
in ``user_scripts/``.


What is included
----------------

The example data are organized into three main folders:

- ``example_data/PICASSO_examples/``
  public two-channel, three-channel, and five-channel PICASSO-family example
  images used by the blind-unmixing tutorials
- ``example_data/synthetic_data/``
  synthetic `TZCYX` data with controlled two-channel bleed-through, used for
  the full-stack unmixing tutorial
- ``example_data/Gockel_Nieves_Rivera_2026/``
  a cropped real microscopy stack used for the helper tutorial on
  registration, filtering, histogram matching, and projection

Each folder contains a dedicated README file with more details about the dataset, 
its provenance, and the license under which it is redistributed. The entire collection 
is licensed under CC BY 4.0.


Dataset-to-tutorial mapping
---------------------------

The main public tutorial mappings are:

- :doc:`usage_unmix_example`
  uses ``example_data/PICASSO_examples/2_color_unmixing_validation.tif``
- :doc:`usage_unmix_full_tzcyx_synthetic_example`
  uses ``example_data/synthetic_data/synthetic_bleedthrough_T9_Z20_C2.tif``
- :doc:`usage_unmix_bidirectional_example`
  uses ``example_data/PICASSO_examples/bidirectional_example.tif``
- :doc:`usage_unmix_picasso_2color_example`
  uses ``example_data/PICASSO_examples/2_color_unmixing_validation.tif``
- :doc:`usage_unmix_picasso_3color_example`
  uses ``example_data/PICASSO_examples/3_color_data.tif``
- :doc:`usage_unmix_picasso_5color_example`
  uses ``example_data/PICASSO_examples/5_color_unmixing_simulation.tif``
- :doc:`usage_filter_and_register_stack`
  uses ``example_data/Gockel_Nieves_Rivera_2026/Gockel_Nieves_Rivera_2026_5D_stack.tif``

Dataset descriptions and provenance
-----------------------------------

``PICASSO_examples``
~~~~~~~~~~~~~~~~~~~~

This folder contains a subset of the public example images shared with the
PICASSO publication.

It includes public two-channel, three-channel, and five-channel example stacks
that are suitable for:

- standard two-channel unmixing examples,
- bidirectional unmixing examples,
- and PICASSO-family blind-unmixing tutorials.

Source:

   Chang, Jae-Byum; Seo, Junyoung; Sim, Yeonbo; Kim, Jeewon; Kim, Hyunwoo;
   Cho, In; et al. (2022). *PICASSO allows ultra-multiplexed fluorescence
   imaging of spatially overlapping proteins without reference spectra
   measurements*. figshare. Figure.
   https://doi.org/10.6084/m9.figshare.19596682.v1

Associated paper:

   Seo, J., Sim, Y., Kim, J. et al. *PICASSO allows ultra-multiplexed
   fluorescence imaging of spatially overlapping proteins without reference
   spectra measurements*. Nature Communications 13, 2475 (2022).
   https://doi.org/10.1038/s41467-022-30168-z

License:

- CC BY 4.0
- https://creativecommons.org/licenses/by/4.0/


``synthetic_data``
~~~~~~~~~~~~~~~~~~

This folder contains synthetic tutorial data created specifically for
`spectral-unmixing`.

The main stack is:

- ``synthetic_bleedthrough_T9_Z20_C2.tif``

It is a synthetic `TZCYX` stack with:

- ``T=9``
- ``Z=20``
- ``C=2``

The stack contains two time-varying 3D Gaussian structures and controlled
bleed-through from channel 0 into channel 1. It is especially useful for the
full-stack tutorial because the construction is known and the behavior of the
unmixing methods can be tested in a controlled setting.

Source:

- generated within this repository using
  ``additional_scripts/generate_synthetic_bleedthrough_stack.py``

License:

- CC BY 4.0
- https://creativecommons.org/licenses/by/4.0/


``Gockel_Nieves_Rivera_2026``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This folder contains a cropped real microscopy stack used in the helper
tutorial for:

- intra-stack z-drift correction,
- across-time registration,
- histogram matching across time,
- filtering,
- and max-z projection.

Source:

   Gockel, N., Nieves-Rivera, N., Musacchio, F., Druart, M., Jaako, K.,
   Fuhrmann, F., Rozkalne, R., Poll, S., Baiba, J., Fuhrmann, M., &
   Le Magueresse, C. (2025). *Example Datasets for Microglial Motility
   Analysis Using the MotilA Pipeline* [Data set]. Zenodo.
   https://doi.org/10.5281/zenodo.15061566

The file used here is a cropped derivative of that public source dataset.

License:

- CC BY-SA 4.0
- https://creativecommons.org/licenses/by-sa/4.0/

Because the stack in this folder is a cropped derivative, it is redistributed
under the same CC BY-SA 4.0 license.


Citation
--------

If you are using the example datasets in your own work, please cite the Zenodo dataset release:

   Musacchio, F. (2026). *Example datasets for the spectral-unmixing pipeline*
   [Data set]. Zenodo. https://doi.org/10.5281/zenodo.20984021  