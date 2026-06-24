"""
Unit tests for filtering and projection helpers.

Author: Fabrizio Musacchio
Date: June 2026
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from spectral_unmixing.filters import apply_filters, max_z_project


class FilterTests(unittest.TestCase):
    def test_apply_filters_preserves_shape_for_singleton_z(self) -> None:
        stack = np.zeros((2, 1, 2, 5, 5), dtype=np.float32)
        stack[:, :, :, 2, 2] = 10.0

        filtered = apply_filters(
            stack,
            filters=["median", "gaussian"],
            median_size=3,
            gaussian_sigma=1.0,
            apply_3d=False,
        )

        self.assertEqual(filtered.shape, stack.shape)
        self.assertEqual(filtered.dtype, np.float32)

    def test_apply_filters_supports_2d_input(self) -> None:
        image = np.zeros((5, 5), dtype=np.float32)
        image[2, 2] = 10.0

        filtered = apply_filters(
            image,
            filters="gaussian",
            gaussian_sigma=1.0,
            apply_3d=False,
        )

        self.assertEqual(filtered.shape, image.shape)
        self.assertGreater(float(filtered[2, 2]), 0.0)

    def test_max_z_project_keeps_t_and_c_with_singleton_z(self) -> None:
        stack = np.zeros((2, 3, 2, 4, 4), dtype=np.float32)
        stack[:, 0, :, :, :] = 1.0
        stack[:, 1, :, :, :] = 3.0
        stack[:, 2, :, :, :] = 2.0

        projected = max_z_project(stack)

        self.assertEqual(projected.shape, (2, 1, 2, 4, 4))
        np.testing.assert_allclose(projected[:, 0, :, :, :], 3.0)


if __name__ == "__main__":
    unittest.main()
