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

from spectral_unmixing.filters import (
    apply_filters,
    match_histograms_across_time,
    max_z_project,
)


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

    def test_max_z_project_uses_requested_zrange(self) -> None:
        stack = np.zeros((1, 4, 1, 3, 3), dtype=np.float32)
        stack[:, 0, :, :, :] = 1.0
        stack[:, 1, :, :, :] = 5.0
        stack[:, 2, :, :, :] = 2.0
        stack[:, 3, :, :, :] = 9.0

        projected = max_z_project(stack, zrange=(1, 3))

        self.assertEqual(projected.shape, (1, 1, 1, 3, 3))
        np.testing.assert_allclose(projected[:, 0, :, :, :], 5.0)

    def test_max_z_project_clamps_out_of_bounds_zrange(self) -> None:
        stack = np.zeros((1, 3, 1, 2, 2), dtype=np.float32)
        stack[:, 0, :, :, :] = 2.0
        stack[:, 1, :, :, :] = 4.0
        stack[:, 2, :, :, :] = 6.0

        projected = max_z_project(stack, zrange=(-10, 99))

        np.testing.assert_allclose(projected[:, 0, :, :, :], 6.0)

    def test_match_histograms_across_time_preserves_shape(self) -> None:
        stack = np.zeros((3, 2, 2, 8, 8), dtype=np.float32)
        stack[0, :, 0, 2:6, 2:6] = 10.0
        stack[1, :, 0, 2:6, 2:6] = 30.0
        stack[2, :, 0, 2:6, 2:6] = 50.0
        stack[:, :, 1, 1:7, 1:7] = np.array([5.0, 20.0, 40.0], dtype=np.float32)[:, None, None, None]

        matched = match_histograms_across_time(stack, reference_t=0)

        self.assertEqual(matched.shape, stack.shape)
        self.assertEqual(matched.dtype, np.float32)
        np.testing.assert_allclose(matched[0], stack[0], atol=1e-6)

    def test_match_histograms_across_time_requires_multiple_timepoints(self) -> None:
        stack = np.zeros((1, 2, 2, 8, 8), dtype=np.float32)
        with self.assertRaises(ValueError):
            match_histograms_across_time(stack)


if __name__ == "__main__":
    unittest.main()
