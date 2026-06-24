"""
Unit tests for publication-render helpers.

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

from spectral_unmixing.publication import render_for_publication, subtract_background


class PublicationRenderTests(unittest.TestCase):
    def test_subtract_background_preserves_shape_and_clips_to_zero(self) -> None:
        stack = np.full((2, 1, 2, 16, 16), 5.0, dtype=np.float32)
        stack[:, :, 0, 6:10, 6:10] = 20.0

        corrected = subtract_background(
            stack,
            method="gaussian",
            gaussian_sigma=(2.0, 2.0),
        )

        self.assertEqual(corrected.shape, stack.shape)
        self.assertEqual(corrected.dtype, np.float32)
        self.assertGreater(float(np.max(corrected[:, :, 0, :, :])), 0.0)
        self.assertGreaterEqual(float(np.min(corrected)), 0.0)

    def test_render_for_publication_returns_zero_to_one_stack(self) -> None:
        stack = np.zeros((3, 1, 2, 24, 24), dtype=np.float32)
        stack[:, :, 0, 8:16, 10:14] = np.array([10.0, 14.0, 20.0], dtype=np.float32)[:, None, None, None]
        stack[:, :, 1, 7:18, 7:18] = np.array([5.0, 7.0, 12.0], dtype=np.float32)[:, None, None, None]

        rendered = render_for_publication(
            stack,
            background_method="gaussian",
            gaussian_sigma=(2.0, 2.0),
            denoise_method="median",
            median_size=(3, 3),
            apply_unsharp_mask=True,
            unsharp_radius=(0.8, 0.8),
            unsharp_amount=(0.6, 0.6),
            lower_percentile=(0.0, 0.0),
            upper_percentile=(99.5, 99.5),
            gamma=(0.9, 0.9),
        )

        self.assertEqual(rendered.shape, stack.shape)
        self.assertEqual(rendered.dtype, np.float32)
        self.assertGreaterEqual(float(np.min(rendered)), 0.0)
        self.assertLessEqual(float(np.max(rendered)), 1.0)

    def test_render_for_publication_accepts_scalar_parameters(self) -> None:
        stack = np.zeros((1, 1, 2, 12, 12), dtype=np.float32)
        stack[:, :, :, 4:8, 4:8] = 10.0

        rendered = render_for_publication(
            stack,
            background_method="white_tophat",
            white_tophat_radius=3,
            denoise_method="none",
            apply_unsharp_mask=False,
            lower_percentile=0.0,
            upper_percentile=100.0,
            gamma=1.0,
        )

        self.assertEqual(rendered.shape, stack.shape)
        self.assertAlmostEqual(float(np.max(rendered)), 1.0, places=5)


if __name__ == "__main__":
    unittest.main()
