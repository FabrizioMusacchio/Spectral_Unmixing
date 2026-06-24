"""
Unit tests for alpha estimation helpers.

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

from spectral_unmixing import estimate_alpha_from_volume


class EstimateAlphaTests(unittest.TestCase):
    def test_estimate_alpha_matches_expected_ratio(self) -> None:
        source = np.zeros((2, 6, 6), dtype=np.float32)
        source[:, 2:5, 2:5] = np.array(
            [
                [80.0, 90.0, 100.0],
                [85.0, 95.0, 105.0],
                [88.0, 98.0, 108.0],
            ],
            dtype=np.float32,
        )
        target = 0.2 * source

        alpha = estimate_alpha_from_volume(
            source,
            target,
            signal_percentile=90.0,
            background_percentile=0.0,
            min_mask_voxels=4,
        )

        self.assertAlmostEqual(alpha, 0.2, places=6)

    def test_estimate_alpha_rejects_shape_mismatch(self) -> None:
        with self.assertRaises(ValueError):
            estimate_alpha_from_volume(np.zeros((2, 2)), np.zeros((2, 3)))

    def test_estimate_alpha_rejects_empty_signal_mask(self) -> None:
        source = np.zeros((2, 4, 4), dtype=np.float32)
        target = np.zeros_like(source)

        with self.assertRaises(ValueError):
            estimate_alpha_from_volume(
                source,
                target,
                signal_percentile=99.0,
                background_percentile=0.0,
            )


if __name__ == "__main__":
    unittest.main()
