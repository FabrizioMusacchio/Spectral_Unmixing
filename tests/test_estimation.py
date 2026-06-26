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
    def _make_random_signals(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        rng = np.random.default_rng(123)
        source = rng.gamma(shape=2.0, scale=10.0, size=(3, 24, 24)).astype(np.float32)
        target_true = rng.gamma(shape=1.8, scale=6.0, size=(3, 24, 24)).astype(np.float32)
        target_true[:, :8, :8] = 0.0
        target_measured = target_true + 0.28 * source
        return source, target_true, target_measured.astype(np.float32)

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

    def test_linear_fit_estimates_expected_alpha(self) -> None:
        rng = np.random.default_rng(3)
        source = rng.gamma(shape=2.0, scale=8.0, size=(3, 20, 20)).astype(np.float32)
        target_measured = 0.28 * source
        alpha = estimate_alpha_from_volume(
            source,
            target_measured,
            method="linear_fit",
            signal_percentile=70.0,
            background_percentile=0.0,
            min_mask_voxels=20,
        )
        self.assertAlmostEqual(alpha, 0.28, delta=0.06)

    def test_corr_min_reduces_correlation(self) -> None:
        source, _, target_measured = self._make_random_signals()
        alpha = estimate_alpha_from_volume(
            source,
            target_measured,
            method="corr_min",
            signal_percentile=60.0,
            background_percentile=0.0,
            alpha_max=1.0,
            min_mask_voxels=20,
        )
        corrected = target_measured - alpha * source
        corr_before = float(np.corrcoef(source.ravel(), target_measured.ravel())[0, 1])
        corr_after = float(np.corrcoef(source.ravel(), corrected.ravel())[0, 1])
        self.assertGreaterEqual(alpha, 0.0)
        self.assertLess(abs(corr_after), abs(corr_before))

    def test_mi_min_reduces_mutual_information(self) -> None:
        from spectral_unmixing import mutual_information_1d

        source, _, target_measured = self._make_random_signals()
        alpha = estimate_alpha_from_volume(
            source,
            target_measured,
            method="mi_min",
            signal_percentile=60.0,
            background_percentile=0.0,
            alpha_max=1.0,
            mi_bins=32,
            min_mask_voxels=20,
            max_alpha_voxels=None,
        )
        corrected = target_measured - alpha * source
        mi_before = mutual_information_1d(source.ravel(), target_measured.ravel(), bins=32)
        mi_after = mutual_information_1d(source.ravel(), corrected.ravel(), bins=32)
        self.assertGreaterEqual(alpha, 0.0)
        self.assertLess(mi_after, mi_before)

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
