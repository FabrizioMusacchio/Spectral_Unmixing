"""
Unit tests for intra-stack Z-drift correction helpers.

Author: Fabrizio Musacchio
Date: June 2026
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.ndimage import shift as ndi_shift

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from spectral_unmixing.registration import correct_intra_stack_z_drift


class IntraStackRegistrationTests(unittest.TestCase):
    def _base_plane(self) -> np.ndarray:
        plane = np.zeros((48, 48), dtype=np.float32)
        plane[12:22, 14:26] = 10.0
        plane[28:35, 30:40] = 6.0
        plane[18:22, 34:38] = 12.0
        return plane

    def _build_stack_with_slice_shifts(self) -> tuple[np.ndarray, np.ndarray]:
        base = self._base_plane()
        reference = np.zeros((2, 5, 2, 48, 48), dtype=np.float32)

        for t in range(reference.shape[0]):
            for z in range(reference.shape[1]):
                reference[t, z, 0, :, :] = 0.5 * base
                reference[t, z, 1, :, :] = base

        shifted = reference.copy()
        shifts = {
            (0, 2): (2.0, -3.0),
            (1, 2): (2.0, -3.0),
        }

        for (t, z), shift_yx in shifts.items():
            for c in range(shifted.shape[2]):
                shifted[t, z, c, :, :] = ndi_shift(
                    shifted[t, z, c, :, :],
                    shift=shift_yx,
                    order=1,
                    mode="constant",
                    cval=0.0,
                )

        return reference, shifted

    def test_neighbor_mode_corrects_shifted_slices_using_channel_one(self) -> None:
        reference, shifted = self._build_stack_with_slice_shifts()

        corrected = correct_intra_stack_z_drift(
            shifted,
            registration_channel=1,
            method="phase_cross_correlation",
            reference_mode="neighbor",
            neighbor_window_size=3,
            verbose=False,
        )

        before_error = np.mean(np.abs(shifted[:, :, 1, :, :] - reference[:, :, 1, :, :]))
        after_error = np.mean(np.abs(corrected[:, :, 1, :, :] - reference[:, :, 1, :, :]))

        self.assertLess(after_error, before_error)
        self.assertLess(after_error, 0.35)
        self.assertEqual(corrected.shape, shifted.shape)

    def test_singleton_z_returns_copy_without_error(self) -> None:
        stack = np.random.default_rng(42).random((2, 1, 2, 16, 16), dtype=np.float32)
        corrected = correct_intra_stack_z_drift(
            stack,
            registration_channel=1,
            verbose=False,
        )

        np.testing.assert_allclose(corrected, stack)

    def test_even_neighbor_window_size_raises(self) -> None:
        stack = np.random.default_rng(123).random((1, 3, 2, 16, 16), dtype=np.float32)
        with self.assertRaises(ValueError):
            correct_intra_stack_z_drift(
                stack,
                registration_channel=0,
                neighbor_window_size=4,
                verbose=False,
            )


if __name__ == "__main__":
    unittest.main()
