"""
Unit tests for registration helpers.

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

from spectral_unmixing.registration import register_stack


class RegistrationTests(unittest.TestCase):
    def _build_shifted_stack(self) -> np.ndarray:
        reference = np.zeros((2, 2, 32, 32), dtype=np.float32)
        reference[:, :, 12:20, 11:19] = 10.0

        stack = np.zeros((3, 2, 2, 32, 32), dtype=np.float32)
        stack[0] = reference
        stack[1] = np.asarray(
            [
                [
                    ndi_shift(reference[z, c], shift=(2.0, -3.0), order=1, mode="constant", cval=0.0)
                    for c in range(reference.shape[1])
                ]
                for z in range(reference.shape[0])
            ],
            dtype=np.float32,
        )
        stack[2] = np.asarray(
            [
                [
                    ndi_shift(reference[z, c], shift=(-1.0, 4.0), order=1, mode="constant", cval=0.0)
                    for c in range(reference.shape[1])
                ]
                for z in range(reference.shape[0])
            ],
            dtype=np.float32,
        )
        return stack

    def test_phase_cross_correlation_registers_shifted_stack(self) -> None:
        stack = self._build_shifted_stack()
        registered = register_stack(
            stack,
            registration_channel=0,
            method="phase_cross_correlation",
            verbose=False,
        )

        np.testing.assert_allclose(registered[0], stack[0], atol=1e-5)
        self.assertLess(np.mean(np.abs(registered[1] - stack[0])), 0.25)
        self.assertLess(np.mean(np.abs(registered[2] - stack[0])), 0.25)

    def test_pystackreg_registers_shifted_stack(self) -> None:
        stack = self._build_shifted_stack()
        registered = register_stack(
            stack,
            registration_channel=0,
            method="pystackreg",
            verbose=False,
        )

        self.assertLess(np.mean(np.abs(registered[1] - stack[0])), 0.5)
        self.assertLess(np.mean(np.abs(registered[2] - stack[0])), 0.5)

    def test_registration_requires_multiple_timepoints(self) -> None:
        stack = np.zeros((1, 2, 2, 16, 16), dtype=np.float32)
        with self.assertRaises(ValueError):
            register_stack(stack, registration_channel=0, verbose=False)


if __name__ == "__main__":
    unittest.main()
