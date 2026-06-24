"""
Unit tests for stack-level bleed-through correction.

Author: Fabrizio Musacchio
Date: June 2026
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from spectral_unmixing.unmixing import unmix_ch0_from_ch1


class UnmixingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.stack = np.zeros((2, 3, 2, 4, 5), dtype=np.float32)
        base_pattern = np.array(
            [
                [6.0, 7.0, 8.0, 9.0, 10.0],
                [7.0, 8.0, 9.0, 10.0, 11.0],
                [8.0, 9.0, 10.0, 11.0, 12.0],
                [9.0, 10.0, 11.0, 12.0, 13.0],
            ],
            dtype=np.float32,
        )
        self.stack[:, :, 0, :, :] = base_pattern
        self.stack[:, :, 1, :, :] = 3.0 + 0.2 * self.stack[:, :, 0, :, :]
        self.metadata = {"axes": "TZCYX", "shape": self.stack.shape}

    def test_fixed_alpha_corrects_only_target_channel(self) -> None:
        written = {}

        def fake_load_stack_with_omio(_input_path):
            return self.stack.copy(), self.metadata.copy()

        def fake_write_stack_with_omio(output_path, stack, metadata):
            written["output_path"] = str(output_path)
            written["stack"] = stack.copy()
            written["metadata"] = metadata.copy()
            return Path(output_path)

        with (
            patch(
                "spectral_unmixing.unmixing.load_stack_with_omio",
                side_effect=fake_load_stack_with_omio,
            ),
            patch(
                "spectral_unmixing.unmixing.write_stack_with_omio",
                side_effect=fake_write_stack_with_omio,
            ),
        ):
            report = unmix_ch0_from_ch1(
                input_path="input.tif",
                output_path="output.tif",
                alpha=0.2,
                alpha_mode="fixed",
            )

        np.testing.assert_allclose(
            written["stack"][:, :, 0, :, :], self.stack[:, :, 0, :, :]
        )
        np.testing.assert_allclose(written["stack"][:, :, 1, :, :], 3.0)
        self.assertEqual(report["alpha"], 0.2)
        self.assertEqual(report["alpha_mode"], "fixed")

    def test_per_t_reports_one_alpha_per_timepoint(self) -> None:
        written = {}

        varying_stack = self.stack.copy()
        varying_stack[1, :, 1, :, :] = 5.0 + 0.4 * varying_stack[1, :, 0, :, :]

        def fake_load_stack_with_omio(_input_path):
            return varying_stack.copy(), self.metadata.copy()

        def fake_write_stack_with_omio(output_path, stack, metadata):
            written["stack"] = stack.copy()
            return Path(output_path)

        with (
            patch(
                "spectral_unmixing.unmixing.load_stack_with_omio",
                side_effect=fake_load_stack_with_omio,
            ),
            patch(
                "spectral_unmixing.unmixing.write_stack_with_omio",
                side_effect=fake_write_stack_with_omio,
            ),
        ):
            report = unmix_ch0_from_ch1(
                input_path="input.tif",
                output_path="output.tif",
                alpha_mode="per_t",
                signal_percentile=50.0,
                background_percentile=0.0,
            )

        np.testing.assert_allclose(written["stack"][0, :, 1, :, :], 3.0, atol=1e-6)
        np.testing.assert_allclose(written["stack"][1, :, 1, :, :], 5.0, atol=1e-6)
        self.assertEqual(len(report["alpha_values"]), 2)

    def test_refuses_to_overwrite_input_path(self) -> None:
        with self.assertRaises(ValueError):
            unmix_ch0_from_ch1(
                input_path="same.tif",
                output_path="same.tif",
                alpha=0.2,
                alpha_mode="fixed",
            )


if __name__ == "__main__":
    unittest.main()
