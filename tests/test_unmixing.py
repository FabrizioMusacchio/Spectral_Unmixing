"""
Unit tests for stack-level bleed-through correction.

Author: Fabrizio Musacchio
Date: June 2026
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from spectral_unmixing.unmixing import report_path_from_output_path, unmix


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
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.tif"

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
                returned_output = unmix(
                    input_path="input.tif",
                    output_path=output_path,
                    alpha=0.2,
                    alpha_mode="fixed",
                    verbose=False,
                )

            np.testing.assert_allclose(
                written["stack"][:, :, 0, :, :], self.stack[:, :, 0, :, :]
            )
            np.testing.assert_allclose(written["stack"][:, :, 1, :, :], 3.0)
            self.assertEqual(returned_output, output_path)

            report_path = report_path_from_output_path(output_path)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["alpha"], 0.2)
            self.assertEqual(report["alpha_mode"], "fixed")
            self.assertEqual(report["size_t"], 2)
            self.assertEqual(report["size_z"], 3)

    def test_per_t_reports_one_alpha_per_timepoint(self) -> None:
        written = {}
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.tif"

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
                unmix(
                    input_path="input.tif",
                    output_path=output_path,
                    alpha_mode="per_t",
                    signal_percentile=50.0,
                    background_percentile=0.0,
                    verbose=False,
                )

            np.testing.assert_allclose(written["stack"][0, :, 1, :, :], 3.0, atol=1e-6)
            np.testing.assert_allclose(written["stack"][1, :, 1, :, :], 5.0, atol=1e-6)

            report_path = report_path_from_output_path(output_path)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(len(report["alpha_values"]), 2)

    def test_refuses_to_overwrite_input_path(self) -> None:
        with self.assertRaises(ValueError):
            unmix(
                input_path="same.tif",
                output_path="same.tif",
                alpha=0.2,
                alpha_mode="fixed",
            )

    def test_singleton_t_and_z_stack_is_supported(self) -> None:
        singleton_stack = np.zeros((1, 1, 2, 4, 5), dtype=np.float32)
        singleton_stack[0, 0, 0, :, :] = np.array(
            [
                [4.0, 5.0, 6.0, 7.0, 8.0],
                [5.0, 6.0, 7.0, 8.0, 9.0],
                [6.0, 7.0, 8.0, 9.0, 10.0],
                [7.0, 8.0, 9.0, 10.0, 11.0],
            ],
            dtype=np.float32,
        )
        singleton_stack[0, 0, 1, :, :] = 2.0 + 0.25 * singleton_stack[0, 0, 0, :, :]
        metadata = {"axes": "TZCYX", "shape": singleton_stack.shape}
        written = {}

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "singleton_output.tif"

            def fake_load_stack_with_omio(_input_path):
                return singleton_stack.copy(), metadata.copy()

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
                unmix(
                    input_path="input.tif",
                    output_path=output_path,
                    alpha=0.25,
                    alpha_mode="fixed",
                    verbose=False,
                )

            np.testing.assert_allclose(written["stack"][0, 0, 1, :, :], 2.0, atol=1e-6)
            report = json.loads(
                report_path_from_output_path(output_path).read_text(encoding="utf-8")
            )
            self.assertEqual(report["size_t"], 1)
            self.assertEqual(report["size_z"], 1)
            self.assertFalse(report["has_multiple_t"])
            self.assertFalse(report["has_multiple_z"])


if __name__ == "__main__":
    unittest.main()
