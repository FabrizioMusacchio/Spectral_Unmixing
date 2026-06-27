"""
Unit tests for stack-level bleed-through correction.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
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

from spectral_unmixing.unmixing import report_path_from_output_path, unmix, unmix_picasso

# %% CLASS DEFINITION
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
                    method="manual",
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
            self.assertEqual(report["method_effective"], "manual")
            self.assertEqual(report["size_t"], 2)
            self.assertEqual(report["size_z"], 3)

    def test_fixed_alpha_with_non_manual_method_still_marks_user_provided_alpha(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.tif"

            def fake_load_stack_with_omio(_input_path):
                return self.stack.copy(), self.metadata.copy()

            def fake_write_stack_with_omio(output_path, stack, metadata):
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
                    alpha=0.2,
                    alpha_mode="fixed",
                    method="mean_ratio",
                    verbose=False,
                )

            report = json.loads(
                report_path_from_output_path(output_path).read_text(encoding="utf-8")
            )
            self.assertEqual(report["method"], "mean_ratio")
            self.assertEqual(report["method_effective"], "manual")
            self.assertEqual(report["alpha_source"], "user_provided")

    def test_default_alpha_mode_uses_fixed_when_alpha_is_provided(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "default_fixed_output.tif"

            def fake_load_stack_with_omio(_input_path):
                return self.stack.copy(), self.metadata.copy()

            def fake_write_stack_with_omio(output_path, stack, metadata):
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
                    alpha=0.2,
                    method="mean_ratio",
                    verbose=False,
                )

            report = json.loads(
                report_path_from_output_path(output_path).read_text(encoding="utf-8")
            )
            self.assertEqual(report["alpha_mode"], "fixed")
            self.assertIsNone(report["alpha_mode_requested"])
            self.assertTrue(report["alpha_mode_was_defaulted"])
            self.assertEqual(report["method"], "mean_ratio")
            self.assertEqual(report["method_effective"], "manual")

    def test_default_alpha_mode_uses_reference_t_for_non_manual_without_alpha(self) -> None:
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
            output_path = Path(tmpdir) / "default_reference_output.tif"

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
                    method="mean_ratio",
                    signal_percentile=50.0,
                    background_percentile=0.0,
                    min_mask_voxels=8,
                    verbose=False,
                )

            np.testing.assert_allclose(written["stack"][0, 0, 1, :, :], 2.0, atol=1e-6)
            report = json.loads(
                report_path_from_output_path(output_path).read_text(encoding="utf-8")
            )
            self.assertEqual(report["alpha_mode"], "reference_t")
            self.assertIsNone(report["alpha_mode_requested"])
            self.assertTrue(report["alpha_mode_was_defaulted"])
            self.assertAlmostEqual(report["alpha"], 0.25, places=6)

    def test_fixed_bidirectional_unmixing_recovers_both_channels(self) -> None:
        true_source = np.array(
            [
                [4.0, 5.0, 6.0, 7.0, 8.0],
                [5.0, 6.0, 7.0, 8.0, 9.0],
                [6.0, 7.0, 8.0, 9.0, 10.0],
                [7.0, 8.0, 9.0, 10.0, 11.0],
            ],
            dtype=np.float32,
        )
        true_target = np.array(
            [
                [2.0, 3.0, 2.0, 3.0, 2.0],
                [3.0, 2.0, 3.0, 2.0, 3.0],
                [2.0, 3.0, 2.0, 3.0, 2.0],
                [3.0, 2.0, 3.0, 2.0, 3.0],
            ],
            dtype=np.float32,
        )
        alpha_forward = 0.2
        alpha_reverse = 0.1

        stack = np.zeros((1, 2, 2, 4, 5), dtype=np.float32)
        stack[:, :, 0, :, :] = true_source
        stack[:, :, 1, :, :] = true_target
        measured_source = stack[:, :, 0, :, :] + alpha_reverse * stack[:, :, 1, :, :]
        measured_target = stack[:, :, 1, :, :] + alpha_forward * stack[:, :, 0, :, :]
        stack[:, :, 0, :, :] = measured_source
        stack[:, :, 1, :, :] = measured_target
        metadata = {"axes": "TZCYX", "shape": stack.shape}
        written = {}

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "bidirectional_output.tif"

            def fake_load_stack_with_omio(_input_path):
                return stack.copy(), metadata.copy()

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
                    alpha=alpha_forward,
                    alpha_reverse=alpha_reverse,
                    alpha_mode="fixed",
                    bidirectional=True,
                    method="manual",
                    method_reverse="manual",
                    verbose=False,
                )

            expected_source = np.broadcast_to(true_source, written["stack"][:, :, 0, :, :].shape)
            expected_target = np.broadcast_to(true_target, written["stack"][:, :, 1, :, :].shape)
            np.testing.assert_allclose(written["stack"][:, :, 0, :, :], expected_source, atol=1e-6)
            np.testing.assert_allclose(written["stack"][:, :, 1, :, :], expected_target, atol=1e-6)
            report = json.loads(
                report_path_from_output_path(output_path).read_text(encoding="utf-8")
            )
            self.assertTrue(report["bidirectional"])
            self.assertAlmostEqual(report["alpha"], alpha_forward)
            self.assertAlmostEqual(report["alpha_reverse"], alpha_reverse)
            self.assertEqual(report["method_reverse_effective"], "manual")

    def test_fixed_bidirectional_inherits_reverse_alpha_when_not_provided(self) -> None:
        true_source = np.full((1, 1, 4, 5), 6.0, dtype=np.float32)
        true_target = np.full((1, 1, 4, 5), 2.0, dtype=np.float32)
        alpha_value = 0.2

        measured_source = true_source + alpha_value * true_target
        measured_target = true_target + alpha_value * true_source
        stack = np.zeros((1, 1, 2, 4, 5), dtype=np.float32)
        stack[:, :, 0, :, :] = measured_source
        stack[:, :, 1, :, :] = measured_target
        metadata = {"axes": "TZCYX", "shape": stack.shape}
        written = {}

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "bidirectional_inherited_output.tif"

            def fake_load_stack_with_omio(_input_path):
                return stack.copy(), metadata.copy()

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
                    alpha=alpha_value,
                    alpha_mode="fixed",
                    bidirectional=True,
                    method="manual",
                    verbose=False,
                )

            np.testing.assert_allclose(written["stack"][:, :, 0, :, :], true_source, atol=1e-6)
            np.testing.assert_allclose(written["stack"][:, :, 1, :, :], true_target, atol=1e-6)
            report = json.loads(
                report_path_from_output_path(output_path).read_text(encoding="utf-8")
            )
            self.assertTrue(report["alpha_reverse_inherited_from_forward"])
            self.assertAlmostEqual(report["alpha_reverse"], alpha_value)

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
            self.assertEqual(report["method"], "mean_ratio")

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

    def test_manual_method_requires_fixed_alpha_mode(self) -> None:
        with self.assertRaises(ValueError):
            unmix(
                input_path="input.tif",
                output_path="other.tif",
                alpha_mode="reference_t",
                method="manual",
            )

    def test_manual_method_requires_alpha_when_alpha_mode_is_none(self) -> None:
        with self.assertRaises(ValueError):
            unmix(
                input_path="input.tif",
                output_path="other.tif",
                alpha_mode=None,
                method="manual",
            )

    def test_explicit_fixed_alpha_mode_still_requires_alpha(self) -> None:
        def fake_load_stack_with_omio(_input_path):
            return self.stack.copy(), self.metadata.copy()

        def fake_write_stack_with_omio(output_path, stack, metadata):
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
            self.assertRaises(ValueError),
        ):
            unmix(
                input_path="input.tif",
                output_path="other.tif",
                alpha_mode="fixed",
                method="mean_ratio",
                verbose=False,
            )

    def test_unmix_picasso_preserves_shape_and_reduces_channel_dependence(self) -> None:
        rng = np.random.default_rng(7)
        fluorophores = np.zeros((3, 1, 24, 24), dtype=np.float32)
        fluorophores[0, 0, 3:9, 3:9] = 12.0
        fluorophores[1, 0, 12:18, 5:12] = 15.0
        fluorophores[2, 0, 8:16, 14:21] = 10.0
        fluorophores += rng.normal(0.0, 0.1, size=fluorophores.shape).astype(np.float32)
        fluorophores = np.clip(fluorophores, 0.0, None)

        mixing = np.array(
            [
                [1.0, 0.22, 0.10],
                [0.15, 1.0, 0.18],
                [0.08, 0.20, 1.0],
            ],
            dtype=np.float32,
        )
        measured = np.einsum("ij,jzyx->izyx", mixing, fluorophores, optimize=True)
        stack = measured[np.newaxis, :, :, :, :]
        stack = np.moveaxis(stack, 1, 2)
        metadata = {"axes": "TZCYX", "shape": stack.shape}
        written = {}

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "picasso_output.tif"

            def fake_load_stack_with_omio(_input_path):
                return stack.copy(), metadata.copy()

            def fake_write_stack_with_omio(output_path, saved_stack, metadata):
                written["stack"] = saved_stack.copy()
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
                unmix_picasso(
                    input_path="input.tif",
                    output_path=output_path,
                    channels=[0, 1, 2],
                    alpha_mode="reference_t",
                    alpha_reference_t=0,
                    background_percentile=0.0,
                    mi_bins=16,
                    max_iter=3,
                    max_alpha_voxels=None,
                    verbose=False,
                )

            self.assertEqual(written["stack"].shape, stack.shape)
            self.assertTrue(np.all(np.isfinite(written["stack"])))

            measured_flat = stack[0, 0].reshape(3, -1)
            corrected_flat = written["stack"][0, 0].reshape(3, -1)
            measured_dependence = (
                abs(np.corrcoef(measured_flat[0], measured_flat[1])[0, 1])
                + abs(np.corrcoef(measured_flat[0], measured_flat[2])[0, 1])
                + abs(np.corrcoef(measured_flat[1], measured_flat[2])[0, 1])
            )
            corrected_dependence = (
                abs(np.corrcoef(corrected_flat[0], corrected_flat[1])[0, 1])
                + abs(np.corrcoef(corrected_flat[0], corrected_flat[2])[0, 1])
                + abs(np.corrcoef(corrected_flat[1], corrected_flat[2])[0, 1])
            )
            self.assertLess(corrected_dependence, measured_dependence)

            report = json.loads(
                report_path_from_output_path(output_path).read_text(encoding="utf-8")
            )
            self.assertEqual(report["method"], "picasso")
            self.assertEqual(report["implementation"], "matlab_3c")
            self.assertEqual(report["channels"], [0, 1, 2])
            self.assertEqual(len(report["unmixing_matrix"]), 3)
            self.assertTrue(all(np.isfinite(np.asarray(report["unmixing_matrix"])).ravel()))
            self.assertEqual(report["iterations_run"], 3)

    def test_unmix_picasso_matlab_n_supports_more_than_three_channels(self) -> None:
        rng = np.random.default_rng(11)
        fluorophores = np.zeros((4, 1, 20, 20), dtype=np.float32)
        fluorophores[0, 0, 2:7, 2:7] = 9.0
        fluorophores[1, 0, 10:16, 3:8] = 11.0
        fluorophores[2, 0, 4:10, 11:17] = 8.0
        fluorophores[3, 0, 12:18, 12:18] = 10.0
        fluorophores += rng.normal(0.0, 0.05, size=fluorophores.shape).astype(np.float32)
        fluorophores = np.clip(fluorophores, 0.0, None)

        mixing = np.array(
            [
                [1.0, 0.18, 0.06, 0.05],
                [0.07, 1.0, 0.14, 0.05],
                [0.05, 0.10, 1.0, 0.16],
                [0.03, 0.05, 0.12, 1.0],
            ],
            dtype=np.float32,
        )
        measured = np.einsum("ij,jzyx->izyx", mixing, fluorophores, optimize=True)
        stack = np.moveaxis(measured[np.newaxis, :, :, :, :], 1, 2)
        metadata = {"axes": "TZCYX", "shape": stack.shape}
        written = {}

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "picasso_matlab_n_output.tif"

            def fake_load_stack_with_omio(_input_path):
                return stack.copy(), metadata.copy()

            def fake_write_stack_with_omio(output_path, saved_stack, metadata):
                written["stack"] = saved_stack.copy()
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
                unmix_picasso(
                    input_path="input.tif",
                    output_path=output_path,
                    channels=[0, 1, 2, 3],
                    implementation="matlab_n",
                    alpha_mode="reference_t",
                    alpha_reference_t=0,
                    background_percentile=0.0,
                    max_iter=2,
                    verbose=False,
                )

            self.assertEqual(written["stack"].shape, stack.shape)
            self.assertTrue(np.all(np.isfinite(written["stack"])))
            report = json.loads(
                report_path_from_output_path(output_path).read_text(encoding="utf-8")
            )
            self.assertEqual(report["implementation"], "matlab_n")
            self.assertEqual(report["channels"], [0, 1, 2, 3])
            self.assertEqual(len(report["unmixing_matrix"]), 4)

    def test_unmix_picasso_source_sink_n_writes_finite_output(self) -> None:
        rng = np.random.default_rng(13)
        fluorophores = np.zeros((4, 1, 18, 18), dtype=np.float32)
        fluorophores[0, 0, 2:8, 2:8] = 7.0
        fluorophores[1, 0, 9:15, 2:8] = 9.0
        fluorophores[2, 0, 3:9, 10:16] = 8.0
        fluorophores[3, 0, 10:16, 10:16] = 10.0
        fluorophores += rng.normal(0.0, 0.05, size=fluorophores.shape).astype(np.float32)
        fluorophores = np.clip(fluorophores, 0.0, None)

        mixing = np.array(
            [
                [1.0, 0.15, 0.00, 0.00],
                [0.00, 1.0, 0.00, 0.00],
                [0.00, 0.00, 1.0, 0.10],
                [0.00, 0.00, 0.00, 1.0],
            ],
            dtype=np.float32,
        )
        measured = np.einsum("ij,jzyx->izyx", mixing, fluorophores, optimize=True)
        stack = np.moveaxis(measured[np.newaxis, :, :, :, :], 1, 2)
        metadata = {"axes": "TZCYX", "shape": stack.shape}
        written = {}

        source_sink_matrix = [
            [1, -1, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, -1],
            [0, 0, 0, 1],
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "picasso_source_sink_output.tif"

            def fake_load_stack_with_omio(_input_path):
                return stack.copy(), metadata.copy()

            def fake_write_stack_with_omio(output_path, saved_stack, metadata):
                written["stack"] = saved_stack.copy()
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
                unmix_picasso(
                    input_path="input.tif",
                    output_path=output_path,
                    channels=[0, 1, 2, 3],
                    implementation="source_sink_n",
                    source_sink_matrix=source_sink_matrix,
                    alpha_mode="reference_t",
                    alpha_reference_t=0,
                    background_percentile=0.0,
                    mi_bins=16,
                    alpha_max=1.0,
                    max_alpha_voxels=None,
                    verbose=False,
                )

            self.assertEqual(written["stack"].shape, stack.shape)
            self.assertTrue(np.all(np.isfinite(written["stack"])))
            report = json.loads(
                report_path_from_output_path(output_path).read_text(encoding="utf-8")
            )
            self.assertEqual(report["implementation"], "source_sink_n")
            self.assertEqual(report["channels"], [0, 1, 2, 3])
            self.assertEqual(report["source_sink_matrix"], source_sink_matrix)
            self.assertEqual(len(report["picasso_estimation"]["alpha_parameters"]), 4)

    def test_unmix_picasso_source_sink_n_can_build_matrix_from_roles(self) -> None:
        rng = np.random.default_rng(17)
        fluorophores = np.zeros((3, 1, 16, 16), dtype=np.float32)
        fluorophores[0, 0, 2:7, 2:7] = 8.0
        fluorophores[1, 0, 8:13, 8:13] = 7.0
        fluorophores[2, 0, 4:10, 10:15] = 6.0
        fluorophores += rng.normal(0.0, 0.05, size=fluorophores.shape).astype(np.float32)
        fluorophores = np.clip(fluorophores, 0.0, None)

        mixing = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.18, 1.0, 0.09],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float32,
        )
        measured = np.einsum("ij,jzyx->izyx", mixing, fluorophores, optimize=True)
        stack = np.moveaxis(measured[np.newaxis, :, :, :, :], 1, 2)
        metadata = {"axes": "TZCYX", "shape": stack.shape}

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "picasso_source_sink_roles_output.tif"

            def fake_load_stack_with_omio(_input_path):
                return stack.copy(), metadata.copy()

            def fake_write_stack_with_omio(output_path, saved_stack, metadata):
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
                unmix_picasso(
                    input_path="input.tif",
                    output_path=output_path,
                    channels=[0, 1, 2],
                    implementation="source_sink_n",
                    sink_channels=[1],
                    neutral_channels=[2],
                    alpha_mode="reference_t",
                    alpha_reference_t=0,
                    background_percentile=0.0,
                    mi_bins=16,
                    alpha_max=1.0,
                    verbose=False,
                )

            report = json.loads(
                report_path_from_output_path(output_path).read_text(encoding="utf-8")
            )
            self.assertEqual(report["sink_channels"], [1])
            self.assertEqual(report["neutral_channels"], [2])
            self.assertEqual(
                report["source_sink_matrix"],
                [
                    [1, -1, 0],
                    [0, 1, 0],
                    [0, 0, 1],
                ],
            )

# %% MAIN
if __name__ == "__main__":
    unittest.main()
# %% END
