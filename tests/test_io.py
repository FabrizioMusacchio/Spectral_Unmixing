"""
Unit tests for OMIO writer path handling.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from spectral_unmixing.io import write_stack_with_omio

# %% CLASS DEFINITION
class WriteStackWithOmioTests(unittest.TestCase):
    def test_moves_omio_reported_ome_tif_to_requested_output_path(self) -> None:
        stack = np.zeros((1, 1, 2, 4, 5), dtype=np.float32)
        metadata = {
            "axes": "TZCYX",
            "shape": stack.shape,
            "PhysicalSizeX": 1.0,
            "PhysicalSizeY": 1.0,
            "PhysicalSizeZ": 1.0,
            "unit": "micron",
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            output_path = tmpdir_path / "requested_output.tif"
            reported_path = tmpdir_path / "source_name.ome.tif"

            class FakeOmio:
                @staticmethod
                def imwrite(*args, **kwargs):
                    reported_path.write_bytes(b"fake-tif")
                    return [str(reported_path)]

            with patch("spectral_unmixing.io.import_omio", return_value=FakeOmio()):
                final_path = write_stack_with_omio(output_path, stack, metadata)

            self.assertEqual(final_path, output_path)
            self.assertTrue(output_path.exists())
            self.assertFalse(reported_path.exists())
            self.assertEqual(output_path.read_bytes(), b"fake-tif")

    def test_replaces_existing_requested_output_with_fresh_omio_result(self) -> None:
        stack = np.zeros((1, 1, 2, 4, 5), dtype=np.float32)
        metadata = {
            "axes": "TZCYX",
            "shape": stack.shape,
            "PhysicalSizeX": 1.0,
            "PhysicalSizeY": 1.0,
            "PhysicalSizeZ": 1.0,
            "unit": "micron",
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            output_path = tmpdir_path / "requested_output.tif"
            reported_path = tmpdir_path / "fresh_write.ome.tif"
            output_path.write_bytes(b"stale-tif")

            class FakeOmio:
                @staticmethod
                def imwrite(*args, **kwargs):
                    reported_path.write_bytes(b"fresh-tif")
                    return [str(reported_path)]

            with patch("spectral_unmixing.io.import_omio", return_value=FakeOmio()):
                final_path = write_stack_with_omio(output_path, stack, metadata)

            self.assertEqual(final_path, output_path)
            self.assertTrue(output_path.exists())
            self.assertFalse(reported_path.exists())
            self.assertEqual(output_path.read_bytes(), b"fresh-tif")

# %% MAIN
if __name__ == "__main__":
    unittest.main()
# %% END
