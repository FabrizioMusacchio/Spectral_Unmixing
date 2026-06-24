"""
Unit tests for napari viewer helpers.

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

import spectral_unmixing.viewer as viewer_module


class FakeLayer:
    def __init__(self, data, name, scale, colormap, blending, opacity, contrast_limits):
        self.data = data
        self.name = name
        self.scale = scale
        self.colormap = colormap
        self.blending = blending
        self.opacity = opacity
        self.contrast_limits = contrast_limits
        self.visible = True


class FakeViewer:
    def __init__(self, title=""):
        self.title = title
        self.layers = []

    def add_image(
        self,
        data,
        name,
        scale,
        colormap,
        blending,
        opacity,
        contrast_limits,
    ):
        layer = FakeLayer(data, name, scale, colormap, blending, opacity, contrast_limits)
        self.layers.append(layer)
        return layer


class FakeNapari:
    def __init__(self):
        self._current_viewer = None

    def current_viewer(self):
        return self._current_viewer

    def Viewer(self, title=""):
        self._current_viewer = FakeViewer(title=title)
        return self._current_viewer


class ViewerTests(unittest.TestCase):
    def setUp(self) -> None:
        viewer_module._VIEWER = None

    def test_show_unmixed_channels_reuses_viewer_and_updates_layers(self) -> None:
        fake_napari = FakeNapari()
        stack = np.zeros((2, 3, 2, 4, 5), dtype=np.float32)
        stack[:, :, 0, :, :] = 5.0
        stack[:, :, 1, :, :] = 2.0
        metadata = {
            "axes": "TZCYX",
            "TimeIncrement": 1.0,
            "PhysicalSizeZ": 2.0,
            "PhysicalSizeY": 3.0,
            "PhysicalSizeX": 4.0,
        }

        updated_stack = stack.copy()
        updated_stack[:, :, 1, :, :] = 7.0

        with patch("spectral_unmixing.viewer.import_napari", return_value=fake_napari):
            with patch(
                "spectral_unmixing.viewer.load_stack_with_omio",
                side_effect=[(stack, metadata), (updated_stack, metadata)],
            ):
                viewer = viewer_module.show_unmixed_channels_in_napari(
                    "dummy_output.tif",
                    layer_prefix="Fixed alpha",
                )
                viewer_again = viewer_module.show_unmixed_channels_in_napari(
                    "dummy_output.tif",
                    layer_prefix="Fixed alpha",
                )

        self.assertIs(viewer, viewer_again)
        self.assertEqual(len(viewer.layers), 2)
        self.assertEqual(viewer.layers[0].name, "Fixed alpha | source C0")
        self.assertEqual(viewer.layers[1].name, "Fixed alpha | target C1")
        np.testing.assert_allclose(viewer.layers[1].data, updated_stack[:, :, 1, :, :])
        self.assertEqual(viewer.layers[0].scale, (1.0, 2.0, 3.0, 4.0))


if __name__ == "__main__":
    unittest.main()
