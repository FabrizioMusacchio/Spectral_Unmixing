"""
Generate figures for the spectral-unmixing preprint.

This script uses externally sourced example microscopy datasets included with
the repository to create figure panels for the manuscript. If the expected
processed example outputs are missing, it recreates them with the same
settings used in the tutorials.

Outputs are written to ``papers/preprint/figures`` and
``papers/preprint/results``.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

import csv
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

# PATH SETUP:
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# The benchmark figures are generated headlessly in the terminal
# On normal user machines this workaround is typically unnecessary, but it is
# harmless and keeps the script robust when napari is imported indirectly by OMIO.
os.environ.setdefault("XDG_CACHE_HOME", str(PROJECT_ROOT / ".cache" / "xdg"))
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / ".cache" / "matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linear_sum_assignment

from spectral_unmixing.estimation import mutual_information_1d
from spectral_unmixing.io import convert_time_encoded_stack_to_channel_stack, load_stack_with_omio
from spectral_unmixing.unmixing import unmix, unmix_picasso

# set global font to Arial:
plt.rcParams["font.family"] = "Arial"
# %% PATHS
PAPER_DIR = PROJECT_ROOT / "papers" / "preprint"
FIGURE_DIR = PAPER_DIR / "figures"
RESULTS_DIR = PAPER_DIR / "results"

EXAMPLE_DIR = PROJECT_ROOT / "example_data" / "PICASSO_examples"
UNMIXED_DIR = EXAMPLE_DIR / "unmixed"

TWO_COLOR_INPUT = EXAMPLE_DIR / "2_color_unmixing_validation.tif"
TWO_COLOR_OUTPUT = UNMIXED_DIR / "2_color_unmixing_validation_unmixed_fixed_alpha.tif"

FIVE_COLOR_INPUT = EXAMPLE_DIR / "5_color_unmixing_simulation.tif"
FIVE_COLOR_GROUND_TRUTH_INPUT = EXAMPLE_DIR / "5_color_unmixing_simulation_ground_truth.tif"
FIVE_COLOR_GROUND_TRUTH_CONVERTED = EXAMPLE_DIR / "5_color_unmixing_simulation_ground_truth_converted.tif"
FIVE_COLOR_MATLAB_N_OUTPUT = UNMIXED_DIR / "5_color_unmixing_simulation_picasso_matlab_n.tif"

GFAP_SOURCE_SINK_INPUT = EXAMPLE_DIR / "GFAP_sink_LMNB1_source.tif"
GFAP_MATLAB_N_OUTPUT = UNMIXED_DIR / "GFAP_sink_LMNB1_source_picasso_matlab_n.tif"
GFAP_SOURCE_SINK_OUTPUT = UNMIXED_DIR / "GFAP_sink_LMNB1_source_picasso_source_sink.tif"

DISPLAY_COLORS: list[tuple[str, tuple[float, float, float]]] = [
    ("cyan", (0.0, 1.0, 1.0)),
    ("magenta", (1.0, 0.0, 1.0)),
    ("yellow", (1.0, 0.9, 0.1)),
    ("orange", (1.0, 0.55, 0.1)),
    ("green", (0.2, 1.0, 0.2)),
    ("red", (1.0, 0.2, 0.2))]

PICASSO_PLOT_COLORS: dict[str, str] = {
    "mixed": "#BAB0AC",
    "matlab_3c": "#4E79A7",
    "matlab_n": "#59A14F",
    "source_sink_n": "#E15759",
    "ground_truth": "#B07AA1"}

DIRECTED_METHOD_COLORS: dict[str, str] = {
    "manual": "#4E79A7",
    "mean_ratio": "#76B7B2",
    "linear_fit": "#59A14F",
    "corr_min": "#F28E2B",
    "mi_min": "#B07AA1"}
# %% PANEL PLOT PARAMETERS
CM2IN = 1.0 / 2.54
PANEL_FONT_SIZE = 10

EXAMPLE_PANEL_PARAMS: dict[str, dict[str, dict[str, object]]] = {
    "figure_public_2color_example": {
        "image_panels": {
            "figsize_cm": (6.5, 6.5),
            "show_scalebar": False,
            "show_scalebar_unit": False,
            "scalebar_position": "lower left"},
        "correlation_sum": {
            "savename": "panel_correlation_sum.png",
            "figsize_cm": (3.05, 5.0),
            "title": "correlation sum",
            "ylabel": "lower is better",
            "legend_show": False,
            "legend_loc": "best",
            "spines": {"top": False, "right": False, "left": False, "bottom": True},
            "xrotation": 20,
            "ylim": None,
            "yticks": None,
            "ytick_length": 3},
        "mutual_information": {
            "savename": "panel_mutual_information.png",
            "figsize_cm": (3.05, 5.0),
            "title": "mutual information",
            "ylabel": "lower is better",
            "legend_show": False,
            "legend_loc": "best",
            "spines": {"top": False, "right": False, "left": False, "bottom": True},
            "xrotation": 20,
            "ylim": None,
            "yticks": [0,0.05, 0.10, 0.15],
            "ytick_length": 3},
        "correlation_sum_all_methods": {
            "savename": "panel_correlation_sum_all_methods.png",
            "figsize_cm": (5.8, 5.1),
            "title": "correlation sum",
            "ylabel": "pairwise sum |r|",
            "legend_show": False,
            "legend_loc": "best",
            "spines": {"top": False, "right": False, "left": False, "bottom": True},
            "xrotation": 30,
            "ylim": None,
            "yticks": None,
            "ytick_length": 3},
        "mutual_information_all_methods": {
            "savename": "panel_mutual_information_all_methods.png",
            "figsize_cm": (5.8, 5.1),
            "title": "mutual information",
            "ylabel": "pairwise MI",
            "legend_show": False,
            "legend_loc": "best",
            "spines": {"top": False, "right": False, "left": False, "bottom": True},
            "xrotation": 30,
            "ylim": None,
            "yticks": None,
            "ytick_length": 3},
        "negative_fraction_all_methods": {
            "savename": "panel_negative_fraction_all_methods.png",
            "figsize_cm": (5.1, 5.1),
            "title": "neg. fraction before clipping     ",
            "ylabel": r"$f_{\mathrm{neg}}$",
            "legend_show": False,
            "legend_loc": "best",
            "spines": {"top": False, "right": False, "left": False, "bottom": True},
            "xrotation": 30,
            "ylim": (0.0, 1.1),
            "yticks": [0, 0.5, 1],
            "ytick_length": 3},
    },
    "figure_public_5color_example": {
        "image_panels": {
            "figsize_cm": (6.5, 6.5),
            "show_scalebar": False,
            "show_scalebar_unit": False,
            "scalebar_position": "lower left"},
        "correlation_sum": {
            "savename": "panel_correlation_sum.png",
            "figsize_cm": (3.05, 5.10),
            "title": "correlation",
            "ylabel": "pairwise sum |r|",
            "legend_show": False,
            "legend_loc": "best",
            "spines": {"top": False, "right": False, "left": False, "bottom": True},
            "xrotation": 35,
            "ylim": None,
            "yticks": None,
            "ytick_length": 3},
        "mutual_information_sum": {
            "savename": "panel_mutual_information_sum.png",
            "figsize_cm": (3.35, 5.10),
            "title": "mutual inform.   ",
            "ylabel": "pairwise MI",
            "legend_show": False,
            "legend_loc": "best",
            "spines": {"top": False, "right": False, "left": False, "bottom": True},
            "xrotation": 35,
            "ylim": None,
            "yticks": None,
            "ytick_length": 3},
        "gt_matched_recovery": {
            "savename": "panel_gt_matched_recovery.png",
            "figsize_cm": (2.65, 5.10),
            "title": "recovery",
            "ylabel": "best-match |r|",
            "legend_show": False,
            "legend_loc": "best",
            "spines": {"top": False, "right": False, "left": False, "bottom": True},
            "xrotation": 35,
            "ylim": (0.0, 1.05),
            "yticks": [0, 0.5, 1.0],
            "ytick_length": 3},
        "direct_gt_channel_correlation": {
            "savename": "panel_direct_channel_correlation.png",
            "figsize_cm": (4.2, 4.60),
            "title": "direct channel correlation       ",
            "ylabel": "Pearson r",
            "legend_show": False,
            "legend_loc": "outer top",
            "spines": {"top": False, "right": False, "left": False, "bottom": True},
            "xrotation": 0,
            "ylim": None,
            "yticks": [0,1],
            "ytick_length": 3},
    },
    "figure_public_source_sink_example": {
        "image_panels": {
            "figsize_cm": (6.5, 6.5),
            "show_scalebar": False,
            "show_scalebar_unit": False,
            "scalebar_position": "lower left"},
        "correlation_sum": {
            "savename": "panel_correlation_sum.png",
            "figsize_cm": (3.35, 5.10),
            "title": "correlation",
            "ylabel": "pairwise sum |r|",
            "legend_show": False,
            "legend_loc": "best",
            "spines": {"top": False, "right": False, "left": False, "bottom": True},
            "xrotation": 35,
            "ylim": None,
            "yticks": None,
            "ytick_length": 3},
        "mutual_information": {
            "savename": "panel_mutual_information.png",
            "figsize_cm": (3.35, 5.10),
            "title": "mutual inform.    ",
            "ylabel": "pairwise MI",
            "legend_show": False,
            "legend_loc": "best",
            "spines": {"top": False, "right": False, "left": False, "bottom": True},
            "xrotation": 35,
            "ylim": None,
            "yticks": None,
            "ytick_length": 3},
    },
}
# %% SMALL HELPERS
def _ensure_output_dirs() -> None:
    """Create the preprint output directories if they do not yet exist."""

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

def _report_has_expected_value(
    output_path: Path,
    key: str,
    expected_value: object,
) -> bool:
    """Return ``True`` when the output sidecar exists and stores the expected value."""

    report_path = output_path.with_suffix(output_path.suffix + ".json")
    if not output_path.exists() or not report_path.exists():
        return False

    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False

    return report.get(key) == expected_value

def _load_sidecar_report(output_path: Path) -> dict[str, object]:
    """Load and return the JSON sidecar report for one generated output stack."""

    report_path = output_path.with_suffix(output_path.suffix + ".json")
    return json.loads(report_path.read_text(encoding="utf-8"))

def _normalize_with_bounds(
    image: np.ndarray,
    *,
    low_value: float,
    high_value: float,
    gamma: float = 0.8,
) -> np.ndarray:
    """Normalize one image to ``[0, 1]`` with explicit display bounds."""

    image = np.asarray(image, dtype=np.float32)
    if high_value <= low_value:
        return np.zeros_like(image, dtype=np.float32)
    normalized = np.clip((image - low_value) / (high_value - low_value), 0.0, 1.0)
    if gamma != 1.0:
        normalized = normalized ** float(gamma)
    return normalized

def _normalize_for_display(
    image: np.ndarray,
    *,
    low: float = 0.5,
    high: float = 99.8,
    gamma: float = 0.8,
) -> np.ndarray:
    """Normalize one image to ``[0, 1]`` for qualitative display."""

    image = np.asarray(image, dtype=np.float32)
    low_value = float(np.percentile(image, low))
    high_value = float(np.percentile(image, high))
    return _normalize_with_bounds(image, low_value=low_value, high_value=high_value, gamma=gamma)

def _compute_display_bounds(
    image: np.ndarray,
    *,
    low: float = 0.5,
    high: float = 99.8,
) -> tuple[float, float]:
    """Return percentile-based lower and upper display bounds."""

    image = np.asarray(image, dtype=np.float32)
    return float(np.percentile(image, low)), float(np.percentile(image, high))

def _colorize_channel(
    image: np.ndarray,
    color: tuple[float, float, float],
    *,
    low_value: float | None = None,
    high_value: float | None = None,
    low: float = 0.5,
    high: float = 99.8,
    gamma: float = 0.8,
) -> np.ndarray:
    """Render one scalar channel as an RGB image with the requested tint."""

    if low_value is None or high_value is None:
        low_value, high_value = _compute_display_bounds(image, low=low, high=high)
    normalized = _normalize_with_bounds(
        image,
        low_value=float(low_value),
        high_value=float(high_value),
        gamma=gamma,
    )
    rgb = np.zeros(normalized.shape + (3,), dtype=np.float32)
    rgb[..., 0] = normalized * float(color[0])
    rgb[..., 1] = normalized * float(color[1])
    rgb[..., 2] = normalized * float(color[2])
    return np.clip(rgb, 0.0, 1.0)

def _cyan_magenta_composite(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Create a simple cyan-magenta RGB composite from two 2D channels."""

    source_n = _normalize_for_display(source)
    target_n = _normalize_for_display(target)

    rgb = np.zeros(source_n.shape + (3,), dtype=np.float32)
    rgb[..., 0] = target_n
    rgb[..., 1] = source_n
    rgb[..., 2] = np.maximum(source_n, target_n)
    return np.clip(rgb, 0.0, 1.0)

def _multichannel_composite(
    channels: np.ndarray,
    *,
    colors: list[tuple[float, float, float]],
    bounds: list[tuple[float, float]] | None = None,
    gamma: float = 0.8,
) -> np.ndarray:
    """Create one RGB composite from a stack of 2D channels."""

    composite = np.zeros(np.asarray(channels[0]).shape + (3,), dtype=np.float32)
    for channel_index, image in enumerate(np.asarray(channels, dtype=np.float32)):
        color = colors[channel_index % len(colors)]
        low_value, high_value = (bounds[channel_index] if bounds is not None else (None, None))
        composite += _colorize_channel(
            image,
            color,
            low_value=low_value,
            high_value=high_value,
            gamma=gamma,
        )
    return np.clip(composite, 0.0, 1.0)

def _pairwise_abs_correlation_sum(channels_first: np.ndarray) -> float:
    """Return the sum of absolute Pearson correlations across channel pairs."""

    flattened = np.asarray(channels_first, dtype=np.float64).reshape(channels_first.shape[0], -1)
    total = 0.0
    for i in range(flattened.shape[0]):
        for j in range(i + 1, flattened.shape[0]):
            total += abs(float(np.corrcoef(flattened[i], flattened[j])[0, 1]))
    return float(total)

def _pairwise_mutual_information_sum(channels_first: np.ndarray, *, bins: int = 32) -> float:
    """Return the sum of pairwise mutual-information estimates across channels."""

    flattened = np.asarray(channels_first, dtype=np.float64).reshape(channels_first.shape[0], -1)
    total = 0.0
    for i in range(flattened.shape[0]):
        for j in range(i + 1, flattened.shape[0]):
            total += float(mutual_information_1d(flattened[i], flattened[j], bins=bins))
    return float(total)

def _best_match_channel_correlation(recovered: np.ndarray, truth: np.ndarray) -> float:
    """Return mean absolute Pearson correlation after Hungarian channel matching."""

    recovered_flat = np.asarray(recovered, dtype=np.float64).reshape(recovered.shape[0], -1)
    truth_flat = np.asarray(truth, dtype=np.float64).reshape(truth.shape[0], -1)

    correlation_matrix = np.zeros((recovered_flat.shape[0], truth_flat.shape[0]), dtype=np.float64)
    for i in range(recovered_flat.shape[0]):
        for j in range(truth_flat.shape[0]):
            correlation_matrix[i, j] = float(np.corrcoef(recovered_flat[i], truth_flat[j])[0, 1])

    row_indices, col_indices = linear_sum_assignment(1.0 - np.abs(correlation_matrix))
    return float(np.mean(np.abs(correlation_matrix[row_indices, col_indices])))

def _direct_channel_correlation_vector(recovered: np.ndarray, truth: np.ndarray) -> np.ndarray:
    """Return one same-index Pearson-correlation value per channel."""

    recovered_flat = np.asarray(recovered, dtype=np.float64).reshape(recovered.shape[0], -1)
    truth_flat = np.asarray(truth, dtype=np.float64).reshape(truth.shape[0], -1)
    if recovered_flat.shape != truth_flat.shape:
        raise ValueError(
            "recovered and truth must have the same channel-first shape. "
            f"Got {recovered.shape!r} and {truth.shape!r}."
        )

    correlations = np.zeros((recovered_flat.shape[0],), dtype=np.float64)
    for channel_index in range(recovered_flat.shape[0]):
        correlations[channel_index] = float(
            np.corrcoef(recovered_flat[channel_index], truth_flat[channel_index])[0, 1]
        )
    return correlations

def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    """Write one CSV table from a list of dictionaries."""

    fieldnames: list[str] = []
    known = set()
    for row in rows:
        for key in row:
            if key not in known:
                fieldnames.append(key)
                known.add(key)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

def _despine_axis(axis) -> None:
    """Remove box-like spines for cleaner publication plots."""

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(False)
    axis.tick_params(axis="y", length=0)

def _panel_directory(figure_path: Path) -> Path:
    """Return and create the panel-output directory for one figure file."""

    panel_dir = FIGURE_DIR / figure_path.stem
    panel_dir.mkdir(parents=True, exist_ok=True)
    return panel_dir

def _panel_pdf_directory(panel_dir: Path) -> Path:
    """Return and create the PDF output directory for one panel directory."""

    pdf_dir = panel_dir / "PDF"
    pdf_dir.mkdir(parents=True, exist_ok=True)
    return pdf_dir

def _cm_to_in(figsize_cm: tuple[float, float]) -> tuple[float, float]:
    """Convert a figure size from centimeters to inches."""

    return float(figsize_cm[0]) * CM2IN, float(figsize_cm[1]) * CM2IN

def _get_panel_params(
    config_root: dict[str, dict[str, dict[str, object]]],
    figure_key: str,
    panel_key: str,
) -> dict[str, object]:
    """Return the configured parameter block for one figure/panel combination."""

    return dict(config_root.get(figure_key, {}).get(panel_key, {}))

def _apply_panel_axis_style(
    axis,
    config: dict[str, object],
    *,
    default_title: str | None = None,
    default_ylabel: str | None = None,
) -> None:
    """Apply configurable styling to one standalone summary axis."""

    title = config.get("title", default_title)
    ylabel = config.get("ylabel", default_ylabel)
    if title:
        axis.set_title(str(title), fontsize=PANEL_FONT_SIZE)
    if ylabel:
        axis.set_ylabel(str(ylabel))

    xrotation = config.get("xrotation")
    if xrotation is not None:
        axis.tick_params(axis="x", labelrotation=float(xrotation))
        if float(xrotation) != 0.0:
            for tick_label in axis.get_xticklabels():
                tick_label.set_horizontalalignment("right")
                tick_label.set_rotation_mode("anchor")

    ylim = config.get("ylim")
    if ylim is not None:
        axis.set_ylim(*ylim)

    yticks = config.get("yticks")
    if yticks is not None:
        axis.set_yticks(list(yticks))

    ytick_length = config.get("ytick_length")
    if ytick_length is not None:
        axis.tick_params(axis="y", length=float(ytick_length))

    spines = config.get("spines", {})
    for side in ("top", "right", "left", "bottom"):
        axis.spines[side].set_visible(bool(spines.get(side, True)))
    if ytick_length is None and not bool(spines.get("left", True)):
        axis.tick_params(axis="y", length=0)

    legend = axis.get_legend()
    show_legend = bool(config.get("legend_show", False))
    if legend is not None:
        legend.remove()
    if show_legend:
        handles, labels = axis.get_legend_handles_labels()
        if handles:
            loc = str(config.get("legend_loc", "best"))
            if loc == "outer top":
                axis.legend(
                    handles,
                    labels,
                    loc="upper left",
                    bbox_to_anchor=(1.02, 1.0),
                    borderaxespad=0.0,
                    frameon=False,
                )
            else:
                axis.legend(handles, labels, loc=loc, frameon=False)

def _choose_scalebar(
    image_width_px: int,
    *,
    pixel_size_x: float | None = None,
    unit: str | None = None,
    fraction: float = 0.2,
) -> tuple[float, str]:
    """Choose one visually reasonable scalebar length in pixels and as label."""

    target_px = max(1.0, image_width_px * fraction)
    if pixel_size_x is not None and np.isfinite(pixel_size_x) and pixel_size_x > 0:
        target_units = target_px * float(pixel_size_x)
        exponent = np.floor(np.log10(target_units)) if target_units > 0 else 0.0
        base = 10.0 ** exponent
        candidates = np.asarray([1.0, 2.0, 5.0, 10.0]) * base
        length_units = float(candidates[np.argmin(np.abs(candidates - target_units))])
        length_px = length_units / float(pixel_size_x)
        label_unit = unit or "um"
        label_unit = "μm" if label_unit == "um" else label_unit
        return float(length_px), f"{length_units:g} {label_unit}"

    length_px = float(max(10, int(round(target_px / 10.0) * 10)))
    return length_px, f"{int(round(length_px))} px"

def _metadata_pixel_size_x_and_unit(metadata: dict[str, object] | None) -> tuple[float | None, str | None]:
    """Extract one horizontal physical pixel size and unit from OMIO metadata."""

    if metadata is None:
        return None, None
    pixel_size_x = metadata.get("PhysicalSizeX")
    unit = metadata.get("PhysicalSizeXUnit")
    try:
        pixel_size_value = float(pixel_size_x) if pixel_size_x is not None else None
    except (TypeError, ValueError):
        pixel_size_value = None
    return pixel_size_value, str(unit) if unit is not None else None

def _add_scalebar_to_image_axis(
    axis,
    image_shape: tuple[int, int],
    config: dict[str, object],
    *,
    pixel_size_x: float | None = None,
    unit: str | None = None,
) -> None:
    """Add one configurable corner scalebar to a standalone image panel."""

    if not bool(config.get("show_scalebar", False)):
        return

    height_px, width_px = image_shape[:2]
    length_px, label = _choose_scalebar(width_px, pixel_size_x=pixel_size_x, unit=unit)
    margin_x = width_px * 0.06
    margin_y = height_px * 0.07
    position = str(config.get("scalebar_position", "lower right")).lower()

    if position not in {"lower left", "lower right", "upper left", "upper right"}:
        raise ValueError(
            "scalebar_position must be one of 'lower left', 'lower right', "
            f"'upper left', or 'upper right', got {position!r}."
        )

    if "left" in position:
        x0 = margin_x
        x1 = x0 + length_px
    else:
        x1 = width_px - margin_x
        x0 = x1 - length_px

    y = margin_y if "upper" in position else height_px - margin_y
    axis.plot([x0, x1], [y, y], color="white", linewidth=3.0, solid_capstyle="butt")
    if bool(config.get("show_scalebar_unit", True)):
        text_y_offset = height_px * 0.04
        text_y = y + text_y_offset if "upper" in position else y - text_y_offset
        text_va = "top" if "upper" in position else "bottom"
        axis.text(
            (x0 + x1) / 2.0,
            text_y,
            label,
            color="white",
            ha="center",
            va=text_va,
            fontsize=10,
        )


@dataclass(slots=True)
class PanelExportSpec:
    """Describe one standalone panel export for a manuscript figure."""

    filename: str
    render_fn: Callable[[object], None] | None = None
    axis: object | None = None
    tight: bool = True
    figsize: tuple[float, float] = (3.4, 3.0)
    full_bleed: bool = False

def _save_axis_panels(
    figure,
    panel_specs: list[PanelExportSpec | tuple[str, object, bool]],
    *,
    figure_path: Path,
    pad_scale: tuple[float, float] = (1.02, 1.04),
) -> None:
    """Save panel images for the requested axes."""

    panel_dir = _panel_directory(figure_path)
    pdf_dir = _panel_pdf_directory(panel_dir)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()

    for spec in panel_specs:
        if isinstance(spec, PanelExportSpec):
            if spec.render_fn is not None:
                panel_figure, panel_axis = plt.subplots(figsize=_cm_to_in(spec.figsize), constrained_layout=False)
                spec.render_fn(panel_axis)
                if spec.full_bleed:
                    panel_axis.set_position([0.0, 0.0, 1.0, 1.0])
                    panel_figure.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
                else:
                    panel_figure.tight_layout(pad=0.3)
                panel_figure.savefig(
                    panel_dir / spec.filename,
                    dpi=300,
                    bbox_inches=None,
                    pad_inches=0.0,
                    facecolor=panel_figure.get_facecolor())
                panel_figure.savefig(
                    pdf_dir / f"{Path(spec.filename).stem}.pdf",
                    bbox_inches=None,
                    pad_inches=0.0,
                    transparent=True)
                plt.close(panel_figure)
                continue

            if spec.axis is None:
                raise ValueError(f"PanelExportSpec for {spec.filename!r} needs either axis or render_fn.")

            filename = spec.filename
            axis = spec.axis
            tight = spec.tight
        else:
            filename, axis, tight = spec

        bbox = axis.get_tightbbox(renderer) if tight else axis.get_window_extent(renderer)
        bbox = bbox.expanded(*pad_scale)
        extent = bbox.transformed(figure.dpi_scale_trans.inverted())
        figure.savefig(
            panel_dir / filename,
            dpi=300,
            bbox_inches=extent,
            facecolor=figure.get_facecolor())
        figure.savefig(
            pdf_dir / f"{Path(filename).stem}.pdf",
            bbox_inches=extent,
            transparent=True)
# %% OUTPUT ENSURERS
def ensure_two_color_output() -> Path:
    """Return the two-color fixed-alpha output, recreating it if required."""

    if TWO_COLOR_OUTPUT.exists():
        return TWO_COLOR_OUTPUT

    UNMIXED_DIR.mkdir(parents=True, exist_ok=True)
    return unmix(
        input_path=TWO_COLOR_INPUT,
        output_path=TWO_COLOR_OUTPUT,
        method="manual",
        alpha=0.62,
        alpha_mode="fixed",
        verbose=False,
    )


def ensure_two_color_method_outputs() -> dict[str, Path]:
    """Return reproducible two-channel outputs for all directed unmixing methods."""

    UNMIXED_DIR.mkdir(parents=True, exist_ok=True)

    method_specs: dict[str, dict[str, object]] = {
        "manual": {
            "output_path": TWO_COLOR_OUTPUT,
            "kwargs": {
                "method": "manual",
                "alpha": 0.62,
                "alpha_mode": "fixed"},
            "report_checks": {
                "method_effective": "manual",
                "alpha": 0.62},
        },
        "mean_ratio": {
            "output_path": UNMIXED_DIR / "2_color_unmixing_validation_unmixed_mean_ratio.tif",
            "kwargs": {
                "method": "mean_ratio"},
            "report_checks": {
                "method_effective": "mean_ratio"},
        },
        "linear_fit": {
            "output_path": UNMIXED_DIR / "2_color_unmixing_validation_unmixed_linear_fit.tif",
            "kwargs": {
                "method": "linear_fit"},
            "report_checks": {
                "method_effective": "linear_fit"},
        },
        "corr_min": {
            "output_path": UNMIXED_DIR / "2_color_unmixing_validation_unmixed_corr_min.tif",
            "kwargs": {
                "method": "corr_min"},
            "report_checks": {
                "method_effective": "corr_min"},
        },
        "mi_min": {
            "output_path": UNMIXED_DIR / "2_color_unmixing_validation_unmixed_mi_min.tif",
            "kwargs": {
                "method": "mi_min"},
            "report_checks": {
                "method_effective": "mi_min"},
        },
    }

    outputs: dict[str, Path] = {}
    for method_name, spec in method_specs.items():
        output_path = Path(spec["output_path"])
        report_path = output_path.with_suffix(output_path.suffix + ".json")
        needs_refresh = not output_path.exists() or not report_path.exists()

        if not needs_refresh:
            try:
                report = json.loads(report_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                needs_refresh = True
            else:
                for key, expected_value in dict(spec["report_checks"]).items():
                    if report.get(key) != expected_value:
                        needs_refresh = True
                        break

        if needs_refresh:
            output_path = unmix(
                input_path=TWO_COLOR_INPUT,
                output_path=output_path,
                verbose=False,
                **dict(spec["kwargs"]))

        outputs[method_name] = output_path

    return outputs

def ensure_five_color_outputs() -> tuple[Path, Path]:
    """Return the five-color ground-truth and MATLAB-N outputs, recreating them if needed."""

    if not FIVE_COLOR_GROUND_TRUTH_CONVERTED.exists():
        convert_time_encoded_stack_to_channel_stack(
            FIVE_COLOR_GROUND_TRUTH_INPUT,
            FIVE_COLOR_GROUND_TRUTH_CONVERTED)

    UNMIXED_DIR.mkdir(parents=True, exist_ok=True)
    matlab_n_path = FIVE_COLOR_MATLAB_N_OUTPUT
    if not _report_has_expected_value(
        matlab_n_path,
        key="pixel_bin_size",
        expected_value=8,
    ):
        matlab_n_path = unmix_picasso(
            input_path=FIVE_COLOR_INPUT,
            output_path=FIVE_COLOR_MATLAB_N_OUTPUT,
            channels=[0, 1, 2, 3, 4],
            implementation="matlab_n",
            background_percentile=1.0,
            mi_bins=64,
            alpha_max=1.0,
            max_iter=50,
            tolerance=1e-4,
            max_alpha_voxels=250_000,
            step_size=0.2,
            qn=100,
            pixel_bin_size=8,
            alpha_clip=0.5,
            random_state=42,
            clip_negative=True,
            output_dtype="float32",
            verbose=False)

    return FIVE_COLOR_GROUND_TRUTH_CONVERTED, matlab_n_path

def ensure_gfap_source_sink_outputs() -> tuple[Path, Path]:
    """Return the GFAP/LMNB1 MATLAB-N and source-sink outputs, recreating them if needed."""

    UNMIXED_DIR.mkdir(parents=True, exist_ok=True)

    matlab_n_path = GFAP_MATLAB_N_OUTPUT
    if not matlab_n_path.exists():
        matlab_n_path = unmix_picasso(
            input_path=GFAP_SOURCE_SINK_INPUT,
            output_path=GFAP_MATLAB_N_OUTPUT,
            channels=[0, 1],
            implementation="matlab_n",
            background_percentile=1.0,
            mi_bins=64,
            alpha_max=1.0,
            max_iter=50,
            tolerance=1e-4,
            max_alpha_voxels=250_000,
            step_size=0.2,
            qn=100,
            pixel_bin_size=16,
            alpha_clip=0.5,
            random_state=42,
            clip_negative=True,
            output_dtype="float32",
            verbose=False)

    source_sink_path = GFAP_SOURCE_SINK_OUTPUT
    if not source_sink_path.exists():
        source_sink_path = unmix_picasso(
            input_path=GFAP_SOURCE_SINK_INPUT,
            output_path=GFAP_SOURCE_SINK_OUTPUT,
            channels=[0, 1],
            implementation="source_sink_n",
            sink_channels=[0],
            neutral_channels=[],
            background_percentile=1.0,
            mi_bins=64,
            alpha_max=1.0,
            max_iter=50,
            tolerance=1e-4,
            max_alpha_voxels=250_000,
            source_sink_optimize_background=True,
            source_sink_max_background=0.2,
            source_sink_n_restarts=6,
            source_sink_joint_optimization=True,
            random_state=0,
            clip_negative=True,
            output_dtype="float32",
            verbose=False)

    return matlab_n_path, source_sink_path
# %% FIGURE GENERATORS
def create_two_color_example_figure() -> list[dict[str, object]]:
    """Create the two-channel biological example figure and return summary metrics."""

    output_path = ensure_two_color_output()
    method_output_paths = ensure_two_color_method_outputs()
    raw_stack, raw_metadata = load_stack_with_omio(TWO_COLOR_INPUT)
    corrected_stack, _ = load_stack_with_omio(output_path)
    method_stacks = {
        method_name: load_stack_with_omio(path)[0]
        for method_name, path in method_output_paths.items()}

    raw_channels = raw_stack[0, 0]
    corrected_channels = corrected_stack[0, 0]

    source_raw = raw_channels[0]
    target_raw = raw_channels[1]
    source_corrected = corrected_channels[0]
    target_corrected = corrected_channels[1]

    source_bounds = _compute_display_bounds(source_raw, low=0.5, high=99.8)
    target_bounds = _compute_display_bounds(
        np.concatenate(
            [np.asarray(target_raw, dtype=np.float32).ravel(),
             np.asarray(target_corrected, dtype=np.float32).ravel(),]),
        low=0.5,
        high=99.8)

    raw_corr_sum = _pairwise_abs_correlation_sum(raw_channels)
    corrected_corr_sum = _pairwise_abs_correlation_sum(corrected_channels)
    raw_mi_sum = _pairwise_mutual_information_sum(raw_channels)
    corrected_mi_sum = _pairwise_mutual_information_sum(corrected_channels)
    corr_values_all_methods = {
        "raw": raw_corr_sum,
        **{method_name: _pairwise_abs_correlation_sum(stack[0, 0])
           for method_name, stack in method_stacks.items()},}
    mi_values_all_methods = {
        "raw": raw_mi_sum,
        **{method_name: _pairwise_mutual_information_sum(stack[0, 0])
           for method_name, stack in method_stacks.items()},}
    pixel_size_x, pixel_unit = _metadata_pixel_size_x_and_unit(raw_metadata)
    figure_key = "figure_public_2color_example"

    figure = plt.figure(figsize=(13, 8), constrained_layout=True)
    grid = figure.add_gridspec(2, 4, wspace=0.15, hspace=0.18)
    panel_specs: list[PanelExportSpec | tuple[str, object, bool]] = []
    image_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "image_panels")
    corr_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "correlation_sum")
    mi_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "mutual_information")
    corr_methods_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "correlation_sum_all_methods")
    mi_methods_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "mutual_information_all_methods")
    negative_methods_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "negative_fraction_all_methods")

    method_alpha_values = {
        method_name: float(_load_sidecar_report(path)["alpha"])
        for method_name, path in method_output_paths.items()}
    negative_fraction_all_methods = {
        method_name: float(
            np.mean(
                (np.asarray(target_raw, dtype=np.float32)
                  - method_alpha_values[method_name] * np.asarray(source_raw, dtype=np.float32)
                ) < 0.0))
        for method_name in method_output_paths}

    def draw_two_color_panel(axis, panel_kind: str, *, standalone: bool = False) -> None:
        """Draw one panel for the external two-color example."""

        if panel_kind == "raw_source":
            axis.imshow(
                _colorize_channel(
                    source_raw,
                    DISPLAY_COLORS[0][1],
                    low_value=source_bounds[0],
                    high_value=source_bounds[1]))
            if not standalone:
                axis.set_title("raw source channel")
        elif panel_kind == "raw_target":
            axis.imshow(
                _colorize_channel(
                    target_raw,
                    DISPLAY_COLORS[1][1],
                    low_value=target_bounds[0],
                    high_value=target_bounds[1]))
            if not standalone:
                axis.set_title("raw target channel")
        elif panel_kind == "raw_composite":
            axis.imshow(
                _multichannel_composite(
                    raw_channels,
                    colors=[DISPLAY_COLORS[0][1], DISPLAY_COLORS[1][1]],
                    bounds=[source_bounds, target_bounds]))
            if not standalone:
                axis.set_title("raw composite")
        elif panel_kind == "corrected_target":
            axis.imshow(
                _colorize_channel(
                    target_corrected,
                    DISPLAY_COLORS[1][1],
                    low_value=target_bounds[0],
                    high_value=target_bounds[1]))
            if not standalone:
                axis.set_title("corrected target channel")
        elif panel_kind == "corrected_composite":
            axis.imshow(
                _multichannel_composite(
                    corrected_channels,
                    colors=[DISPLAY_COLORS[0][1], DISPLAY_COLORS[1][1]],
                    bounds=[source_bounds, target_bounds]))
            if not standalone:
                axis.set_title("corrected composite")
        elif panel_kind == "corr":
            axis.bar(
                ["raw", "corrected"],
                [raw_corr_sum, corrected_corr_sum],
                color=[PICASSO_PLOT_COLORS["mixed"], DIRECTED_METHOD_COLORS["manual"]],
                width=0.72)
            _apply_panel_axis_style(axis, corr_panel_config, default_title="correlation sum", default_ylabel="lower is better")
            return
        else:
            axis.bar(
                ["raw", "corrected"],
                [raw_mi_sum, corrected_mi_sum],
                color=[PICASSO_PLOT_COLORS["mixed"], DIRECTED_METHOD_COLORS["manual"]],
                width=0.72)
            _apply_panel_axis_style(axis, mi_panel_config, default_title="mutual information", default_ylabel="lower is better")
            return
        axis.axis("off")
        if standalone:
            _add_scalebar_to_image_axis(
                axis,
                source_raw.shape,
                image_panel_config,
                pixel_size_x=pixel_size_x,
                unit=pixel_unit)

    method_plot_labels = ["raw", "manual", "mean-ratio", "linear-fit", "corr-min", "MI-min"]
    method_plot_keys = ["raw", "manual", "mean_ratio", "linear_fit", "corr_min", "mi_min"]
    method_plot_colors = [
        PICASSO_PLOT_COLORS["mixed"],
        DIRECTED_METHOD_COLORS["manual"],
        DIRECTED_METHOD_COLORS["mean_ratio"],
        DIRECTED_METHOD_COLORS["linear_fit"],
        DIRECTED_METHOD_COLORS["corr_min"],
        DIRECTED_METHOD_COLORS["mi_min"]]

    def draw_all_method_corr_panel(axis) -> None:
        """Draw the all-method correlation-sum comparison for the two-color example."""

        axis.bar(
            method_plot_labels,
            [corr_values_all_methods[key] for key in method_plot_keys],
            color=method_plot_colors,
            width=0.72)
        _apply_panel_axis_style(
            axis,
            corr_methods_panel_config,
            default_title="correlation sum",
            default_ylabel="lower is better")

    def draw_all_method_mi_panel(axis) -> None:
        """Draw the all-method mutual-information comparison for the two-color example."""

        axis.bar(
            method_plot_labels,
            [mi_values_all_methods[key] for key in method_plot_keys],
            color=method_plot_colors,
            width=0.72)
        _apply_panel_axis_style(
            axis,
            mi_methods_panel_config,
            default_title="mutual information",
            default_ylabel="lower is better")

    def draw_all_method_negative_panel(axis) -> None:
        """Draw the fraction of negative pre-clipping target pixels across methods."""

        negative_labels = ["manual", "mean-ratio", "linear-fit", "corr-min", "MI-min"]
        negative_keys = ["manual", "mean_ratio", "linear_fit", "corr_min", "mi_min"]
        negative_colors = [
            DIRECTED_METHOD_COLORS["manual"],
            DIRECTED_METHOD_COLORS["mean_ratio"],
            DIRECTED_METHOD_COLORS["linear_fit"],
            DIRECTED_METHOD_COLORS["corr_min"],
            DIRECTED_METHOD_COLORS["mi_min"]]
        axis.bar(
            negative_labels,
            [negative_fraction_all_methods[key] for key in negative_keys],
            color=negative_colors,
            width=0.72)
        _apply_panel_axis_style(
            axis,
            negative_methods_panel_config,
            default_title="negative fraction before clipping",
            default_ylabel="fraction of target pixels")

    ax_source = figure.add_subplot(grid[0, 0])
    ax_target_raw = figure.add_subplot(grid[0, 1])
    ax_raw_composite = figure.add_subplot(grid[0, 2:4])
    ax_corr_metric = figure.add_subplot(grid[1, 0])
    ax_mi_metric = figure.add_subplot(grid[1, 1])
    ax_target_corrected = figure.add_subplot(grid[1, 2])
    ax_corrected_composite = figure.add_subplot(grid[1, 3])

    draw_two_color_panel(ax_source, "raw_source")
    draw_two_color_panel(ax_target_raw, "raw_target")
    draw_two_color_panel(ax_raw_composite, "raw_composite")
    draw_two_color_panel(ax_target_corrected, "corrected_target")
    draw_two_color_panel(ax_corrected_composite, "corrected_composite")
    draw_two_color_panel(ax_corr_metric, "corr")
    draw_two_color_panel(ax_mi_metric, "mi")

    figure.suptitle("Two-channel biological microscopy example", fontsize=15)
    figure_path = FIGURE_DIR / "figure_public_2color_example.png"
    _save_axis_panels(
        figure,
        [
            PanelExportSpec("panel_raw_source_channel.png", render_fn=lambda axis: draw_two_color_panel(axis, "raw_source", standalone=True), figsize=tuple(image_panel_config["figsize_cm"]), full_bleed=True),
            PanelExportSpec("panel_raw_target_channel.png", render_fn=lambda axis: draw_two_color_panel(axis, "raw_target", standalone=True), figsize=tuple(image_panel_config["figsize_cm"]), full_bleed=True),
            PanelExportSpec("panel_raw_composite.png", render_fn=lambda axis: draw_two_color_panel(axis, "raw_composite", standalone=True), figsize=tuple(image_panel_config["figsize_cm"]), full_bleed=True),
            PanelExportSpec(str(corr_panel_config["savename"]), render_fn=lambda axis: draw_two_color_panel(axis, "corr"), figsize=tuple(corr_panel_config["figsize_cm"])),
            PanelExportSpec(str(mi_panel_config["savename"]), render_fn=lambda axis: draw_two_color_panel(axis, "mi"), figsize=tuple(mi_panel_config["figsize_cm"])),
            PanelExportSpec(str(corr_methods_panel_config["savename"]), render_fn=draw_all_method_corr_panel, figsize=tuple(corr_methods_panel_config["figsize_cm"])),
            PanelExportSpec(str(mi_methods_panel_config["savename"]), render_fn=draw_all_method_mi_panel, figsize=tuple(mi_methods_panel_config["figsize_cm"])),
            PanelExportSpec(str(negative_methods_panel_config["savename"]), render_fn=draw_all_method_negative_panel, figsize=tuple(negative_methods_panel_config["figsize_cm"])),
            PanelExportSpec("panel_corrected_target_channel.png", render_fn=lambda axis: draw_two_color_panel(axis, "corrected_target", standalone=True), figsize=tuple(image_panel_config["figsize_cm"]), full_bleed=True),
            PanelExportSpec("panel_corrected_composite.png", render_fn=lambda axis: draw_two_color_panel(axis, "corrected_composite", standalone=True), figsize=tuple(image_panel_config["figsize_cm"]), full_bleed=True),
        ],
        figure_path=figure_path)
    figure.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(figure)

    return [
        {
            "example": "2_color_validation",
            "variant": "raw",
            "n_channels": 2,
            "pairwise_abs_corr_sum": _pairwise_abs_correlation_sum(raw_channels),
            "pairwise_mi_sum": _pairwise_mutual_information_sum(raw_channels),
        },
        {
            "example": "2_color_validation",
            "variant": "fixed_alpha",
            "n_channels": 2,
            "pairwise_abs_corr_sum": _pairwise_abs_correlation_sum(corrected_channels),
            "pairwise_mi_sum": _pairwise_mutual_information_sum(corrected_channels),
        }]

def create_five_color_example_figure() -> list[dict[str, object]]:
    """Create the five-channel simulation example figure and return summary metrics."""

    ground_truth_path, matlab_n_path = ensure_five_color_outputs()

    raw_stack, raw_metadata = load_stack_with_omio(FIVE_COLOR_INPUT)
    ground_truth_stack, _ = load_stack_with_omio(ground_truth_path)
    matlab_n_stack, _ = load_stack_with_omio(matlab_n_path)

    raw_channels = raw_stack[0, 0]
    ground_truth_channels = ground_truth_stack[0, 0]
    matlab_n_channels = matlab_n_stack[0, 0]

    stacks = [
        ("ground truth", ground_truth_channels),
        ("measured", raw_channels),
        ("MATLAB-style N", matlab_n_channels)]

    stack_bounds: list[list[tuple[float, float]]] = []
    for _, stack in stacks:
        channel_bounds: list[tuple[float, float]] = []
        for channel_index in range(5):
            channel_values = np.asarray(stack[channel_index], dtype=np.float32)
            channel_bounds.append(
                (float(np.percentile(channel_values, 0.5)),
                 float(np.percentile(channel_values, 99.8))))
        stack_bounds.append(channel_bounds)

    corr_values = [
        _pairwise_abs_correlation_sum(ground_truth_channels),
        _pairwise_abs_correlation_sum(raw_channels),
        _pairwise_abs_correlation_sum(matlab_n_channels)]
    mi_values = [
        _pairwise_mutual_information_sum(ground_truth_channels),
        _pairwise_mutual_information_sum(raw_channels),
        _pairwise_mutual_information_sum(matlab_n_channels)]
    recovery_values = [
        _best_match_channel_correlation(ground_truth_channels, ground_truth_channels),
        _best_match_channel_correlation(raw_channels, ground_truth_channels),
        _best_match_channel_correlation(matlab_n_channels, ground_truth_channels)]
    direct_corr_raw = _direct_channel_correlation_vector(raw_channels, ground_truth_channels)
    direct_corr_matlab_n = _direct_channel_correlation_vector(matlab_n_channels, ground_truth_channels)
    pixel_size_x, pixel_unit = _metadata_pixel_size_x_and_unit(raw_metadata)
    figure_key = "figure_public_5color_example"

    figure = plt.figure(figsize=(18, 12), constrained_layout=True)
    grid = figure.add_gridspec(4, 8, height_ratios=[1.0, 1.0, 1.0, 0.9])
    panel_specs: list[PanelExportSpec | tuple[str, object, bool]] = []
    image_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "image_panels")
    corr_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "correlation_sum")
    mi_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "mutual_information_sum")
    recovery_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "gt_matched_recovery")
    direct_corr_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "direct_gt_channel_correlation")

    def draw_five_color_panel(axis, stack: np.ndarray, row_bounds: list[tuple[float, float]], channel_index: int | None, *, standalone: bool = False) -> None:
        """Draw one channel or composite panel for the external five-color example."""

        if channel_index is None:
            axis.imshow(
                _multichannel_composite(
                    stack,
                    colors=[item[1] for item in DISPLAY_COLORS[:5]],
                    bounds=row_bounds,
                    gamma=0.75))
            if not standalone:
                axis.set_title("composite")
        else:
            axis.imshow(
                _colorize_channel(
                    stack[channel_index],
                    DISPLAY_COLORS[channel_index][1],
                    low_value=row_bounds[channel_index][0],
                    high_value=row_bounds[channel_index][1],
                    gamma=0.75))
            if not standalone:
                axis.set_title(f"channel {channel_index}")
        axis.axis("off")
        if standalone:
            _add_scalebar_to_image_axis(
                axis,
                stack[0].shape,
                image_panel_config,
                pixel_size_x=pixel_size_x,
                unit=pixel_unit)

    def draw_summary_barplot(
        axis,
        title: str,
        ylabel: str,
        labels: list[str],
        values: list[float],
        colors: list[str],
        *,
        ylim: tuple[float, float] | None = None,
        reference: float | None = None,
    ) -> None:
        """Draw one standalone summary bar plot."""

        axis.bar(labels, values, color=colors, width=0.72)
        if reference is not None:
            axis.axhline(
                float(reference),
                color=PICASSO_PLOT_COLORS["ground_truth"],
                linestyle="--",
                linewidth=1.4)
        panel_config = {
            "correlation sum": corr_panel_config,
            "mutual information sum": mi_panel_config,
            "GT-matched recovery": recovery_panel_config}[title]
        effective_config = dict(panel_config)
        if ylim is not None and effective_config.get("ylim") is None:
            effective_config["ylim"] = ylim
        _apply_panel_axis_style(axis, effective_config, default_title=title, default_ylabel=ylabel)

    def draw_direct_corr_panel(axis) -> None:
        """Draw the grouped direct channel-correlation summary for the external example."""

        channel_positions = np.arange(5)
        width = 0.36
        axis.bar(
            channel_positions - width / 2,
            direct_corr_raw,
            width=width,
            color=PICASSO_PLOT_COLORS["mixed"],
            label="measured")
        axis.bar(
            channel_positions + width / 2,
            direct_corr_matlab_n,
            width=width,
            color=PICASSO_PLOT_COLORS["matlab_n"],
            label="MATLAB-N")
        axis.axhline(1.0, color=PICASSO_PLOT_COLORS["ground_truth"], linestyle="--", linewidth=1.4, label="GT reference")
        axis.set_xticks(channel_positions, [str(index) for index in range(5)])
        axis.set_xlabel("channel")
        effective_config = dict(direct_corr_panel_config)
        if effective_config.get("ylim") is None:
            ymin = float(np.min([direct_corr_raw.min(), direct_corr_matlab_n.min()]))
            effective_config["ylim"] = (0.0, 1.05) if ymin >= 0.0 else (min(-0.1, ymin - 0.05), 1.05)
        _apply_panel_axis_style(axis, effective_config, default_title="direct GT channel correlation", default_ylabel="Pearson r")

    for row_index, (row_label, stack) in enumerate(stacks):
        row_slug = row_label.lower().replace(" ", "_").replace("-", "_")
        row_bounds = stack_bounds[row_index]
        for channel_index in range(5):
            axis = figure.add_subplot(grid[row_index, channel_index])
            draw_five_color_panel(axis, stack, row_bounds, channel_index)
            if row_index == 0:
                axis.set_title(f"channel {channel_index}")
            else:
                axis.set_title("")
            if channel_index == 0:
                axis.text(
                    -0.22,
                    0.5,
                    row_label,
                    transform=axis.transAxes,
                    rotation=90,
                    va="center",
                    ha="center",
                    fontsize=11)
            panel_specs.append(
                PanelExportSpec(
                    f"panel_{row_slug}_channel_{channel_index}.png",
                    render_fn=lambda axis, stack=stack.copy(), row_bounds=list(row_bounds), channel_index=channel_index: draw_five_color_panel(axis, stack, row_bounds, channel_index, standalone=True),
                    figsize=tuple(image_panel_config["figsize_cm"]),
                    full_bleed=True))

        composite_axis = figure.add_subplot(grid[row_index, 5])
        draw_five_color_panel(composite_axis, stack, row_bounds, None)
        if row_index == 0:
            composite_axis.set_title("composite")
        else:
            composite_axis.set_title("")
        panel_specs.append(
            PanelExportSpec(
                f"panel_{row_slug}_composite.png",
                render_fn=lambda axis, stack=stack.copy(), row_bounds=list(row_bounds): draw_five_color_panel(axis, stack, row_bounds, None, standalone=True),
                figsize=tuple(image_panel_config["figsize_cm"]),
                full_bleed=True))

    corr_axis = figure.add_subplot(grid[3, 0:2])
    draw_summary_barplot(
        corr_axis,
        "correlation sum",
        "lower is better",
        ["GT", "measured", "MATLAB-N"],
        corr_values,
        [PICASSO_PLOT_COLORS["ground_truth"], PICASSO_PLOT_COLORS["mixed"], PICASSO_PLOT_COLORS["matlab_n"]])

    mi_axis = figure.add_subplot(grid[3, 2:4])
    draw_summary_barplot(
        mi_axis,
        "mutual information sum",
        "lower is better",
        ["GT", "measured", "MATLAB-N"],
        mi_values,
        [PICASSO_PLOT_COLORS["ground_truth"], PICASSO_PLOT_COLORS["mixed"], PICASSO_PLOT_COLORS["matlab_n"]])

    recovery_axis = figure.add_subplot(grid[3, 4:6])
    draw_summary_barplot(
        recovery_axis,
        "GT-matched recovery",
        "higher is better",
        ["measured", "MATLAB-N"],
        recovery_values[1:],
        [PICASSO_PLOT_COLORS["mixed"], PICASSO_PLOT_COLORS["matlab_n"]],
        ylim=(0.0, 1.05),
        reference=1.0)

    per_channel_axis = figure.add_subplot(grid[3, 6:8])
    draw_direct_corr_panel(per_channel_axis)

    figure.suptitle("Five-channel PICASSO simulation image", fontsize=15)
    figure_path = FIGURE_DIR / "figure_public_5color_example.png"
    panel_specs.extend(
        [
            PanelExportSpec(
                str(corr_panel_config["savename"]),
                render_fn=lambda axis: draw_summary_barplot(
                    axis,
                    "correlation sum",
                    "lower is better",
                    ["GT", "measured", "MATLAB-N"],
                    corr_values,
                    [PICASSO_PLOT_COLORS["ground_truth"], PICASSO_PLOT_COLORS["mixed"], PICASSO_PLOT_COLORS["matlab_n"]],
                ), figsize=tuple(corr_panel_config["figsize_cm"])),
            PanelExportSpec(
                str(mi_panel_config["savename"]),
                render_fn=lambda axis: draw_summary_barplot(
                    axis,
                    "mutual information sum",
                    "lower is better",
                    ["GT", "measured", "MATLAB-N"],
                    mi_values,
                    [PICASSO_PLOT_COLORS["ground_truth"], PICASSO_PLOT_COLORS["mixed"], PICASSO_PLOT_COLORS["matlab_n"]],
                ), figsize=tuple(mi_panel_config["figsize_cm"])),
            PanelExportSpec(
                str(recovery_panel_config["savename"]),
                render_fn=lambda axis: draw_summary_barplot(
                    axis,
                    "GT-matched recovery",
                    "higher is better",
                    ["measured", "MATLAB-N"],
                    recovery_values[1:],
                    [PICASSO_PLOT_COLORS["mixed"], PICASSO_PLOT_COLORS["matlab_n"]],
                    ylim=(0.0, 1.05),
                    reference=1.0), figsize=tuple(recovery_panel_config["figsize_cm"])),
            PanelExportSpec(
                str(direct_corr_panel_config["savename"]),
                render_fn=draw_direct_corr_panel,
                figsize=tuple(direct_corr_panel_config["figsize_cm"])),
        ])
    _save_axis_panels(figure, panel_specs, figure_path=figure_path)
    figure.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(figure)

    return [
        {
            "example": "5_color_simulation",
            "variant": "ground_truth",
            "n_channels": 5,
            "pairwise_abs_corr_sum": _pairwise_abs_correlation_sum(ground_truth_channels),
            "pairwise_mi_sum": _pairwise_mutual_information_sum(ground_truth_channels),
            "best_match_corr_to_ground_truth": _best_match_channel_correlation(
                ground_truth_channels,
                ground_truth_channels),
        },
        {
            "example": "5_color_simulation",
            "variant": "raw",
            "n_channels": 5,
            "pairwise_abs_corr_sum": _pairwise_abs_correlation_sum(raw_channels),
            "pairwise_mi_sum": _pairwise_mutual_information_sum(raw_channels),
            "best_match_corr_to_ground_truth": _best_match_channel_correlation(
                raw_channels,
                ground_truth_channels),
            **{f"direct_corr_channel_{channel_index}": float(direct_corr_raw[channel_index])
               for channel_index in range(5)},
        },
        {
            "example": "5_color_simulation",
            "variant": "matlab_n",
            "n_channels": 5,
            "pairwise_abs_corr_sum": _pairwise_abs_correlation_sum(matlab_n_channels),
            "pairwise_mi_sum": _pairwise_mutual_information_sum(matlab_n_channels),
            "best_match_corr_to_ground_truth": _best_match_channel_correlation(
                matlab_n_channels,
                ground_truth_channels),
            **{f"direct_corr_channel_{channel_index}": float(direct_corr_matlab_n[channel_index])
               for channel_index in range(5)},
        },
    ]

def create_gfap_source_sink_example_figure() -> list[dict[str, object]]:
    """Create the realistic two-channel source-sink example figure."""

    matlab_n_path, source_sink_path = ensure_gfap_source_sink_outputs()

    raw_stack, raw_metadata = load_stack_with_omio(GFAP_SOURCE_SINK_INPUT)
    matlab_n_stack, _ = load_stack_with_omio(matlab_n_path)
    source_sink_stack, _ = load_stack_with_omio(source_sink_path)

    raw_channels = raw_stack[0, 0]
    matlab_n_channels = matlab_n_stack[0, 0]
    source_sink_channels = source_sink_stack[0, 0]

    source_raw = raw_channels[1]
    sink_raw = raw_channels[0]
    matlab_n_sink = matlab_n_channels[0]
    source_sink_sink = source_sink_channels[0]

    source_bounds = _compute_display_bounds(source_raw, low=0.5, high=99.8)
    sink_bounds_raw = _compute_display_bounds(sink_raw, low=0.5, high=99.8)
    sink_bounds_matlab_n = _compute_display_bounds(matlab_n_sink, low=0.5, high=99.8)
    sink_bounds_source_sink = _compute_display_bounds(source_sink_sink, low=0.5, high=99.8)

    corr_values = [
        _pairwise_abs_correlation_sum(raw_channels),
        _pairwise_abs_correlation_sum(matlab_n_channels),
        _pairwise_abs_correlation_sum(source_sink_channels)]
    mi_values = [
        _pairwise_mutual_information_sum(raw_channels),
        _pairwise_mutual_information_sum(matlab_n_channels),
        _pairwise_mutual_information_sum(source_sink_channels)]
    pixel_size_x, pixel_unit = _metadata_pixel_size_x_and_unit(raw_metadata)
    figure_key = "figure_public_source_sink_example"

    figure = plt.figure(figsize=(13, 8), constrained_layout=True)
    grid = figure.add_gridspec(2, 4, wspace=0.15, hspace=0.18)
    panel_specs: list[PanelExportSpec | tuple[str, object, bool]] = []
    image_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "image_panels")
    corr_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "correlation_sum")
    mi_panel_config = _get_panel_params(EXAMPLE_PANEL_PARAMS, figure_key, "mutual_information")

    def draw_source_sink_panel(axis, panel_kind: str, *, standalone: bool = False) -> None:
        """Draw one panel for the directional source-sink example."""

        if panel_kind == "raw_source":
            axis.imshow(
                _colorize_channel(
                    source_raw,
                    DISPLAY_COLORS[1][1],
                    low_value=source_bounds[0],
                    high_value=source_bounds[1]))
            if not standalone:
                axis.set_title("raw source channel")
        elif panel_kind == "raw_sink":
            axis.imshow(
                _colorize_channel(
                    sink_raw,
                    DISPLAY_COLORS[0][1],
                    low_value=sink_bounds_raw[0],
                    high_value=sink_bounds_raw[1]))
            if not standalone:
                axis.set_title("raw sink channel")
        elif panel_kind == "raw_composite":
            axis.imshow(
                _multichannel_composite(
                    raw_channels,
                    colors=[DISPLAY_COLORS[0][1], DISPLAY_COLORS[1][1]],
                    bounds=[sink_bounds_raw, source_bounds]))
            if not standalone:
                axis.set_title("raw composite")
        elif panel_kind == "matlab_n_sink":
            axis.imshow(
                _colorize_channel(
                    matlab_n_sink,
                    DISPLAY_COLORS[0][1],
                    low_value=sink_bounds_matlab_n[0],
                    high_value=sink_bounds_matlab_n[1]))
            if not standalone:
                axis.set_title("MATLAB-N sink channel")
        elif panel_kind == "matlab_n_composite":
            axis.imshow(
                _multichannel_composite(
                    np.stack([matlab_n_sink, source_raw], axis=0),
                    colors=[DISPLAY_COLORS[0][1], DISPLAY_COLORS[1][1]],
                    bounds=[sink_bounds_matlab_n, source_bounds]))
            if not standalone:
                axis.set_title("MATLAB-N composite")
        elif panel_kind == "source_sink_n_sink":
            axis.imshow(
                _colorize_channel(
                    source_sink_sink,
                    DISPLAY_COLORS[0][1],
                    low_value=sink_bounds_source_sink[0],
                    high_value=sink_bounds_source_sink[1]))
            if not standalone:
                axis.set_title("source-sink-N sink channel")
        elif panel_kind == "source_sink_n_composite":
            axis.imshow(
                _multichannel_composite(
                    source_sink_channels,
                    colors=[DISPLAY_COLORS[0][1], DISPLAY_COLORS[1][1]],
                    bounds=[sink_bounds_source_sink, source_bounds]))
            if not standalone:
                axis.set_title("source-sink-N composite")
        elif panel_kind == "corr":
            axis.bar(
                ["raw", "MATLAB-N", "source-sink-N"],
                corr_values,
                color=[PICASSO_PLOT_COLORS["mixed"], PICASSO_PLOT_COLORS["matlab_n"], PICASSO_PLOT_COLORS["source_sink_n"]],
                width=0.72)
            _apply_panel_axis_style(axis, corr_panel_config, default_title="correlation sum", default_ylabel="lower is better")
            return
        else:
            axis.bar(
                ["raw", "MATLAB-N", "source-sink-N"],
                mi_values,
                color=[PICASSO_PLOT_COLORS["mixed"], PICASSO_PLOT_COLORS["matlab_n"], PICASSO_PLOT_COLORS["source_sink_n"]],
                width=0.72)
            _apply_panel_axis_style(axis, mi_panel_config, default_title="mutual information", default_ylabel="lower is better")
            return
        axis.axis("off")
        if standalone:
            _add_scalebar_to_image_axis(
                axis,
                source_raw.shape,
                image_panel_config,
                pixel_size_x=pixel_size_x,
                unit=pixel_unit)

    ax_source_raw = figure.add_subplot(grid[0, 0])
    ax_sink_raw = figure.add_subplot(grid[0, 1])
    ax_raw_composite = figure.add_subplot(grid[0, 2])
    ax_matlab_n_sink = figure.add_subplot(grid[1, 0])
    ax_source_sink_sink = figure.add_subplot(grid[1, 1])
    ax_source_sink_composite = figure.add_subplot(grid[1, 2])

    metrics_grid = grid[:, 3].subgridspec(2, 1, hspace=0.35)
    ax_corr_metric = figure.add_subplot(metrics_grid[0, 0])
    ax_mi_metric = figure.add_subplot(metrics_grid[1, 0])

    draw_source_sink_panel(ax_source_raw, "raw_source")
    draw_source_sink_panel(ax_sink_raw, "raw_sink")
    draw_source_sink_panel(ax_raw_composite, "raw_composite")
    draw_source_sink_panel(ax_matlab_n_sink, "matlab_n_sink")
    draw_source_sink_panel(ax_source_sink_sink, "source_sink_n_sink")
    draw_source_sink_panel(ax_source_sink_composite, "source_sink_n_composite")
    draw_source_sink_panel(ax_corr_metric, "corr")
    draw_source_sink_panel(ax_mi_metric, "mi")

    figure.suptitle("Directional source-sink example", fontsize=15)
    figure_path = FIGURE_DIR / "figure_public_source_sink_example.png"
    _save_axis_panels(
        figure,
        [
            PanelExportSpec("panel_raw_source_channel.png", render_fn=lambda axis: draw_source_sink_panel(axis, "raw_source", standalone=True), figsize=tuple(image_panel_config["figsize_cm"]), full_bleed=True),
            PanelExportSpec("panel_raw_sink_channel.png", render_fn=lambda axis: draw_source_sink_panel(axis, "raw_sink", standalone=True), figsize=tuple(image_panel_config["figsize_cm"]), full_bleed=True),
            PanelExportSpec("panel_raw_composite.png", render_fn=lambda axis: draw_source_sink_panel(axis, "raw_composite", standalone=True), figsize=tuple(image_panel_config["figsize_cm"]), full_bleed=True),
            PanelExportSpec("panel_matlab_n_sink_channel.png", render_fn=lambda axis: draw_source_sink_panel(axis, "matlab_n_sink", standalone=True), figsize=tuple(image_panel_config["figsize_cm"]), full_bleed=True),
            PanelExportSpec("panel_matlab_n_composite.png", render_fn=lambda axis: draw_source_sink_panel(axis, "matlab_n_composite", standalone=True), figsize=tuple(image_panel_config["figsize_cm"]), full_bleed=True),
            PanelExportSpec("panel_source_sink_n_sink_channel.png", render_fn=lambda axis: draw_source_sink_panel(axis, "source_sink_n_sink", standalone=True), figsize=tuple(image_panel_config["figsize_cm"]), full_bleed=True),
            PanelExportSpec("panel_source_sink_n_composite.png", render_fn=lambda axis: draw_source_sink_panel(axis, "source_sink_n_composite", standalone=True), figsize=tuple(image_panel_config["figsize_cm"]), full_bleed=True),
            PanelExportSpec(str(corr_panel_config["savename"]), render_fn=lambda axis: draw_source_sink_panel(axis, "corr"), figsize=tuple(corr_panel_config["figsize_cm"])),
            PanelExportSpec(str(mi_panel_config["savename"]), render_fn=lambda axis: draw_source_sink_panel(axis, "mi"), figsize=tuple(mi_panel_config["figsize_cm"])),
        ], figure_path=figure_path)
    figure.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(figure)

    return [
        {
            "example": "GFAP_sink_LMNB1_source",
            "variant": "raw",
            "n_channels": 2,
            "pairwise_abs_corr_sum": _pairwise_abs_correlation_sum(raw_channels),
            "pairwise_mi_sum": _pairwise_mutual_information_sum(raw_channels),
        },
        {
            "example": "GFAP_sink_LMNB1_source",
            "variant": "matlab_n",
            "n_channels": 2,
            "pairwise_abs_corr_sum": _pairwise_abs_correlation_sum(matlab_n_channels),
            "pairwise_mi_sum": _pairwise_mutual_information_sum(matlab_n_channels),
        },
        {
            "example": "GFAP_sink_LMNB1_source",
            "variant": "source_sink_n",
            "n_channels": 2,
            "pairwise_abs_corr_sum": _pairwise_abs_correlation_sum(source_sink_channels),
            "pairwise_mi_sum": _pairwise_mutual_information_sum(source_sink_channels),
        },
    ]

# %% MAIN FUNCTION
def main() -> None:
    """Generate all externally sourced example figures and their summary metrics."""

    _ensure_output_dirs()
    metric_rows: list[dict[str, object]] = []
    metric_rows.extend(create_two_color_example_figure())
    metric_rows.extend(create_five_color_example_figure())
    metric_rows.extend(create_gfap_source_sink_example_figure())
    _write_csv(RESULTS_DIR / "public_example_dependence_metrics.csv", metric_rows)

    print(f"Wrote example microscopy figures to: {FIGURE_DIR}")
    print(f"Wrote example microscopy metrics to: {RESULTS_DIR / 'public_example_dependence_metrics.csv'}")
# %% MAIN
if __name__ == "__main__":
    main()
# %% END
