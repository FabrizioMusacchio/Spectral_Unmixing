"""
Command-line interface for spectral unmixing.

Author: Fabrizio Musacchio
Date: June 2026
"""
# %% IMPORTS
from __future__ import annotations

import argparse

from .unmixing import unmix

# %% CLI PARSER AND ENTRY POINT
def build_parser() -> argparse.ArgumentParser:
    """
    Build the command-line parser for the two-channel unmixing workflow.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser for the package CLI entry point.
    """

    parser = argparse.ArgumentParser(
        description="Remove bleed-through from one microscopy channel into another."
    )
    parser.add_argument("input_path", help="Input microscopy stack path readable by OMIO.")
    parser.add_argument("output_path", help="Output stack path to be written via OMIO.")
    parser.add_argument(
        "--alpha",
        type=float,
        default=None,
        help="Fixed bleed-through coefficient. Required for alpha_mode='fixed'.",
    )
    parser.add_argument(
        "--alpha-mode",
        choices=["fixed", "reference_t", "per_t"],
        default="fixed",
        help="How alpha should be obtained.",
    )
    parser.add_argument(
        "--alpha-reference-t",
        type=int,
        default=0,
        help="Reference time point used when alpha_mode='reference_t'.",
    )
    parser.add_argument(
        "--method",
        choices=["manual", "mean_ratio", "linear_fit", "corr_min", "mi_min"],
        default="mean_ratio",
        help="Method used to estimate alpha for non-fixed alpha modes.",
    )
    parser.add_argument(
        "--bidirectional",
        action="store_true",
        help="Enable bidirectional two-channel unmixing with a 2x2 mixing model.",
    )
    parser.add_argument(
        "--alpha-reverse",
        type=float,
        default=None,
        help="Optional reverse-direction fixed bleed-through coefficient.",
    )
    parser.add_argument(
        "--method-reverse",
        choices=["manual", "mean_ratio", "linear_fit", "corr_min", "mi_min"],
        default=None,
        help="Optional reverse-direction alpha-estimation method.",
    )
    parser.add_argument(
        "--source-channel",
        type=int,
        default=0,
        help="Channel that bleeds into the target channel.",
    )
    parser.add_argument(
        "--target-channel",
        type=int,
        default=1,
        help="Channel to be corrected.",
    )
    parser.add_argument(
        "--signal-percentile",
        type=float,
        default=99.0,
        help="Percentile threshold for the bright source-signal mask.",
    )
    parser.add_argument(
        "--target-low-percentile",
        type=float,
        default=None,
        help="Optional low-target percentile used to restrict the alpha mask.",
    )
    parser.add_argument(
        "--background-percentile",
        type=float,
        default=1.0,
        help="Percentile used for rough background subtraction.",
    )
    parser.add_argument(
        "--signal-percentile-reverse",
        type=float,
        default=None,
        help="Optional reverse-direction source-mask percentile.",
    )
    parser.add_argument(
        "--target-low-percentile-reverse",
        type=float,
        default=None,
        help="Optional reverse-direction low-target percentile.",
    )
    parser.add_argument(
        "--background-percentile-reverse",
        type=float,
        default=None,
        help="Optional reverse-direction background percentile.",
    )
    parser.add_argument(
        "--alpha-max",
        type=float,
        default=1.0,
        help="Upper bound used by optimization-based alpha methods.",
    )
    parser.add_argument(
        "--alpha-max-reverse",
        type=float,
        default=None,
        help="Optional reverse-direction optimization bound.",
    )
    parser.add_argument(
        "--mi-bins",
        type=int,
        default=64,
        help="Number of histogram bins used for mutual-information estimation.",
    )
    parser.add_argument(
        "--mi-bins-reverse",
        type=int,
        default=None,
        help="Optional reverse-direction histogram bin count for mutual information.",
    )
    parser.add_argument(
        "--max-alpha-voxels",
        type=int,
        default=500000,
        help="Maximum number of voxels used for alpha estimation after subsampling.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=0,
        help="Random seed used when alpha-estimation voxels are subsampled.",
    )
    parser.add_argument(
        "--no-preprocess-alpha-inputs",
        action="store_true",
        help="Disable optional background subtraction and clipping for alpha estimation.",
    )
    parser.add_argument(
        "--no-clip-negative",
        action="store_true",
        help="Disable clipping of negative corrected target values.",
    )
    parser.add_argument(
        "--output-dtype",
        default="float32",
        help="Output dtype used when saving the corrected stack.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress terminal progress output.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """
    Run the package command-line interface.

    Parameters
    ----------
    argv : list of str or None, optional
        Optional argument vector. If ``None``, arguments are taken from
        ``sys.argv``.

    Returns
    -------
    int
        Process-style exit code. ``0`` indicates success.
    """

    parser = build_parser()
    args = parser.parse_args(argv)
    unmix(
        input_path=args.input_path,
        output_path=args.output_path,
        alpha=args.alpha,
        alpha_mode=args.alpha_mode,
        alpha_reference_t=args.alpha_reference_t,
        method=args.method,
        bidirectional=args.bidirectional,
        alpha_reverse=args.alpha_reverse,
        method_reverse=args.method_reverse,
        source_channel=args.source_channel,
        target_channel=args.target_channel,
        signal_percentile=args.signal_percentile,
        target_low_percentile=args.target_low_percentile,
        background_percentile=args.background_percentile,
        signal_percentile_reverse=args.signal_percentile_reverse,
        target_low_percentile_reverse=args.target_low_percentile_reverse,
        background_percentile_reverse=args.background_percentile_reverse,
        preprocess_alpha_inputs=not args.no_preprocess_alpha_inputs,
        alpha_max=args.alpha_max,
        alpha_max_reverse=args.alpha_max_reverse,
        mi_bins=args.mi_bins,
        mi_bins_reverse=args.mi_bins_reverse,
        max_alpha_voxels=args.max_alpha_voxels,
        random_state=args.random_state,
        clip_negative=not args.no_clip_negative,
        output_dtype=args.output_dtype,
        verbose=not args.quiet,
    )
    return 0
# %% END
