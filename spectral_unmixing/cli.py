"""
Command-line interface for spectral unmixing.

Author: Fabrizio Musacchio
Date: June 2026
"""

from __future__ import annotations

import argparse

from .unmixing import unmix


def build_parser() -> argparse.ArgumentParser:
    """Build the package CLI parser."""

    parser = argparse.ArgumentParser(
        description="Remove bleed-through from one microscopy channel into another."
    )
    parser.add_argument("input_path", help="Input TIFF stack path.")
    parser.add_argument("output_path", help="Output TIFF stack path.")
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
        "--background-percentile",
        type=float,
        default=1.0,
        help="Percentile used for rough background subtraction.",
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
    """Run the spectral unmixing CLI."""

    parser = build_parser()
    args = parser.parse_args(argv)
    unmix(
        input_path=args.input_path,
        output_path=args.output_path,
        alpha=args.alpha,
        alpha_mode=args.alpha_mode,
        alpha_reference_t=args.alpha_reference_t,
        source_channel=args.source_channel,
        target_channel=args.target_channel,
        signal_percentile=args.signal_percentile,
        background_percentile=args.background_percentile,
        clip_negative=not args.no_clip_negative,
        output_dtype=args.output_dtype,
        verbose=not args.quiet,
    )
    return 0
