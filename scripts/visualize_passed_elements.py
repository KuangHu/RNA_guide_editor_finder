#!/usr/bin/env python3
"""
Batch visualization of IS elements with passed alignments.

Reads passed_alignments_annotated.json, groups by IS element, loads the full
element from samples/{sample}/is_extraction/is_elements.json, and generates
per-element PNG diagrams.

Usage:
    python scripts/visualize_passed_elements.py \
        --annotated-alignments /path/to/passed_alignments_annotated.json \
        --samples-dir /path/to/samples/ \
        --output-dir /path/to/output_pngs/ \
        --max-length 2000

Author: Kuang Hu
Date: 2026-02-10
"""

import argparse
import json
import logging
import os
import sys
from collections import defaultdict

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from modules.is_element_visualizer import ISElementVisualizer
from modules.is_element_genbank import ISElementGenBank

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)


def load_annotated_alignments(path: str) -> dict:
    """Load passed_alignments_annotated.json."""
    with open(path) as f:
        return json.load(f)


def group_alignments_by_is_id(alignments: list) -> dict:
    """Group alignment dicts by is_id."""
    grouped = defaultdict(list)
    for aln in alignments:
        grouped[aln["is_id"]].append(aln)
    return dict(grouped)


def load_element_from_sample(is_id: str, sample: str, samples_dir: str) -> dict | None:
    """Load a specific IS element from its sample's is_elements.json."""
    json_path = os.path.join(samples_dir, sample, "is_extraction", "is_elements.json")
    if not os.path.exists(json_path):
        logger.warning("Missing is_elements.json for sample %s", sample)
        return None
    with open(json_path) as f:
        elements = json.load(f)
    for el in elements:
        if el.get("is_id") == is_id:
            return el
    logger.warning("Element %s not found in %s", is_id, json_path)
    return None


def collect_orf_annotations(alignments: list) -> list:
    """Collect unique ORF annotations from alignment list."""
    seen = set()
    annotations = []
    for aln in alignments:
        for ann in aln.get("orf_annotations", []):
            orf_id = ann.get("orf_id", "")
            if orf_id not in seen:
                seen.add(orf_id)
                annotations.append(ann)
    return annotations


def main():
    parser = argparse.ArgumentParser(
        description="Generate PNG diagrams for IS elements with passed alignments",
    )
    parser.add_argument(
        "--annotated-alignments", required=True,
        help="Path to passed_alignments_annotated.json",
    )
    parser.add_argument(
        "--samples-dir", required=True,
        help="Path to samples directory containing {sample}/is_extraction/is_elements.json",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Directory to write PNG files",
    )
    parser.add_argument(
        "--max-length", type=int, default=2000,
        help="Maximum IS element length to visualize (default: 2000)",
    )
    parser.add_argument(
        "--dpi", type=int, default=150,
        help="PNG resolution (default: 150)",
    )
    args = parser.parse_args()

    # Load alignments
    logger.info("Loading annotated alignments from %s", args.annotated_alignments)
    data = load_annotated_alignments(args.annotated_alignments)
    alignments = data["alignments"]
    logger.info("Loaded %d passed alignments", len(alignments))

    # Group by IS element
    grouped = group_alignments_by_is_id(alignments)
    logger.info("Found %d unique IS elements with passed alignments", len(grouped))

    # Prepare output
    os.makedirs(args.output_dir, exist_ok=True)
    visualizer = ISElementVisualizer()
    genbank_exporter = ISElementGenBank()

    # Process each element
    processed = 0
    skipped_length = 0
    skipped_missing = 0

    for is_id, alns in sorted(grouped.items()):
        sample = alns[0]["sample"]

        # Load full element data
        element = load_element_from_sample(is_id, sample, args.samples_dir)
        if element is None:
            skipped_missing += 1
            continue

        # Check length filter
        is_length = element.get("is_element", {}).get("length", 0)
        if is_length > args.max_length:
            skipped_length += 1
            continue

        # Collect ORF annotations from alignments
        orf_annotations = collect_orf_annotations(alns)

        # Generate PNG
        png_path = os.path.join(args.output_dir, f"{is_id}.png")
        try:
            visualizer.save_element_png(
                element, alns, png_path,
                dpi=args.dpi,
                orf_annotations=orf_annotations or None,
            )
        except Exception:
            logger.exception("Failed to visualize %s", is_id)

        # Generate GenBank
        gbk_path = os.path.join(args.output_dir, f"{is_id}.gbk")
        try:
            genbank_exporter.save_genbank(
                element, alns, gbk_path,
                orf_annotations=orf_annotations or None,
            )
        except Exception:
            logger.exception("Failed to export GenBank for %s", is_id)

        processed += 1
        if processed % 50 == 0:
            logger.info("Progress: %d elements processed", processed)

    logger.info(
        "Done. Visualized: %d, Skipped (length > %d): %d, Skipped (missing data): %d",
        processed, args.max_length, skipped_length, skipped_missing,
    )


if __name__ == "__main__":
    main()
