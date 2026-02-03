"""
Demo script for the Short Alignment Finder module.

This script demonstrates how to use the ShortAlignmentFinder to identify
direct repeats and inverted repeats in DNA sequences.

Use cases:
- Finding repeats in transposon upstream regions
- Identifying potential regulatory elements
- Analyzing repeat structure in genomic sequences

Author: Kuang Hu
Date: 2026-01-26
"""

import sys
import os
import json

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules.short_alignment_finder import ShortAlignmentFinder, find_short_alignments


def example_1_basic_usage():
    """Example 1: Basic usage with perfect direct repeats."""
    print("=" * 80)
    print("Example 1: Finding perfect direct repeats")
    print("=" * 80)

    # A sequence with a 15bp direct repeat separated by 20bp
    sequence = (
        "ATGCATGCATGCATGC"      # 16bp background
        "TTTAAACCCGGGAAA"        # 15bp repeat 1 (pos 16-30)
        "NNNNNNNNNNNNNNNNNNNNN"  # 21bp gap
        "TTTAAACCCGGGAAA"        # 15bp repeat 2 (pos 52-66)
        "GCTAGCTAGCTAGCTA"       # 16bp background
    )

    print(f"\nSequence length: {len(sequence)}bp")
    print(f"Sequence: {sequence}\n")

    # Find repeats with IS110-like parameters (min 9bp, no mismatches)
    finder = ShortAlignmentFinder(min_length=9, max_mismatches=0)
    results = finder.find_alignments(sequence)

    print(f"Found {len(results)} alignment(s):\n")

    for i, result in enumerate(results, 1):
        print(f"Alignment {i}:")
        print(f"  Position 1: {result['pos1']}")
        print(f"  Position 2: {result['pos2']}")
        print(f"  Distance: {result['distance']}bp")
        print(f"  Length: {result['length']}bp")
        print(f"  Sequence 1: {result['seq1']}")
        print(f"  Sequence 2: {result['seq2']}")
        print(f"  Mismatches: {result['mismatches']}")
        print(f"  Orientation: {result['orientation']}")
        print()


def example_2_with_mismatches():
    """Example 2: Finding repeats with mismatches allowed."""
    print("=" * 80)
    print("Example 2: Finding repeats with mismatches")
    print("=" * 80)

    # A 12bp repeat with 1 mismatch
    sequence = (
        "ATGCATGCATGC"
        "ACGTACGTACGT"  # 12bp repeat 1
        "NNNNNNNN"
        "ACGTACGTTCGT"  # 12bp repeat 2 (1 mismatch at position 8: A->T)
        "GCTAGCTAGCTA"
    )

    print(f"\nSequence: {sequence}")
    print("Note: There's 1 mismatch at position 8 (A->T)\n")

    # Try with no mismatches allowed
    finder_strict = ShortAlignmentFinder(min_length=9, max_mismatches=0)
    results_strict = finder_strict.find_alignments(sequence)

    print(f"With max_mismatches=0: Found {len(results_strict)} long alignment(s)")
    for r in results_strict:
        if r['length'] >= 10:
            print(f"  - {r['length']}bp at positions {r['pos1']}, {r['pos2']}")

    # Try with 1 mismatch allowed
    finder_relaxed = ShortAlignmentFinder(min_length=9, max_mismatches=1)
    results_relaxed = finder_relaxed.find_alignments(sequence)

    print(f"\nWith max_mismatches=1: Found {len(results_relaxed)} alignment(s)")
    for r in results_relaxed:
        if r['length'] >= 10:
            print(f"  - {r['length']}bp at positions {r['pos1']}, {r['pos2']}")
            print(f"    Mismatches: {r['mismatches']} at positions {r['mismatch_positions']}")
            print(f"    Seq1: {r['seq1']}")
            print(f"    Seq2: {r['seq2']}")
    print()


def example_3_inverted_repeats():
    """Example 3: Finding inverted repeats (reverse complement matches)."""
    print("=" * 80)
    print("Example 3: Finding inverted repeats")
    print("=" * 80)

    # ACGTACGTACGT has reverse complement ACGTACGTACGT (palindrome!)
    # Let's use a non-palindromic example
    sequence = (
        "ATGCATGCATGC"
        "ACGTACGTAC"    # Forward sequence
        "NNNNNNNN"
        "GTACGTACGT"    # Reverse complement of ACGTACGTAC
        "GCTAGCTAGCTA"
    )

    print(f"\nSequence: {sequence}")
    print("Note: Contains both direct and inverted repeats\n")

    # Find only direct repeats
    finder_forward = ShortAlignmentFinder(min_length=9, check_forward=True,
                                         check_revcomp=False)
    results_forward = finder_forward.find_alignments(sequence)

    print(f"Direct repeats only: {len(results_forward)} found")
    for r in results_forward:
        print(f"  - {r['length']}bp at positions {r['pos1']}, {r['pos2']} ({r['orientation']})")

    # Find only inverted repeats
    finder_inverted = ShortAlignmentFinder(min_length=9, check_forward=False,
                                          check_revcomp=True)
    results_inverted = finder_inverted.find_alignments(sequence)

    print(f"\nInverted repeats only: {len(results_inverted)} found")
    for r in results_inverted:
        print(f"  - {r['length']}bp at positions {r['pos1']}, {r['pos2']} ({r['orientation']})")
        print(f"    Seq1: {r['seq1']}")
        print(f"    Seq2: {r['seq2']} (appears as this in sequence)")

    # Find both
    finder_both = ShortAlignmentFinder(min_length=9, check_forward=True,
                                      check_revcomp=True)
    results_both = finder_both.find_alignments(sequence)

    print(f"\nBoth types: {len(results_both)} found")
    print()


def example_4_distance_constraints():
    """Example 4: Using distance constraints."""
    print("=" * 80)
    print("Example 4: Distance constraints")
    print("=" * 80)

    # Sequence with repeats at different distances
    sequence = (
        "AAAAAAAAAA"   # pos 0-9
        "NN"           # 2bp gap
        "AAAAAAAAAA"   # pos 12-21 (distance from pos 0 = 12)
        "NNNNNN"       # 6bp gap
        "AAAAAAAAAA"   # pos 28-37 (distance from pos 0 = 28)
        "NNNNNNNNNN"   # 10bp gap
        "AAAAAAAAAA"   # pos 48-57 (distance from pos 0 = 48)
    )

    print(f"\nSequence: {sequence}")
    print("Contains repeats at distances: 12bp, 28bp, 48bp from first repeat\n")

    # Only find nearby repeats (distance 10-30)
    finder_nearby = ShortAlignmentFinder(min_length=9, min_distance=10,
                                        max_distance=30)
    results_nearby = finder_nearby.find_alignments(sequence)

    print(f"With distance range [10, 30]: {len(results_nearby)} alignment(s)")
    for r in results_nearby:
        print(f"  - Positions {r['pos1']} to {r['pos2']}, distance = {r['distance']}bp")

    # Only find distant repeats (distance > 35)
    finder_distant = ShortAlignmentFinder(min_length=9, min_distance=35)
    results_distant = finder_distant.find_alignments(sequence)

    print(f"\nWith min_distance=35: {len(results_distant)} alignment(s)")
    for r in results_distant:
        print(f"  - Positions {r['pos1']} to {r['pos2']}, distance = {r['distance']}bp")
    print()


def example_5_auto_extension():
    """Example 5: Auto-extension feature."""
    print("=" * 80)
    print("Example 5: Auto-extension of alignments")
    print("=" * 80)

    # A long 20bp repeat
    sequence = (
        "ATGC"
        "ACGTACGTACGTACGTACGT"  # 20bp repeat 1
        "NNNN"
        "ACGTACGTACGTACGTACGT"  # 20bp repeat 2
        "GCTA"
    )

    print(f"\nSequence: {sequence}")
    print("Contains a 20bp perfect repeat\n")

    # Even though we set min_length=9, the algorithm should extend to 20bp
    finder = ShortAlignmentFinder(min_length=9, max_mismatches=0)
    results = finder.find_alignments(sequence)

    print(f"With min_length=9, found {len(results)} alignment(s):")
    for r in results:
        print(f"  - Length: {r['length']}bp (extended from minimum 9bp)")
        print(f"    Positions: {r['pos1']}, {r['pos2']}")
        print(f"    Sequence: {r['seq1']}")

    print("\nNote: The algorithm automatically extended the 9bp seed matches")
    print("      to the full 20bp alignment, and consolidated overlapping matches.")
    print()


def example_6_json_output():
    """Example 6: Exporting results as JSON."""
    print("=" * 80)
    print("Example 6: JSON output for downstream analysis")
    print("=" * 80)

    sequence = "ACGTACGTAC" + "NNNN" + "ACGTACGTAC"

    results = find_short_alignments(sequence, min_length=9)

    # Convert to JSON
    json_output = json.dumps(results, indent=2)

    print("\nJSON output:")
    print(json_output)
    print("\nThis format is ideal for:")
    print("  - Saving to files")
    print("  - Integration with other tools")
    print("  - Downstream analysis pipelines")
    print()


def example_7_is110_workflow():
    """Example 7: Typical IS110 analysis workflow."""
    print("=" * 80)
    print("Example 7: IS110 transposon analysis workflow")
    print("=" * 80)

    # Simulate an IS110 upstream region
    upstream_sequence = (
        "ATGCTAGCTAGCTAGCTAGC"      # Random background
        "TTTAAACCCGGG"               # 12bp repeat 1
        "NNNNNNNNNNNNNNNNNNNNNNNNN"  # 25bp spacer (typical IS110)
        "TTTAAACCCGGG"               # 12bp repeat 2
        "NNNNNNNNNNNN"               # More spacer
        "GTACGTACGTAC"               # 12bp repeat 3 (different sequence)
        "NNNNNNNNNNNNNNNNNNNNNNNNN"  # 25bp spacer
        "GTACGTACGTAC"               # 12bp repeat 4
        "GCTAGCTAGCTAGCTAGCTA"       # Random background
    )

    print("\nAnalyzing IS110 upstream region...")
    print(f"Sequence length: {len(upstream_sequence)}bp\n")

    # IS110-specific parameters
    finder = ShortAlignmentFinder(
        min_length=9,        # IS110 guide RNAs typically 9-15bp
        max_mismatches=1,    # Allow 1 mismatch for biological variation
        min_distance=20,     # Repeats typically separated by spacers
        max_distance=200,    # Within a reasonable genomic window
        check_forward=True,
        check_revcomp=True   # Check both orientations
    )

    results = finder.find_alignments(upstream_sequence)

    print(f"Found {len(results)} potential guide RNA pairs:\n")

    for i, result in enumerate(results, 1):
        print(f"Pair {i}:")
        print(f"  Positions: {result['pos1']} and {result['pos2']}")
        print(f"  Length: {result['length']}bp")
        print(f"  Distance: {result['distance']}bp")
        print(f"  Orientation: {result['orientation']}")
        print(f"  Sequence 1: {result['seq1']}")
        print(f"  Sequence 2: {result['seq2']}")

        if result['mismatches'] > 0:
            print(f"  Mismatches: {result['mismatches']} at positions {result['mismatch_positions']}")

        print()

    print("Analysis complete!")
    print()


def main():
    """Run all examples."""
    print("\n")
    print("*" * 80)
    print("SHORT ALIGNMENT FINDER - DEMONSTRATION")
    print("*" * 80)
    print()

    example_1_basic_usage()
    example_2_with_mismatches()
    example_3_inverted_repeats()
    example_4_distance_constraints()
    example_5_auto_extension()
    example_6_json_output()
    example_7_is110_workflow()

    print("*" * 80)
    print("All examples completed!")
    print("*" * 80)
    print()


if __name__ == '__main__':
    main()
