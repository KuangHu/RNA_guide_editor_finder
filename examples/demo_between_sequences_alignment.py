"""
Demo script for finding alignments between two different sequences.

This script demonstrates how to use the ShortAlignmentFinder to identify
alignments between two distinct DNA sequences, such as:
- Flanking region vs non-coding region of transposon
- Upstream region vs downstream region
- Any two distinct genomic regions

Use cases:
- Finding guide RNA sequences shared between different genomic regions
- Identifying regulatory elements that appear in multiple contexts
- Comparing different parts of transposon structure

Author: Kuang Hu
Date: 2026-01-26
"""

import sys
import os
import json

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules.short_alignment_finder import (
    ShortAlignmentFinder,
    find_alignments_between_sequences
)


def example_1_basic_between_sequences():
    """Example 1: Basic usage - finding alignments between two sequences."""
    print("=" * 80)
    print("Example 1: Finding alignments between two sequences")
    print("=" * 80)

    # Flanking region (e.g., 50bp upstream of transposon)
    flanking_region = (
        "ATGCTAGCTAGC"      # Random upstream sequence
        "TTTAAACCCGGG"       # 12bp potential guide RNA
        "GCTAGCTAGCTA"       # More upstream
    )

    # Non-coding region within transposon
    noncoding_region = (
        "NNNNNNNNNNNN"       # Some transposon sequence
        "TTTAAACCCGGG"       # Same 12bp guide RNA
        "NNNNNNNNNNNN"       # More transposon sequence
    )

    print(f"\nFlanking region ({len(flanking_region)}bp):")
    print(f"  {flanking_region}")
    print(f"\nNon-coding region ({len(noncoding_region)}bp):")
    print(f"  {noncoding_region}\n")

    # Find alignments between the two regions
    finder = ShortAlignmentFinder(min_length=9, max_mismatches=0)
    results = finder.find_alignments_between(flanking_region, noncoding_region)

    print(f"Found {len(results)} alignment(s) between the two sequences:\n")

    for i, result in enumerate(results, 1):
        print(f"Alignment {i}:")
        print(f"  Position in flanking: {result['pos1']}")
        print(f"  Position in noncoding: {result['pos2']}")
        print(f"  Length: {result['length']}bp")
        print(f"  Sequence (flanking): {result['seq1']}")
        print(f"  Sequence (noncoding): {result['seq2']}")
        print(f"  Mismatches: {result['mismatches']}")
        print(f"  Orientation: {result['orientation']}")
        print()


def example_2_with_mismatches():
    """Example 2: Finding alignments with mismatches between sequences."""
    print("=" * 80)
    print("Example 2: Between-sequences alignment with mismatches")
    print("=" * 80)

    flanking = "ATGCATGC" + "ACGTACGTACGT" + "NNNN"
    noncoding = "NNNN" + "ACGTACGTTCGT" + "GCTA"  # 1 mismatch at position 8

    print(f"\nFlanking:  {flanking}")
    print(f"Noncoding: {noncoding}")
    print("Note: There's 1 mismatch in the 12bp potential match\n")

    # Try with no mismatches
    finder_strict = ShortAlignmentFinder(min_length=9, max_mismatches=0)
    results_strict = finder_strict.find_alignments_between(flanking, noncoding)

    print(f"With max_mismatches=0: Found {len(results_strict)} alignment(s)")
    for r in results_strict:
        if r['length'] >= 10:
            print(f"  - {r['length']}bp match")

    # Try with 1 mismatch allowed
    finder_relaxed = ShortAlignmentFinder(min_length=9, max_mismatches=1)
    results_relaxed = finder_relaxed.find_alignments_between(flanking, noncoding)

    print(f"\nWith max_mismatches=1: Found {len(results_relaxed)} alignment(s)")
    for r in results_relaxed:
        if r['length'] >= 10:
            print(f"  - {r['length']}bp match at positions {r['pos1']}, {r['pos2']}")
            print(f"    Flanking:  {r['seq1']}")
            print(f"    Noncoding: {r['seq2']}")
            print(f"    Mismatches: {r['mismatches']} at positions {r['mismatch_positions']}")
    print()


def example_3_inverted_repeats_between():
    """Example 3: Finding inverted repeats between two sequences."""
    print("=" * 80)
    print("Example 3: Inverted repeats between two sequences")
    print("=" * 80)

    # Flanking has forward motif
    flanking = "ATGCATGC" + "ACGTACGTAC" + "NNNN"
    # Noncoding has reverse complement: RC(ACGTACGTAC) = GTACGTACGT
    noncoding = "NNNN" + "GTACGTACGT" + "GCTA"

    print(f"\nFlanking:  {flanking}")
    print(f"Noncoding: {noncoding}")
    print("Note: Noncoding contains the reverse complement of the flanking motif\n")

    # Find only inverted repeats
    finder = ShortAlignmentFinder(min_length=9, check_forward=False,
                                  check_revcomp=True)
    results = finder.find_alignments_between(flanking, noncoding)

    print(f"Found {len(results)} inverted repeat(s):\n")

    for r in results:
        if r['length'] >= 10:
            print(f"Inverted repeat ({r['length']}bp):")
            print(f"  Flanking position {r['pos1']}: {r['seq1']}")
            print(f"  Noncoding position {r['pos2']}: {r['seq2']}")
            print(f"  Orientation: {r['orientation']}")
            print()


def example_4_realistic_transposon_scenario():
    """Example 4: Realistic IS110 transposon analysis scenario."""
    print("=" * 80)
    print("Example 4: IS110 transposon guide RNA identification")
    print("=" * 80)

    # Simulate a realistic scenario:
    # - Flanking region: 100bp upstream of transposon
    # - Non-coding region: internal non-coding region of transposon
    # - Looking for guide RNA pairs (typically 9-15bp)

    flanking_100bp = (
        "ATGCTAGCTAGCTAGCTAGCTAGC"    # 24bp upstream
        "TTTAAACCCGGGAAA"              # 15bp guide RNA candidate
        "GCTAGCTAGCTAGCTAGCTAGCTA"    # 24bp
        "NNNNNNNNNNNN"                 # 12bp
        "GTACGTACGTACGTAC"             # 16bp
        "ATGCATGC"                     # 8bp
    )  # Total: ~99bp

    noncoding_region = (
        "NNNNNNNNNNNNNNNN"             # 16bp transposon sequence
        "TTTAAACCCGGGAAA"              # 15bp guide RNA (perfect match)
        "NNNNNNNNNNNNNNNNNNNNNNNNNNN" # 27bp spacer
        "GTACGTACGTACGTAC"             # 16bp different sequence
        "NNNNNNNN"                     # 8bp
        "TTTGAACCCGGGAAA"              # 15bp similar but 1 mismatch
        "NNNNNNNN"                     # 8bp
    )

    print(f"\nFlanking region length: {len(flanking_100bp)}bp")
    print(f"Non-coding region length: {len(noncoding_region)}bp")
    print("\nSearching for guide RNA candidates...\n")

    # IS110 parameters
    finder = ShortAlignmentFinder(
        min_length=9,          # Minimum guide RNA length
        max_mismatches=1,      # Allow 1 mismatch for biological variation
        check_forward=True,
        check_revcomp=True     # Check both orientations
    )

    results = finder.find_alignments_between(flanking_100bp, noncoding_region)

    # Filter for high-quality matches
    high_quality = [r for r in results if r['length'] >= 12 and r['mismatches'] <= 1]

    print(f"Found {len(high_quality)} high-quality guide RNA candidate(s):\n")

    for i, result in enumerate(high_quality, 1):
        print(f"Candidate {i}:")
        print(f"  Length: {result['length']}bp")
        print(f"  Flanking position: {result['pos1']}")
        print(f"  Noncoding position: {result['pos2']}")
        print(f"  Orientation: {result['orientation']}")
        print(f"  Flanking sequence:  {result['seq1']}")
        print(f"  Noncoding sequence: {result['seq2']}")

        if result['mismatches'] > 0:
            print(f"  Mismatches: {result['mismatches']} at positions {result['mismatch_positions']}")
        else:
            print(f"  Perfect match!")

        print()

    print("Analysis complete!")
    print()


def example_5_multiple_regions_comparison():
    """Example 5: Comparing multiple sequence pairs."""
    print("=" * 80)
    print("Example 5: Batch comparison of multiple sequence pairs")
    print("=" * 80)

    # Multiple transposon elements with their flanking and non-coding regions
    transposons = [
        {
            'id': 'IS110_001',
            'flanking': 'ATGCATGC' + 'TTTAAACCCGGG' + 'GCTAGCTA',
            'noncoding': 'NNNNNNNN' + 'TTTAAACCCGGG' + 'NNNNNNNN'
        },
        {
            'id': 'IS110_002',
            'flanking': 'GCTAGCTA' + 'GTACGTACGTAC' + 'ATGCATGC',
            'noncoding': 'NNNNNNNN' + 'GTACGTACGTAC' + 'NNNNNNNN'
        },
        {
            'id': 'IS110_003',
            'flanking': 'ATGCATGC' + 'ACGTACGTACGT' + 'GCTAGCTA',
            'noncoding': 'NNNNNNNN' + 'ACGTGCGTACGT' + 'NNNNNNNN'  # 1 mismatch
        },
    ]

    print("\nAnalyzing multiple transposon elements...\n")

    finder = ShortAlignmentFinder(min_length=9, max_mismatches=1)

    for transposon in transposons:
        results = finder.find_alignments_between(
            transposon['flanking'],
            transposon['noncoding']
        )

        # Get best match
        if results:
            best = max(results, key=lambda x: (x['length'], -x['mismatches']))
            print(f"{transposon['id']}:")
            print(f"  Best match: {best['length']}bp with {best['mismatches']} mismatch(es)")
            print(f"  Sequence: {best['seq1']}")
        else:
            print(f"{transposon['id']}: No matches found")

    print()


def example_6_convenience_function():
    """Example 6: Using the convenience function."""
    print("=" * 80)
    print("Example 6: Convenience function for quick analysis")
    print("=" * 80)

    flanking = "ATGCACGTACGTACGT"
    noncoding = "NNNNACGTACGTACGT"

    print(f"\nFlanking:  {flanking}")
    print(f"Noncoding: {noncoding}\n")

    # One-line function call
    results = find_alignments_between_sequences(
        flanking,
        noncoding,
        min_length=9,
        max_mismatches=0
    )

    print(f"Found {len(results)} alignment(s) using convenience function:\n")

    for r in results:
        print(f"  - {r['length']}bp match at positions {r['pos1']}, {r['pos2']}")

    print()


def example_7_json_export():
    """Example 7: Exporting results for downstream analysis."""
    print("=" * 80)
    print("Example 7: JSON export for bioinformatics pipelines")
    print("=" * 80)

    flanking = "ATGCATGC" + "TTTAAACCCGGG" + "GCTAGCTA"
    noncoding = "NNNNNNNN" + "TTTAAACCCGGG" + "NNNNNNNN"

    results = find_alignments_between_sequences(flanking, noncoding, min_length=9)

    # Create a structured output for a pipeline
    pipeline_output = {
        'analysis': 'transposon_guide_rna_finder',
        'sequences': {
            'flanking_length': len(flanking),
            'noncoding_length': len(noncoding)
        },
        'parameters': {
            'min_length': 9,
            'max_mismatches': 0
        },
        'results': results
    }

    # Convert to JSON
    json_output = json.dumps(pipeline_output, indent=2)

    print("\nJSON output for downstream analysis:\n")
    print(json_output)
    print("\nThis format can be:")
    print("  - Saved to files for later analysis")
    print("  - Integrated into bioinformatics pipelines")
    print("  - Parsed by other tools and scripts")
    print()


def main():
    """Run all examples."""
    print("\n")
    print("*" * 80)
    print("BETWEEN-SEQUENCES ALIGNMENT FINDER - DEMONSTRATION")
    print("*" * 80)
    print()

    example_1_basic_between_sequences()
    example_2_with_mismatches()
    example_3_inverted_repeats_between()
    example_4_realistic_transposon_scenario()
    example_5_multiple_regions_comparison()
    example_6_convenience_function()
    example_7_json_export()

    print("*" * 80)
    print("All examples completed!")
    print("*" * 80)
    print()


if __name__ == '__main__':
    main()
