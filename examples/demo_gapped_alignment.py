#!/usr/bin/env python3
"""
Demo: Gapped Alignment Extension for CAST Homing Spacer Detection

This script demonstrates how to use the gapped alignment extension feature
to find CAST homing spacer-like patterns (e.g., 17/20, 27/32 matches with gaps).

The workflow:
1. Find ungapped seed alignments (exact or near-exact matches)
2. Extend seeds with gaps allowed
3. Filter by identity threshold (e.g., >80%, >85%)

Author: Kuang Hu
Date: 2026-01-27
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules.short_alignment_finder import (
    ShortAlignmentFinder,
    find_alignments_between_sequences,
    extend_alignment_with_gaps
)


def demo_basic_gapped_extension():
    """Demo 1: Basic gapped alignment extension"""
    print("=" * 70)
    print("DEMO 1: Basic Gapped Alignment Extension")
    print("=" * 70)

    # Sequences with a 3bp insertion in seq2
    seq1 = "ACGTACGTACGTACGT"
    seq2 = "ACGTACGTNNNACGTACGT"

    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print()

    # Step 1: Find ungapped seeds
    print("Step 1: Finding ungapped seed alignments...")
    finder = ShortAlignmentFinder(min_length=8, max_mismatches=0)
    seeds = finder.find_alignments_between(seq1, seq2)

    print(f"Found {len(seeds)} ungapped seeds")
    for i, seed in enumerate(seeds, 1):
        print(f"  Seed {i}: pos1={seed['pos1']}, pos2={seed['pos2']}, "
              f"length={seed['length']}, seq='{seed['seq1']}'")
    print()

    # Step 2: Extend with gaps
    if seeds:
        print("Step 2: Extending best seed with gaps allowed...")
        gapped = finder.extend_alignment_with_gaps(
            seq1, seq2, seeds[0],
            min_identity=0.75,  # 75% identity minimum
            max_extension=20
        )

        if gapped:
            print(f"✓ Gapped extension successful!")
            print(f"  Position in seq1: {gapped['pos1']}-{gapped['end1']}")
            print(f"  Position in seq2: {gapped['pos2']}-{gapped['end2']}")
            print(f"  Alignment length: {gapped['alignment_length']}")
            print(f"  Matches: {gapped['matches']}")
            print(f"  Mismatches: {gapped['mismatches']}")
            print(f"  Gaps: {gapped['gaps']}")
            print(f"  Identity: {gapped['identity']:.1%} ({gapped['matches']}/{gapped['alignment_length']})")
            print()
            print("  Alignment visualization:")
            print(f"    Seq1: {gapped['seq1_aligned']}")
            print(f"          {gapped['alignment_string']}")
            print(f"    Seq2: {gapped['seq2_aligned']}")
        else:
            print("✗ Extension did not meet identity threshold")
    print()


def demo_real_cast_examples():
    """Demo 2: Real CAST homing spacer examples"""
    print("=" * 70)
    print("DEMO 2: Real CAST Homing Spacer Examples")
    print("=" * 70)

    # Real examples from CAST homing spacer and att site alignments
    examples = [
        {
            'name': 'Example 1',
            'seq1': "AAACGAACCACCUUCAUUCGAGUACAAGAGCA",  # 33bp (has U)
            'seq2': "AAGCGTACAACGTTTATCCGAGTTCAAGAGCA",  # 32bp
        },
        {
            'name': 'Example 2',
            'seq1': "AAAUUUACAUCUCCAGGCUCAGCCAAAAAGC",  # 32bp (has U)
            'seq2': "AAATTAACGTCGCCAGGTTCAGCTAAAAAAC",  # 31bp
        },
        {
            'name': 'Example 3',
            'seq1': "AGGACAGGAAGAAAACACCCAAGUUGGG",  # 28bp (has U)
            'seq2': "AGGACCGGAAGGTAGCAGCCAAGGCGGG",  # 28bp
        },
        {
            'name': 'Example 4',
            'seq1': "AAAAACAAUCGAUUCGAUGUUUAAGUUAA",  # 29bp (has U)
            'seq2': "AAAAAGAACCGTTTTGATGTCTATGTAAA",  # 29bp
        },
    ]

    for example in examples:
        print(f"\n{example['name']}:")
        # Convert U to T for DNA analysis
        seq1 = example['seq1'].replace('U', 'T')
        seq2 = example['seq2'].replace('U', 'T')

        print(f"  Spacer: {example['seq1']} ({len(seq1)}bp)")
        print(f"  AttSite: {example['seq2']} ({len(seq2)}bp)")

        # Find seeds with some mismatches allowed
        seeds = find_alignments_between_sequences(seq1, seq2, min_length=8, max_mismatches=2)

        if seeds:
            # Try different extension strategies
            best_result = None
            best_identity = 0

            for seed in seeds[:3]:  # Try top 3 seeds
                gapped = extend_alignment_with_gaps(
                    seq1, seq2, seed,
                    min_identity=0.70,  # Lower threshold to see results
                    max_extension=40
                )
                if gapped and gapped['identity'] > best_identity:
                    best_result = gapped
                    best_identity = gapped['identity']

            if best_result:
                ratio = f"{best_result['matches']}/{best_result['alignment_length']}"
                print(f"  ✓ Best alignment: {ratio} ({best_result['identity']:.1%})")
                print(f"    Matches:    {best_result['matches']}")
                print(f"    Mismatches: {best_result['mismatches']}")
                print(f"    Gaps:       {best_result['gaps']}")
                print()
                print("    Alignment:")
                print(f"      {best_result['seq1_aligned']}")
                print(f"      {best_result['alignment_string']}")
                print(f"      {best_result['seq2_aligned']}")
            else:
                print(f"  ✗ No extension met threshold")
        else:
            print(f"  ✗ No seeds found")
    print()


def demo_alignment_analysis():
    """Demo 3: Detailed alignment analysis"""
    print("=" * 70)
    print("DEMO 3: Detailed Alignment Analysis")
    print("=" * 70)

    # Use first real example for detailed analysis
    seq1 = "AAACGAACCACCTTCATTCGAGTACAAGAGCA"  # U->T conversion
    seq2 = "AAGCGTACAACGTTTATCCGAGTTCAAGAGCA"

    print(f"Analyzing alignment between:")
    print(f"  Seq1: {seq1} ({len(seq1)}bp)")
    print(f"  Seq2: {seq2} ({len(seq2)}bp)")
    print()

    # Find all possible seeds with different parameters
    print("Searching for seeds with different mismatch tolerances:")
    for max_mm in [0, 1, 2]:
        seeds = find_alignments_between_sequences(
            seq1, seq2,
            min_length=8,
            max_mismatches=max_mm
        )
        print(f"  {max_mm} mismatches: {len(seeds)} seeds found")

    print()

    # Use seeds with 2 mismatches for better coverage
    seeds = find_alignments_between_sequences(seq1, seq2, min_length=8, max_mismatches=2)

    if seeds:
        print(f"Extending top {min(3, len(seeds))} seeds:")
        for i, seed in enumerate(seeds[:3], 1):
            print(f"\n  Seed {i}: pos1={seed['pos1']}, pos2={seed['pos2']}, "
                  f"length={seed['length']}, mismatches={seed['mismatches']}")

            gapped = extend_alignment_with_gaps(
                seq1, seq2, seed,
                min_identity=0.70,
                max_extension=40
            )

            if gapped:
                ratio = f"{gapped['matches']}/{gapped['alignment_length']}"
                print(f"    Extended: {ratio} ({gapped['identity']:.1%}), "
                      f"{gapped['gaps']} gaps, {gapped['mismatches']} mismatches")
            else:
                print(f"    Could not extend (below threshold)")
    print()


def demo_reverse_complement():
    """Demo 4: Gapped extension with reverse complement"""
    print("=" * 70)
    print("DEMO 4: Gapped Extension with Reverse Complement")
    print("=" * 70)

    seq1 = "ACGTACGTACGTACGT"
    # Reverse complement of seq1 with 2bp gap
    seq2 = "ACGTACGNNACGTACGT"

    from utils.parsers import reverse_complement
    seq2_rc = reverse_complement(seq1)
    seq2_with_gap = seq2_rc[:7] + "NN" + seq2_rc[9:]  # Insert 2bp gap

    print(f"Sequence 1:    {seq1}")
    print(f"Sequence 2 RC: {seq2_with_gap}")
    print()

    print("Finding reverse complement alignments...")
    finder = ShortAlignmentFinder(min_length=7, max_mismatches=0,
                                  check_forward=False, check_revcomp=True)
    seeds = finder.find_alignments_between(seq1, seq2_with_gap)

    print(f"Found {len(seeds)} RC seeds")

    if seeds:
        gapped = finder.extend_alignment_with_gaps(
            seq1, seq2_with_gap, seeds[0],
            min_identity=0.80,
            max_extension=20
        )

        if gapped:
            print(f"✓ Extended RC alignment: {gapped['matches']}/{gapped['alignment_length']} ({gapped['identity']:.1%})")
            print()
            print("  Alignment:")
            print(f"    {gapped['seq1_aligned']}")
            print(f"    {gapped['alignment_string']}")
            print(f"    {gapped['seq2_aligned']}")
        else:
            print("✗ Did not meet identity threshold")
    print()


def demo_parameter_tuning():
    """Demo 5: Effect of different parameters"""
    print("=" * 70)
    print("DEMO 5: Parameter Tuning")
    print("=" * 70)

    seq1 = "ACGTACGTACGTACGTACGT"
    seq2 = "ACGTACGNNACGTACGTACGT"

    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print()

    seeds = find_alignments_between_sequences(seq1, seq2, min_length=8)

    if not seeds:
        print("No seeds found")
        return

    # Test different identity thresholds
    print("Testing different identity thresholds:")
    for threshold in [0.70, 0.75, 0.80, 0.85, 0.90]:
        gapped = extend_alignment_with_gaps(
            seq1, seq2, seeds[0],
            min_identity=threshold,
            max_extension=25
        )

        if gapped:
            print(f"  {threshold:.0%}: ✓ {gapped['matches']}/{gapped['alignment_length']} "
                  f"({gapped['identity']:.1%})")
        else:
            print(f"  {threshold:.0%}: ✗ Below threshold")

    print()

    # Test different gap penalties
    print("Testing different gap penalties (80% identity):")
    for gap_penalty in [-5, -3, -1]:
        gapped = extend_alignment_with_gaps(
            seq1, seq2, seeds[0],
            min_identity=0.80,
            gap_open_penalty=gap_penalty,
            gap_extend_penalty=-1,
            max_extension=25
        )

        if gapped:
            print(f"  Gap open={gap_penalty}: "
                  f"{gapped['matches']}/{gapped['alignment_length']} "
                  f"({gapped['gaps']} gaps, {gapped['identity']:.1%})")
        else:
            print(f"  Gap open={gap_penalty}: Below threshold")

    print()


def main():
    """Run all demos"""
    print("\n")
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 15 + "GAPPED ALIGNMENT EXTENSION DEMO" + " " * 22 + "║")
    print("║" + " " * 12 + "Finding CAST Homing Spacer Patterns" + " " * 21 + "║")
    print("╚" + "=" * 68 + "╝")
    print()

    demo_basic_gapped_extension()
    demo_real_cast_examples()
    demo_alignment_analysis()
    demo_reverse_complement()
    demo_parameter_tuning()

    print("=" * 70)
    print("All demos completed!")
    print("=" * 70)
    print()
    print("Key takeaways:")
    print("  • Start with ungapped seeds for fast initial detection")
    print("  • Extend seeds with gaps to find longer patterns")
    print("  • Real CAST examples show 70-85% identity is typical")
    print("  • Allow 2-3 mismatches in seeds for better coverage")
    print("  • Adjust identity threshold based on your requirements:")
    print("    - 85% for stringent matches")
    print("    - 70-80% for more permissive matches")
    print("  • Tune gap penalties to control gap vs. mismatch preference")
    print()


if __name__ == "__main__":
    main()
