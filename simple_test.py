#!/usr/bin/env python3
"""
Simple Two-Step Test: Seed Finding → Gapped Extension

Step 1: Find short ungapped seeds (≥9bp, ≤1 mismatch)
Step 2: Extend seeds with gaps to longer alignments
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from modules.short_alignment_finder import (
    find_alignments_between_sequences,
    extend_alignment_with_gaps
)


def test_two_step(name, seq1, seq2):
    """Simple two-step test"""
    print(f"\n{'=' * 70}")
    print(f"{name}")
    print(f"{'=' * 70}")

    seq1 = seq1.replace('U', 'T')
    seq2 = seq2.replace('U', 'T')

    print(f"Spacer:  {seq1} ({len(seq1)}bp)")
    print(f"AttSite: {seq2} ({len(seq2)}bp)\n")

    # STEP 1: Find seeds
    print("STEP 1: Find short ungapped seeds (≥9bp, ≤1 mismatch)")
    seeds = find_alignments_between_sequences(
        seq1, seq2,
        min_length=9,
        max_mismatches=1
    )

    if not seeds:
        print("  ❌ No seeds found\n")
        return False

    print(f"  ✓ Found {len(seeds)} seed(s):")
    for i, seed in enumerate(seeds, 1):
        pattern = f"{seed['length']-seed['mismatches']}/{seed['length']}"
        print(f"    Seed {i}: {pattern} at pos1={seed['pos1']}, pos2={seed['pos2']}")

    # STEP 2: Extend with gaps
    print(f"\nSTEP 2: Extend seed to longer alignment (with gaps allowed)")

    result = extend_alignment_with_gaps(
        seq1, seq2, seeds[0],
        min_identity=0.75,
        max_extension=50
    )

    if not result:
        print("  ❌ Extension failed\n")
        return False

    print(f"  ✓ Extended to {result['alignment_length']}bp:")
    print(f"\n  {result['seq1_aligned']}")
    print(f"  {result['alignment_string']}")
    print(f"  {result['seq2_aligned']}\n")
    print(f"  Pattern:  {result['matches']}/{result['alignment_length']} ({result['identity']:.1%})")
    print(f"  Length:   {result['alignment_length']}bp (≥25bp: {'✅' if result['alignment_length']>=25 else '❌'})")
    print(f"  Identity: {result['identity']:.1%} (≥75%: {'✅' if result['identity']>=0.75 else '❌'})\n")

    return True


def main():
    print("\n" + "=" * 70)
    print("TWO-STEP PROCESS: Seed Finding → Gapped Extension")
    print("=" * 70)

    examples = [
        ("Example 1", "AAACGAACCACCUUCAUUCGAGUACAAGAGCA", "AAGCGTACAACGTTTATCCGAGTTCAAGAGCA"),
        ("Example 2", "AAAUUUACAUCUCCAGGCUCAGCCAAAAAGC", "AAATTAACGTCGCCAGGTTCAGCTAAAAAAC"),
        ("Example 3", "AGGACAGGAAGAAAACACCCAAGUUGGG", "AGGACCGGAAGGTAGCAGCCAAGGCGGG"),
        ("Example 4", "AAAAACAAUCGAUUCGAUGUUUAAGUUAA", "AAAAAGAACCGTTTTGATGTCTATGTAAA"),
    ]

    passed = 0
    for name, seq1, seq2 in examples:
        if test_two_step(name, seq1, seq2):
            passed += 1

    print("=" * 70)
    print(f"RESULT: {passed}/{len(examples)} examples successfully extended")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
