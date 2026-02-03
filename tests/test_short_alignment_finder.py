"""
Unit tests for the Short Alignment Finder module.

Tests cover:
- Direct repeat finding (forward matches)
- Inverted repeat finding (reverse complement matches)
- Auto-extension of matches
- Mismatch tolerance
- Distance constraints
- Edge cases

Author: Kuang Hu
Date: 2026-01-26
"""

import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules.short_alignment_finder import (
    ShortAlignmentFinder,
    find_short_alignments,
    find_alignments_between_sequences
)


class TestBasicFunctionality:
    """Test basic repeat finding functionality."""

    def test_perfect_direct_repeat(self):
        """Test finding a perfect direct repeat."""
        # Sequence with a 10bp perfect direct repeat
        sequence = "ACGTACGTAC" + "NNNN" + "ACGTACGTAC"
        #           0-9              14-23

        finder = ShortAlignmentFinder(min_length=9, max_mismatches=0)
        results = finder.find_alignments(sequence)

        assert len(results) == 1
        assert results[0]['pos1'] == 0
        assert results[0]['pos2'] == 14
        assert results[0]['length'] == 10
        assert results[0]['distance'] == 14
        assert results[0]['seq1'] == "ACGTACGTAC"
        assert results[0]['seq2'] == "ACGTACGTAC"
        assert results[0]['mismatches'] == 0
        assert results[0]['mismatch_positions'] == []
        assert results[0]['orientation'] == 'forward'

    def test_inverted_repeat(self):
        """Test finding an inverted repeat (reverse complement)."""
        # ACGTACGTAC has reverse complement GTACGTACGT
        # Use GGGG as separator instead of NNNN to avoid extension issues
        sequence = "ACGTACGTAC" + "GGGG" + "GTACGTACGT"
        #           0-9              14-23

        finder = ShortAlignmentFinder(min_length=9, max_mismatches=0,
                                     check_forward=False, check_revcomp=True)
        results = finder.find_alignments(sequence)

        assert len(results) >= 1
        # Find the best reverse complement match
        rc_matches = [r for r in results if r['orientation'] == 'reverse_complement']
        assert len(rc_matches) >= 1
        best_match = max(rc_matches, key=lambda x: x['length'])
        assert best_match['length'] >= 9
        assert best_match['seq1'] in sequence
        assert best_match['orientation'] == 'reverse_complement'

    def test_both_orientations(self):
        """Test finding both direct and inverted repeats."""
        # Has both a direct repeat and an inverted repeat
        sequence = "AAAAAAAAAA" + "NN" + "AAAAAAAAAA" + "NN" + "TTTTTTTTTT"
        #           0-9              12-21            24-33

        finder = ShortAlignmentFinder(min_length=9, max_mismatches=0)
        results = finder.find_alignments(sequence)

        # Should find both the direct repeat and inverted repeat
        assert len(results) >= 2

        # Check we have both orientations
        orientations = [r['orientation'] for r in results]
        assert 'forward' in orientations
        assert 'reverse_complement' in orientations


class TestAutoExtension:
    """Test the auto-extension feature."""

    def test_extension_consolidation(self):
        """Test that overlapping matches are consolidated to the longest one."""
        # A 12bp repeat: ACGTACGTACGT
        # The algorithm should find many 9bp matches that overlap,
        # but should consolidate them into one 12bp match
        # Note: ACGTACGTACGT is self-reverse-complementary, so disable revcomp check
        sequence = "ACGTACGTACGT" + "NNNN" + "ACGTACGTACGT"

        finder = ShortAlignmentFinder(min_length=9, max_mismatches=0,
                                     check_forward=True, check_revcomp=False)
        results = finder.find_alignments(sequence)

        # Should get exactly 1 result (the 12bp match), not multiple 9bp matches
        assert len(results) == 1
        assert results[0]['length'] == 12
        assert results[0]['seq1'] == "ACGTACGTACGT"

    def test_extension_with_mismatches(self):
        """Test extension continues even when encountering mismatches."""
        # 10bp repeat with 1 mismatch in the middle
        sequence = "ACGTACGTAC" + "NNNN" + "ACGTGCGTAC"
        #           ACGTACGTAC              ACGTGCGTAC
        #               ^                       ^ (mismatch at position 4)

        finder = ShortAlignmentFinder(min_length=8, max_mismatches=1)
        results = finder.find_alignments(sequence)

        assert len(results) >= 1
        best = max(results, key=lambda x: x['length'])
        assert best['length'] == 10
        assert best['mismatches'] == 1
        assert 4 in best['mismatch_positions']


class TestMismatchHandling:
    """Test mismatch tolerance."""

    def test_no_mismatches_allowed(self):
        """Test that matches with mismatches are rejected when max_mismatches=0."""
        sequence = "ACGTACGTAC" + "NNNN" + "ACGTGCGTAC"  # 1 mismatch

        finder = ShortAlignmentFinder(min_length=9, max_mismatches=0)
        results = finder.find_alignments(sequence)

        # Should find shorter perfect matches, but not the full 10bp with mismatch
        for result in results:
            assert result['mismatches'] == 0

    def test_one_mismatch_allowed(self):
        """Test finding matches with 1 mismatch."""
        sequence = "ACGTACGTAC" + "NNNN" + "ACGTGCGTAC"  # 1 mismatch at pos 4

        finder = ShortAlignmentFinder(min_length=9, max_mismatches=1)
        results = finder.find_alignments(sequence)

        # Should find the 10bp match with 1 mismatch
        assert len(results) >= 1
        match_found = any(r['length'] == 10 and r['mismatches'] == 1 for r in results)
        assert match_found

    def test_mismatch_positions_recorded(self):
        """Test that mismatch positions are correctly recorded."""
        # Two mismatches at positions 2 and 6 (0-indexed)
        # ACGTACGTAC vs ACATACATAC:
        #   pos 0: A=A, pos 1: C=C, pos 2: G!=A, pos 3: T=T,
        #   pos 4: A=A, pos 5: C=C, pos 6: G!=A, pos 7: T=T,
        #   pos 8: A=A, pos 9: C=C
        sequence = "ACGTACGTAC" + "NNNN" + "ACATACATAC"

        finder = ShortAlignmentFinder(min_length=9, max_mismatches=2)
        results = finder.find_alignments(sequence)

        match_with_2_mm = [r for r in results if r['mismatches'] == 2]
        assert len(match_with_2_mm) >= 1
        assert 2 in match_with_2_mm[0]['mismatch_positions']
        assert 6 in match_with_2_mm[0]['mismatch_positions']


class TestDistanceConstraints:
    """Test min and max distance constraints."""

    def test_min_distance(self):
        """Test minimum distance constraint."""
        # Two repeats at different distances
        sequence = "AAAAAAAAAA" + "NN" + "AAAAAAAAAA" + "NNNNNNNN" + "AAAAAAAAAA"
        #           0-9              12-21                          30-39

        # Only allow matches with distance >= 20
        finder = ShortAlignmentFinder(min_length=9, min_distance=20)
        results = finder.find_alignments(sequence)

        # Should only find the match between pos1=0 and pos2=30 (distance=30)
        # NOT the match between pos1=0 and pos2=12 (distance=12)
        for result in results:
            assert result['distance'] >= 20

    def test_max_distance(self):
        """Test maximum distance constraint."""
        sequence = "AAAAAAAAAA" + "NN" + "AAAAAAAAAA" + "NNNNNNNN" + "AAAAAAAAAA"
        #           0-9              12-21                          30-39

        # Only allow matches with distance <= 15
        finder = ShortAlignmentFinder(min_length=9, max_distance=15)
        results = finder.find_alignments(sequence)

        # Should only find the match between pos1=0 and pos2=12 (distance=12)
        # NOT the match between pos1=0 and pos2=30 (distance=30)
        for result in results:
            assert result['distance'] <= 15

    def test_distance_range(self):
        """Test both min and max distance together."""
        sequence = "A" * 10 + "N" * 5 + "A" * 10 + "N" * 20 + "A" * 10

        # Only matches with distance between 12 and 25
        finder = ShortAlignmentFinder(min_length=9, min_distance=12, max_distance=25)
        results = finder.find_alignments(sequence)

        for result in results:
            assert 12 <= result['distance'] <= 25


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_sequence_too_short(self):
        """Test with sequence shorter than min_length."""
        sequence = "ACGT"  # Only 4bp

        finder = ShortAlignmentFinder(min_length=9)
        results = finder.find_alignments(sequence)

        assert len(results) == 0

    def test_no_repeats(self):
        """Test with sequence that has no repeats."""
        sequence = "ACGTACGTACGTACGT" * 2  # Continuous pattern, no true repeats

        finder = ShortAlignmentFinder(min_length=30)  # Looking for 30bp repeats
        results = finder.find_alignments(sequence)

        # May find some repeats due to the pattern, but test should run without error
        assert isinstance(results, list)

    def test_lowercase_sequence(self):
        """Test that lowercase sequences are handled correctly."""
        sequence = "acgtacgtac" + "nnnn" + "acgtacgtac"

        finder = ShortAlignmentFinder(min_length=9)
        results = finder.find_alignments(sequence)

        assert len(results) == 1
        assert results[0]['seq1'] == "ACGTACGTAC"  # Should be uppercase

    def test_convenience_function(self):
        """Test the convenience function works correctly."""
        sequence = "ACGTACGTAC" + "NNNN" + "ACGTACGTAC"

        results = find_short_alignments(sequence, min_length=9)

        assert len(results) == 1
        assert results[0]['length'] == 10


class TestParameterValidation:
    """Test parameter validation."""

    def test_invalid_min_length(self):
        """Test that invalid min_length raises error."""
        with pytest.raises(ValueError):
            ShortAlignmentFinder(min_length=0)

    def test_negative_mismatches(self):
        """Test that negative max_mismatches raises error."""
        with pytest.raises(ValueError):
            ShortAlignmentFinder(max_mismatches=-1)

    def test_invalid_distance_range(self):
        """Test that max_distance < min_distance raises error."""
        with pytest.raises(ValueError):
            ShortAlignmentFinder(min_distance=100, max_distance=50)

    def test_no_orientation_selected(self):
        """Test that at least one orientation must be checked."""
        with pytest.raises(ValueError):
            ShortAlignmentFinder(check_forward=False, check_revcomp=False)


class TestRealWorldScenarios:
    """Test with real-world-like scenarios."""

    def test_is110_like_sequence(self):
        """Test with IS110-like parameters and sequence."""
        # Simulate a transposon upstream region with repeats
        sequence = (
            "ATGCATGCATGCATGC"  # Some sequence
            "TTTAAACCCGGGAAA"   # 15bp repeat 1
            "NNNNNNNNNNNNNNNNNNNNN"  # Gap
            "TTTAAACCCGGGAAA"   # 15bp repeat 2 (perfect match)
            "GCTAGCTAGCTAGCTA"  # More sequence
        )

        # IS110 parameters: min 9bp, no mismatches
        finder = ShortAlignmentFinder(min_length=9, max_mismatches=0)
        results = finder.find_alignments(sequence)

        # Should find the 15bp repeat
        long_matches = [r for r in results if r['length'] >= 15]
        assert len(long_matches) >= 1

    def test_multiple_overlapping_repeats(self):
        """Test with multiple overlapping repeats."""
        # A sequence with nested repeats
        sequence = "ATATATATAT" + "NN" + "ATATATATAT" + "NN" + "ATATATATAT"

        finder = ShortAlignmentFinder(min_length=8, max_mismatches=0)
        results = finder.find_alignments(sequence)

        # Should find multiple matches between different positions
        assert len(results) >= 3  # At least 3 pairs of repeats

        # All should be valid (within distance constraints, etc.)
        for result in results:
            assert result['length'] >= 8
            assert result['mismatches'] == 0


class TestBetweenSequences:
    """Test finding alignments between two different sequences."""

    def test_basic_between_sequences(self):
        """Test finding alignments between two sequences."""
        # Flanking region with a 10bp motif
        seq1 = "ATGCATGC" + "ACGTACGTAC" + "NNNN"
        # Non-coding region with the same 10bp motif
        seq2 = "NNNN" + "ACGTACGTAC" + "GCTAGCTA"

        finder = ShortAlignmentFinder(min_length=9)
        results = finder.find_alignments_between(seq1, seq2)

        assert len(results) >= 1
        # Find the 10bp match
        match = [r for r in results if r['length'] == 10][0]
        assert match['pos1'] == 8  # Position in seq1
        assert match['pos2'] == 4  # Position in seq2
        assert match['seq1'] == "ACGTACGTAC"
        assert match['seq2'] == "ACGTACGTAC"
        assert match['mismatches'] == 0
        assert match['orientation'] == 'forward'

    def test_between_sequences_reverse_complement(self):
        """Test finding reverse complement matches between sequences."""
        # Sequence 1 has forward motif
        seq1 = "ATGCATGC" + "ACGTACGTAC" + "GGGG"
        # Sequence 2 has reverse complement: RC(ACGTACGTAC) = GTACGTACGT
        seq2 = "GGGG" + "GTACGTACGT" + "GCTAGCTA"

        finder = ShortAlignmentFinder(min_length=9, check_forward=False,
                                     check_revcomp=True)
        results = finder.find_alignments_between(seq1, seq2)

        assert len(results) >= 1
        # Find reverse complement matches
        rc_matches = [r for r in results if r['orientation'] == 'reverse_complement']
        assert len(rc_matches) >= 1
        best_match = max(rc_matches, key=lambda x: x['length'])
        assert best_match['length'] >= 9
        assert best_match['orientation'] == 'reverse_complement'

    def test_between_sequences_with_mismatches(self):
        """Test between-sequences with mismatches allowed."""
        seq1 = "ATGC" + "ACGTACGTAC" + "NNNN"
        seq2 = "NNNN" + "ACGTGCGTAC" + "GCTA"  # 1 mismatch at position 4

        # With no mismatches allowed
        finder_strict = ShortAlignmentFinder(min_length=9, max_mismatches=0)
        results_strict = finder_strict.find_alignments_between(seq1, seq2)
        # Should find shorter perfect matches
        for r in results_strict:
            assert r['mismatches'] == 0

        # With 1 mismatch allowed
        finder_relaxed = ShortAlignmentFinder(min_length=9, max_mismatches=1)
        results_relaxed = finder_relaxed.find_alignments_between(seq1, seq2)
        # Should find the 10bp match with 1 mismatch
        long_match = [r for r in results_relaxed if r['length'] == 10]
        assert len(long_match) >= 1
        assert long_match[0]['mismatches'] == 1
        assert 4 in long_match[0]['mismatch_positions']

    def test_between_sequences_no_distance_constraint(self):
        """Test that distance constraints don't apply between sequences."""
        # Two sequences with matches at any position
        seq1 = "AAAAAAAAAA" + "NNNN" + "CCCCCCCCCC"
        seq2 = "AAAAAAAAAA" + "NNNNNNNNNN" + "CCCCCCCCCC"

        # Even with distance constraints set, they shouldn't affect between-sequence search
        finder = ShortAlignmentFinder(min_length=9, min_distance=100, max_distance=200)
        results = finder.find_alignments_between(seq1, seq2)

        # Should find matches regardless of their positions
        assert len(results) >= 2  # Should find both the A and C matches

    def test_convenience_function_between_sequences(self):
        """Test the convenience function for between-sequences."""
        seq1 = "ATGCACGTACGTAC"
        seq2 = "NNNNACGTACGTAC"

        results = find_alignments_between_sequences(seq1, seq2, min_length=9)

        assert len(results) >= 1
        assert results[0]['length'] >= 9

    def test_between_sequences_transposon_use_case(self):
        """Test realistic transposon flanking vs non-coding scenario."""
        # Flanking region (e.g., 100bp upstream of transposon)
        flanking = (
            "ATGCTAGCTAGCTAGC"      # Random upstream sequence (16bp)
            "TTTAAACCCGGG"           # 12bp potential guide RNA sequence
            "ACTAGCTAGCTAGCTA"       # More upstream (16bp) - starts with A, not G
        )

        # Non-coding region of transposon - use different chars to prevent over-extension
        noncoding = (
            "CCCCCCCCCCCC"           # Some sequence (12bp)
            "TTTAAACCCGGG"           # Same 12bp guide RNA sequence
            "ACACACACACAC"           # More sequence (different from flanking)
            "GTACGTACGTAC"           # Different sequence
            "CCCCCCCC"
        )

        # Find guide RNA pairs between flanking and non-coding regions
        finder = ShortAlignmentFinder(min_length=9, max_mismatches=0,
                                     check_forward=True, check_revcomp=False)
        results = finder.find_alignments_between(flanking, noncoding)

        # Should find the 12bp guide RNA sequence
        guide_matches = [r for r in results if r['length'] >= 12]
        assert len(guide_matches) >= 1

        # Verify the match contains the guide sequence
        match = guide_matches[0]
        assert "TTTAAACCCGGG" in match['seq1']
        assert "TTTAAACCCGGG" in match['seq2']
        assert match['orientation'] == 'forward'

    def test_between_sequences_auto_extension(self):
        """Test that auto-extension works between sequences."""
        # 15bp match, should extend from 9bp seed
        seq1 = "ATGC" + "ACGTACGTACGTACG" + "NNNN"
        seq2 = "NNNN" + "ACGTACGTACGTACG" + "GCTA"

        finder = ShortAlignmentFinder(min_length=9)
        results = finder.find_alignments_between(seq1, seq2)

        # Should find the full 15bp match, not multiple 9bp matches
        max_length = max(r['length'] for r in results)
        assert max_length == 15

    def test_between_empty_sequences(self):
        """Test with one or both sequences too short."""
        seq1 = "ACGT"  # Only 4bp
        seq2 = "ACGTACGTACGT"  # 12bp

        finder = ShortAlignmentFinder(min_length=9)
        results = finder.find_alignments_between(seq1, seq2)

        # seq1 is too short for 9bp matches
        assert len(results) == 0

    def test_between_sequences_both_orientations(self):
        """Test finding both forward and reverse complement between sequences."""
        seq1 = "AAAAAAAAAA" + "NN" + "GGGGGGGGGG"
        # AAAAAAAAAA has RC = TTTTTTTTTT
        # GGGGGGGGGG has RC = CCCCCCCCCC
        seq2 = "AAAAAAAAAA" + "NN" + "TTTTTTTTTT" + "NN" + "CCCCCCCCCC"

        finder = ShortAlignmentFinder(min_length=9, check_forward=True,
                                     check_revcomp=True)
        results = finder.find_alignments_between(seq1, seq2)

        # Should find both forward match (A-A) and RC matches (A-T, G-C)
        orientations = [r['orientation'] for r in results]
        assert 'forward' in orientations
        assert 'reverse_complement' in orientations


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
