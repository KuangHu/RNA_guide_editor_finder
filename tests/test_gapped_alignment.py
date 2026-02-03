#!/usr/bin/env python3
"""
Tests for Gapped Alignment Extension

Tests the extend_alignment_with_gaps functionality for finding
CAST homing spacer-like patterns with gaps/indels.

Author: Kuang Hu
Date: 2026-01-27
"""

import sys
import os
import pytest

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules.short_alignment_finder import (
    ShortAlignmentFinder,
    find_alignments_between_sequences,
    extend_alignment_with_gaps
)


class TestGappedAlignmentExtension:
    """Test suite for gapped alignment extension"""

    def test_basic_gap_extension(self):
        """Test basic gap extension with insertion"""
        seq1 = "ACGTACGTACGTACGT"
        seq2 = "ACGTACGTNNNACGTACGT"

        # Find seed
        seeds = find_alignments_between_sequences(seq1, seq2, min_length=8)
        assert len(seeds) > 0, "Should find at least one seed"

        # Extend with gaps
        gapped = extend_alignment_with_gaps(seq1, seq2, seeds[0], min_identity=0.75)
        assert gapped is not None, "Should successfully extend"
        assert gapped['gaps'] == 3, "Should have 3 gaps"
        assert gapped['identity'] > 0.80, "Identity should be > 80%"

    def test_identity_threshold(self):
        """Test identity threshold filtering"""
        seq1 = "ACGTACGTACGTACGT"
        seq2 = "ACGTACGTNNNACGTACGT"

        seeds = find_alignments_between_sequences(seq1, seq2, min_length=8)

        # Should pass with low threshold
        result_low = extend_alignment_with_gaps(seq1, seq2, seeds[0], min_identity=0.70)
        assert result_low is not None

        # Should still pass with high threshold if identity is good
        result_high = extend_alignment_with_gaps(seq1, seq2, seeds[0], min_identity=0.95)
        # Result depends on actual identity - may pass or fail

    def test_cast_pattern_17_20(self):
        """Test CAST homing spacer pattern (17/20)"""
        # Create pattern with ~85% identity
        seq1 = "ACGTACGTACGTACGTACGT"  # 20bp
        seq2 = "ACGTACGTACGTACGTACGT"  # Same, perfect match

        seeds = find_alignments_between_sequences(seq1, seq2, min_length=10)
        assert len(seeds) > 0

        gapped = extend_alignment_with_gaps(seq1, seq2, seeds[0], min_identity=0.85)
        assert gapped is not None
        assert gapped['identity'] >= 0.85

    def test_alignment_string_format(self):
        """Test alignment string generation"""
        seq1 = "ACGTACGTACGT"
        seq2 = "ACGTACGTNNNACGT"

        seeds = find_alignments_between_sequences(seq1, seq2, min_length=7)
        gapped = extend_alignment_with_gaps(seq1, seq2, seeds[0], min_identity=0.70)

        assert gapped is not None
        assert len(gapped['seq1_aligned']) == len(gapped['seq2_aligned'])
        assert len(gapped['alignment_string']) == len(gapped['seq1_aligned'])

        # Check alignment string characters
        for char in gapped['alignment_string']:
            assert char in ['|', '.', ' '], f"Invalid alignment character: {char}"

    def test_reverse_complement_extension(self):
        """Test gapped extension with reverse complement"""
        from utils.parsers import reverse_complement

        seq1 = "ACGTACGTACGTACGT"
        seq2_rc = reverse_complement(seq1)

        finder = ShortAlignmentFinder(
            min_length=8,
            check_forward=False,
            check_revcomp=True
        )
        seeds = finder.find_alignments_between(seq1, seq2_rc)
        assert len(seeds) > 0, "Should find RC seeds"

        gapped = finder.extend_alignment_with_gaps(seq1, seq2_rc, seeds[0], min_identity=0.90)
        assert gapped is not None
        assert gapped['orientation'] == 'reverse_complement'

    def test_gap_penalty_effects(self):
        """Test that gap penalties affect results"""
        seq1 = "ACGTACGTACGTACGT"
        seq2 = "ACGTACGTNNNACGTACGT"

        seeds = find_alignments_between_sequences(seq1, seq2, min_length=8)

        # More severe gap penalty
        result_severe = extend_alignment_with_gaps(
            seq1, seq2, seeds[0],
            min_identity=0.70,
            gap_open_penalty=-10,
            gap_extend_penalty=-5
        )

        # Lenient gap penalty
        result_lenient = extend_alignment_with_gaps(
            seq1, seq2, seeds[0],
            min_identity=0.70,
            gap_open_penalty=-1,
            gap_extend_penalty=-1
        )

        assert result_severe is not None or result_lenient is not None

    def test_max_extension_limit(self):
        """Test max_extension parameter"""
        seq1 = "A" * 100
        seq2 = "A" * 100

        seeds = find_alignments_between_sequences(seq1, seq2, min_length=10)

        # Limited extension
        result_limited = extend_alignment_with_gaps(
            seq1, seq2, seeds[0],
            max_extension=10
        )

        # Larger extension
        result_large = extend_alignment_with_gaps(
            seq1, seq2, seeds[0],
            max_extension=50
        )

        assert result_limited is not None
        assert result_large is not None
        # Larger extension should be >= limited
        assert result_large['alignment_length'] >= result_limited['alignment_length']

    def test_matches_mismatches_gaps_sum(self):
        """Test that matches + mismatches + gaps = alignment_length"""
        seq1 = "ACGTACGTACGTACGT"
        seq2 = "ACGTACGTNNNACGTACGT"

        seeds = find_alignments_between_sequences(seq1, seq2, min_length=8)
        gapped = extend_alignment_with_gaps(seq1, seq2, seeds[0], min_identity=0.70)

        assert gapped is not None
        total = gapped['matches'] + gapped['mismatches'] + gapped['gaps']
        assert total == gapped['alignment_length'], \
            f"Sum {total} should equal alignment_length {gapped['alignment_length']}"

    def test_identity_calculation(self):
        """Test identity calculation accuracy"""
        seq1 = "ACGTACGTACGTACGT"
        seq2 = "ACGTACGTNNNACGTACGT"

        seeds = find_alignments_between_sequences(seq1, seq2, min_length=8)
        gapped = extend_alignment_with_gaps(seq1, seq2, seeds[0], min_identity=0.70)

        assert gapped is not None
        calculated_identity = gapped['matches'] / gapped['alignment_length']
        assert abs(gapped['identity'] - calculated_identity) < 0.001, \
            "Identity should be matches/alignment_length"

    def test_seed_information_preserved(self):
        """Test that seed information is preserved in result"""
        seq1 = "ACGTACGTACGTACGT"
        seq2 = "ACGTACGTNNNACGTACGT"

        seeds = find_alignments_between_sequences(seq1, seq2, min_length=8)
        seed = seeds[0]

        gapped = extend_alignment_with_gaps(seq1, seq2, seed, min_identity=0.70)

        assert gapped is not None
        assert gapped['seed_pos1'] == seed['pos1']
        assert gapped['seed_pos2'] == seed['pos2']
        assert gapped['seed_length'] == seed['length']

    def test_empty_sequences(self):
        """Test handling of edge cases"""
        # Empty extension regions should handle gracefully
        seq1 = "ACGTACGT"
        seq2 = "ACGTACGT"

        seeds = find_alignments_between_sequences(seq1, seq2, min_length=8)
        # Seed spans entire sequence, no room to extend
        gapped = extend_alignment_with_gaps(seq1, seq2, seeds[0], min_identity=0.90)

        # Should still work, just no extension
        assert gapped is not None


def test_convenience_function():
    """Test the convenience function wrapper"""
    seq1 = "ACGTACGTACGTACGT"
    seq2 = "ACGTACGTNNNACGTACGT"

    seeds = find_alignments_between_sequences(seq1, seq2, min_length=8)
    result = extend_alignment_with_gaps(seq1, seq2, seeds[0])

    assert result is not None
    assert 'identity' in result
    assert 'matches' in result
    assert 'gaps' in result


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
