"""
Short Alignment Finder Module

This module finds short DNA alignments (direct repeats and inverted repeats)
within a single DNA sequence. Useful for identifying repeat elements,
palindromic sequences, and potential regulatory motifs in transposons.

Author: Kuang Hu
Date: 2026-01-26
"""

from typing import List, Dict, Optional, Tuple
import sys
import os

# Add parent directory to path to import utils
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils.parsers import reverse_complement


class ShortAlignmentFinder:
    """
    Find short DNA alignments (repeats) in a DNA sequence.

    This class identifies both direct repeats (forward-forward matches) and
    inverted repeats (forward-reverse complement matches) with configurable
    parameters for length, mismatches, and distance constraints.

    Features:
    - Configurable minimum alignment length
    - Allow 0 or more mismatches
    - Distance constraints (min/max) between alignments
    - Auto-extension: extends matches to their maximum length
    - De-duplication: reports longest match, not overlapping shorter ones
    - Tracks mismatch positions for quality assessment

    Attributes:
        min_length (int): Minimum alignment length to report
        max_mismatches (int): Maximum number of mismatches allowed
        min_distance (int): Minimum distance between alignment positions (pos2 - pos1)
        max_distance (int or None): Maximum distance between positions (None = no limit)
        check_forward (bool): Whether to check for direct repeats
        check_revcomp (bool): Whether to check for inverted repeats
    """

    def __init__(self,
                 min_length: int = 9,
                 max_mismatches: int = 0,
                 min_distance: int = 0,
                 max_distance: Optional[int] = None,
                 check_forward: bool = True,
                 check_revcomp: bool = True):
        """
        Initialize the ShortAlignmentFinder.

        Parameters:
            min_length: Minimum alignment length to report (default: 9bp for IS110)
            max_mismatches: Maximum mismatches allowed (default: 0)
            min_distance: Minimum distance between positions, pos2-pos1 (default: 0)
            max_distance: Maximum distance between positions (default: None = no limit)
            check_forward: Find direct repeats (default: True)
            check_revcomp: Find inverted repeats (default: True)

        Raises:
            ValueError: If parameters are invalid
        """
        if min_length < 1:
            raise ValueError("min_length must be at least 1")
        if max_mismatches < 0:
            raise ValueError("max_mismatches cannot be negative")
        if min_distance < 0:
            raise ValueError("min_distance cannot be negative")
        if max_distance is not None and max_distance < min_distance:
            raise ValueError("max_distance must be >= min_distance")
        if not check_forward and not check_revcomp:
            raise ValueError("At least one of check_forward or check_revcomp must be True")

        self.min_length = min_length
        self.max_mismatches = max_mismatches
        self.min_distance = min_distance
        self.max_distance = max_distance
        self.check_forward = check_forward
        self.check_revcomp = check_revcomp

    def find_alignments(self, sequence: str) -> List[Dict]:
        """
        Find all short alignments in the sequence.

        Algorithm:
        1. Find all seed matches of min_length with <= max_mismatches
        2. Extend each seed match bidirectionally while staying within mismatch limit
        3. Consolidate overlapping matches, keeping the longest ones

        Parameters:
            sequence: DNA sequence to search (string, A/T/G/C)

        Returns:
            List of alignment dictionaries, each containing:
                - pos1: Start position of first occurrence (0-based)
                - pos2: Start position of second occurrence (0-based)
                - length: Length of the alignment
                - distance: Distance between positions (pos2 - pos1)
                - seq1: Sequence at pos1
                - seq2: Sequence at pos2 (as it appears in the input sequence)
                - mismatches: Number of mismatches
                - mismatch_positions: List of positions (relative to alignment) with mismatches
                - orientation: 'forward' or 'reverse_complement'

        Example:
            >>> finder = ShortAlignmentFinder(min_length=9, max_mismatches=1)
            >>> results = finder.find_alignments("ACGTACGTACGTNNNNACGTACGTACGT")
            >>> results[0]
            {
                'pos1': 0,
                'pos2': 16,
                'length': 12,
                'distance': 16,
                'seq1': 'ACGTACGTACGT',
                'seq2': 'ACGTACGTACGT',
                'mismatches': 0,
                'mismatch_positions': [],
                'orientation': 'forward'
            }
        """
        sequence = sequence.upper()

        # Find all seed matches
        seed_matches = self._find_seed_matches(sequence)

        # Extend each seed match to maximum length
        extended_matches = []
        for seed in seed_matches:
            extended = self._extend_match(sequence, seed)
            if extended:  # Only add if extension was successful
                extended_matches.append(extended)

        # Consolidate overlapping matches (keep longest)
        consolidated = self._consolidate_matches(extended_matches)

        # Sort by pos1, then pos2
        consolidated.sort(key=lambda x: (x['pos1'], x['pos2']))

        return consolidated

    def _find_seed_matches(self, sequence: str) -> List[Dict]:
        """
        Find all seed matches of min_length in the sequence.

        A seed match is a pair of positions (i, j) where:
        - The subsequences starting at i and j are at least min_length long
        - They match with <= max_mismatches
        - j - i is within [min_distance, max_distance]

        Returns:
            List of seed match dictionaries
        """
        seeds = []
        seq_len = len(sequence)

        # For each starting position
        for i in range(seq_len - self.min_length + 1):
            kmer = sequence[i:i + self.min_length]

            # Calculate search range for position j
            j_min = i + max(self.min_distance, self.min_length)  # Ensure no overlap
            j_max = seq_len - self.min_length + 1
            if self.max_distance is not None:
                j_max = min(j_max, i + self.max_distance + 1)

            # Search forward matches
            if self.check_forward:
                for j in range(j_min, j_max):
                    target = sequence[j:j + self.min_length]
                    mismatches, mismatch_pos = self._count_mismatches(kmer, target)

                    if mismatches <= self.max_mismatches:
                        seeds.append({
                            'pos1': i,
                            'pos2': j,
                            'length': self.min_length,
                            'orientation': 'forward',
                            'mismatches': mismatches,
                            'mismatch_positions': mismatch_pos
                        })

            # Search reverse complement matches
            if self.check_revcomp:
                kmer_rc = reverse_complement(kmer)
                for j in range(j_min, j_max):
                    target = sequence[j:j + self.min_length]
                    mismatches, mismatch_pos = self._count_mismatches(kmer_rc, target)

                    if mismatches <= self.max_mismatches:
                        seeds.append({
                            'pos1': i,
                            'pos2': j,
                            'length': self.min_length,
                            'orientation': 'reverse_complement',
                            'mismatches': mismatches,
                            'mismatch_positions': mismatch_pos
                        })

        return seeds

    def _extend_match(self, sequence: str, seed: Dict) -> Optional[Dict]:
        """
        Extend a seed match bidirectionally as far as possible.

        Extension continues while:
        - Total mismatches <= max_mismatches
        - We don't run off the sequence ends
        - For forward matches: we don't cause overlap (pos2 >= pos1 + length)

        Parameters:
            sequence: The full DNA sequence
            seed: Seed match dictionary

        Returns:
            Extended match dictionary, or None if extension fails
        """
        pos1 = seed['pos1']
        pos2 = seed['pos2']
        length = seed['length']
        orientation = seed['orientation']

        # Get the comparison sequence for pos2
        if orientation == 'forward':
            get_seq2 = lambda start, end: sequence[start:end]
        else:  # reverse_complement
            get_seq2 = lambda start, end: reverse_complement(sequence[start:end])

        # Extend right
        while True:
            # Check boundaries
            if pos1 + length >= len(sequence) or pos2 + length >= len(sequence):
                break

            # Get next bases
            base1 = sequence[pos1 + length]
            if orientation == 'forward':
                base2 = sequence[pos2 + length]
            else:
                # For reverse complement, we need the complement of the base before pos2
                # Actually, when extending right on pos1, we extend left on the RC of pos2
                # This is getting complex. Let me rethink this.
                #
                # When we have a reverse complement match:
                # pos1: ACGT...
                # pos2: ...ACGT (which is RC of what's at pos1)
                #
                # The sequence at pos2 is stored forward in the sequence string
                # So if we extend pos1 to the right by 1, we need to check if
                # the RC of that extended base matches the extension of pos2
                #
                # Actually, let me simplify: just compare the full extended sequences
                extended_seq1 = sequence[pos1:pos1 + length + 1]
                extended_seq2_region = sequence[pos2:pos2 + length + 1]

                if orientation == 'forward':
                    extended_seq2 = extended_seq2_region
                else:
                    extended_seq2 = reverse_complement(extended_seq2_region)

                mismatches, _ = self._count_mismatches(extended_seq1, extended_seq2)

                if mismatches <= self.max_mismatches:
                    length += 1
                else:
                    break
                continue

            # For forward orientation, simple comparison
            mismatch = (base1 != base2)

            # Check if adding this position would exceed mismatch limit
            test_seq1 = sequence[pos1:pos1 + length + 1]
            test_seq2 = sequence[pos2:pos2 + length + 1]
            test_mismatches, _ = self._count_mismatches(test_seq1, test_seq2)

            if test_mismatches <= self.max_mismatches:
                length += 1
            else:
                break

        # Extend left
        while True:
            # Check boundaries
            if pos1 == 0 or pos2 == 0:
                break

            # Check if adding this position would exceed mismatch limit
            test_seq1 = sequence[pos1 - 1:pos1 + length]
            if orientation == 'forward':
                test_seq2 = sequence[pos2 - 1:pos2 + length]
            else:
                test_seq2 = reverse_complement(sequence[pos2 - 1:pos2 + length])

            test_mismatches, _ = self._count_mismatches(test_seq1, test_seq2)

            if test_mismatches <= self.max_mismatches:
                pos1 -= 1
                pos2 -= 1
                length += 1
            else:
                break

        # Final check: ensure minimum length is met after extension
        if length < self.min_length:
            return None

        # Get final sequences and mismatch info
        seq1 = sequence[pos1:pos1 + length]
        seq2_raw = sequence[pos2:pos2 + length]

        if orientation == 'forward':
            seq2_compare = seq2_raw
        else:
            seq2_compare = reverse_complement(seq2_raw)

        mismatches, mismatch_positions = self._count_mismatches(seq1, seq2_compare)

        # Ensure we're within distance constraints after extension
        distance = pos2 - pos1
        if distance < self.min_distance:
            return None
        if self.max_distance is not None and distance > self.max_distance:
            return None

        return {
            'pos1': pos1,
            'pos2': pos2,
            'length': length,
            'distance': distance,
            'seq1': seq1,
            'seq2': seq2_raw,  # Store as it appears in sequence
            'mismatches': mismatches,
            'mismatch_positions': mismatch_positions,
            'orientation': orientation
        }

    def _consolidate_matches(self, matches: List[Dict]) -> List[Dict]:
        """
        Consolidate overlapping matches, keeping the longest ones.

        Two matches are considered overlapping if they have:
        - The same orientation
        - Overlapping regions at pos1 OR pos2

        When matches overlap, we keep the longest one. If lengths are equal,
        we keep the one with fewer mismatches.

        Parameters:
            matches: List of match dictionaries

        Returns:
            Consolidated list of matches
        """
        if not matches:
            return []

        # Sort by length (descending), then by mismatches (ascending)
        matches_sorted = sorted(matches,
                               key=lambda x: (-x['length'], x['mismatches']))

        consolidated = []

        for match in matches_sorted:
            # Check if this match overlaps with any already consolidated match
            is_redundant = False

            for existing in consolidated:
                if self._matches_overlap(match, existing):
                    is_redundant = True
                    break

            if not is_redundant:
                consolidated.append(match)

        return consolidated

    def _matches_overlap(self, match1: Dict, match2: Dict) -> bool:
        """
        Check if two matches overlap (indicating one is likely a subset of the other).

        Matches overlap if they have the same orientation and their regions overlap
        at both pos1 and pos2 positions.

        Parameters:
            match1, match2: Match dictionaries

        Returns:
            True if matches overlap
        """
        # Must be same orientation to overlap
        if match1['orientation'] != match2['orientation']:
            return False

        # Check if regions overlap at pos1
        pos1_overlap = self._regions_overlap(
            match1['pos1'], match1['pos1'] + match1['length'],
            match2['pos1'], match2['pos1'] + match2['length']
        )

        # Check if regions overlap at pos2
        pos2_overlap = self._regions_overlap(
            match1['pos2'], match1['pos2'] + match1['length'],
            match2['pos2'], match2['pos2'] + match2['length']
        )

        # Both must overlap for matches to be considered overlapping
        return pos1_overlap and pos2_overlap

    def _regions_overlap(self, start1: int, end1: int,
                        start2: int, end2: int) -> bool:
        """
        Check if two genomic regions overlap.

        Parameters:
            start1, end1: First region [start1, end1)
            start2, end2: Second region [start2, end2)

        Returns:
            True if regions overlap
        """
        return not (end1 <= start2 or end2 <= start1)

    def _count_mismatches(self, seq1: str, seq2: str) -> Tuple[int, List[int]]:
        """
        Count mismatches between two sequences of equal length.

        Parameters:
            seq1, seq2: DNA sequences to compare

        Returns:
            Tuple of (mismatch_count, list of mismatch positions)
        """
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be equal length")

        mismatches = 0
        mismatch_positions = []

        for i, (base1, base2) in enumerate(zip(seq1, seq2)):
            if base1 != base2:
                mismatches += 1
                mismatch_positions.append(i)

        return mismatches, mismatch_positions

    def find_alignments_between(self, sequence1: str, sequence2: str) -> List[Dict]:
        """
        Find all short alignments between two different sequences.

        This is useful for finding alignments between:
        - Flanking region vs non-coding region of transposon
        - Upstream region vs downstream region
        - Any two distinct genomic regions

        Algorithm:
        1. Find all seed matches between seq1 and seq2 with <= max_mismatches
        2. Extend each seed match bidirectionally while staying within mismatch limit
        3. Consolidate overlapping matches, keeping the longest ones

        Parameters:
            sequence1: First DNA sequence (e.g., flanking region)
            sequence2: Second DNA sequence (e.g., non-coding region)

        Returns:
            List of alignment dictionaries, each containing:
                - pos1: Position in sequence1 (0-based)
                - pos2: Position in sequence2 (0-based)
                - length: Length of the alignment
                - seq1: Sequence from sequence1
                - seq2: Sequence from sequence2 (as it appears in sequence2)
                - mismatches: Number of mismatches
                - mismatch_positions: List of positions with mismatches
                - orientation: 'forward' or 'reverse_complement'

        Note:
            - min_distance and max_distance are ignored (not applicable between different sequences)
            - pos1 refers to position in sequence1
            - pos2 refers to position in sequence2

        Example:
            >>> finder = ShortAlignmentFinder(min_length=9, max_mismatches=1)
            >>> seq1 = "ATGCACGTACGTACGT"  # Flanking region
            >>> seq2 = "NNNNACGTACGTACGTNNNN"  # Non-coding region
            >>> results = finder.find_alignments_between(seq1, seq2)
            >>> results[0]
            {
                'pos1': 4,
                'pos2': 4,
                'length': 12,
                'seq1': 'ACGTACGTACGT',
                'seq2': 'ACGTACGTACGT',
                'mismatches': 0,
                'mismatch_positions': [],
                'orientation': 'forward'
            }
        """
        sequence1 = sequence1.upper()
        sequence2 = sequence2.upper()

        # Find all seed matches between the two sequences
        seed_matches = self._find_seed_matches_between(sequence1, sequence2)

        # Extend each seed match to maximum length
        extended_matches = []
        for seed in seed_matches:
            extended = self._extend_match_between(sequence1, sequence2, seed)
            if extended:
                extended_matches.append(extended)

        # Consolidate overlapping matches (keep longest)
        consolidated = self._consolidate_matches_between(extended_matches)

        # Sort by pos1, then pos2
        consolidated.sort(key=lambda x: (x['pos1'], x['pos2']))

        return consolidated

    def _find_seed_matches_between(self, sequence1: str, sequence2: str) -> List[Dict]:
        """
        Find all seed matches of min_length between two sequences.

        Returns:
            List of seed match dictionaries
        """
        seeds = []
        len1 = len(sequence1)
        len2 = len(sequence2)

        # For each starting position in sequence1
        for i in range(len1 - self.min_length + 1):
            kmer = sequence1[i:i + self.min_length]

            # Search forward matches in sequence2
            if self.check_forward:
                for j in range(len2 - self.min_length + 1):
                    target = sequence2[j:j + self.min_length]
                    mismatches, mismatch_pos = self._count_mismatches(kmer, target)

                    if mismatches <= self.max_mismatches:
                        seeds.append({
                            'pos1': i,
                            'pos2': j,
                            'length': self.min_length,
                            'orientation': 'forward',
                            'mismatches': mismatches,
                            'mismatch_positions': mismatch_pos
                        })

            # Search reverse complement matches in sequence2
            if self.check_revcomp:
                kmer_rc = reverse_complement(kmer)
                for j in range(len2 - self.min_length + 1):
                    target = sequence2[j:j + self.min_length]
                    mismatches, mismatch_pos = self._count_mismatches(kmer_rc, target)

                    if mismatches <= self.max_mismatches:
                        seeds.append({
                            'pos1': i,
                            'pos2': j,
                            'length': self.min_length,
                            'orientation': 'reverse_complement',
                            'mismatches': mismatches,
                            'mismatch_positions': mismatch_pos
                        })

        return seeds

    def _extend_match_between(self, sequence1: str, sequence2: str, seed: Dict) -> Optional[Dict]:
        """
        Extend a seed match between two sequences bidirectionally.

        Parameters:
            sequence1: First sequence
            sequence2: Second sequence
            seed: Seed match dictionary

        Returns:
            Extended match dictionary, or None if extension fails
        """
        pos1 = seed['pos1']
        pos2 = seed['pos2']
        length = seed['length']
        orientation = seed['orientation']

        # Extend right
        while True:
            # Check boundaries
            if pos1 + length >= len(sequence1) or pos2 + length >= len(sequence2):
                break

            # Get extended sequences and check mismatch limit
            test_seq1 = sequence1[pos1:pos1 + length + 1]
            test_seq2_raw = sequence2[pos2:pos2 + length + 1]

            if orientation == 'forward':
                test_seq2 = test_seq2_raw
            else:
                test_seq2 = reverse_complement(test_seq2_raw)

            test_mismatches, _ = self._count_mismatches(test_seq1, test_seq2)

            if test_mismatches <= self.max_mismatches:
                length += 1
            else:
                break

        # Extend left
        while True:
            # Check boundaries
            if pos1 == 0 or pos2 == 0:
                break

            # Get extended sequences and check mismatch limit
            test_seq1 = sequence1[pos1 - 1:pos1 + length]
            test_seq2_raw = sequence2[pos2 - 1:pos2 + length]

            if orientation == 'forward':
                test_seq2 = test_seq2_raw
            else:
                test_seq2 = reverse_complement(test_seq2_raw)

            test_mismatches, _ = self._count_mismatches(test_seq1, test_seq2)

            if test_mismatches <= self.max_mismatches:
                pos1 -= 1
                pos2 -= 1
                length += 1
            else:
                break

        # Final check: ensure minimum length is met
        if length < self.min_length:
            return None

        # Get final sequences and mismatch info
        seq1 = sequence1[pos1:pos1 + length]
        seq2_raw = sequence2[pos2:pos2 + length]

        if orientation == 'forward':
            seq2_compare = seq2_raw
        else:
            seq2_compare = reverse_complement(seq2_raw)

        mismatches, mismatch_positions = self._count_mismatches(seq1, seq2_compare)

        return {
            'pos1': pos1,
            'pos2': pos2,
            'length': length,
            'seq1': seq1,
            'seq2': seq2_raw,  # Store as it appears in sequence2
            'mismatches': mismatches,
            'mismatch_positions': mismatch_positions,
            'orientation': orientation
        }

    def _consolidate_matches_between(self, matches: List[Dict]) -> List[Dict]:
        """
        Consolidate overlapping matches between two sequences.

        Similar to _consolidate_matches but for pairwise comparisons.
        """
        if not matches:
            return []

        # Sort by length (descending), then by mismatches (ascending)
        matches_sorted = sorted(matches,
                               key=lambda x: (-x['length'], x['mismatches']))

        consolidated = []

        for match in matches_sorted:
            # Check if this match overlaps with any already consolidated match
            is_redundant = False

            for existing in consolidated:
                if self._matches_overlap_between(match, existing):
                    is_redundant = True
                    break

            if not is_redundant:
                consolidated.append(match)

        return consolidated

    def _matches_overlap_between(self, match1: Dict, match2: Dict) -> bool:
        """
        Check if two matches between sequences overlap.

        Matches overlap if they have the same orientation and their regions overlap
        in both sequence1 and sequence2.
        """
        # Must be same orientation to overlap
        if match1['orientation'] != match2['orientation']:
            return False

        # Check if regions overlap in sequence1
        seq1_overlap = self._regions_overlap(
            match1['pos1'], match1['pos1'] + match1['length'],
            match2['pos1'], match2['pos1'] + match2['length']
        )

        # Check if regions overlap in sequence2
        seq2_overlap = self._regions_overlap(
            match1['pos2'], match1['pos2'] + match1['length'],
            match2['pos2'], match2['pos2'] + match2['length']
        )

        # Both must overlap for matches to be considered overlapping
        return seq1_overlap and seq2_overlap

    def extend_alignment_with_gaps(self,
                                    sequence1: str,
                                    sequence2: str,
                                    seed_match: Dict,
                                    min_identity: float = 0.80,
                                    max_extension: int = 50,
                                    match_score: int = 2,
                                    mismatch_penalty: int = -1,
                                    gap_open_penalty: int = -3,
                                    gap_extend_penalty: int = -1) -> Optional[Dict]:
        """
        Extend an alignment with gaps allowed (for CAST homing spacer-like patterns).

        Takes an existing ungapped alignment as a seed and extends it bidirectionally
        using semi-global alignment, allowing insertions and deletions. Useful for
        finding patterns like 17/20 or 27/32 matches with gaps.

        Algorithm:
        1. Use seed alignment as anchor point
        2. Extend 5' (left) and 3' (right) using dynamic programming
        3. Allow gaps with configurable penalties
        4. Track alignment path and calculate identity
        5. Return if meets minimum identity threshold

        Parameters:
            sequence1: First DNA sequence (or full sequence if from find_alignments)
            sequence2: Second DNA sequence (or full sequence if from find_alignments)
            seed_match: Existing alignment dictionary (from find_alignments or find_alignments_between)
                Must contain: 'pos1', 'pos2', 'length', 'orientation'
            min_identity: Minimum identity fraction (matches/total_length) required (default: 0.80)
            max_extension: Maximum bases to extend in each direction (default: 50)
            match_score: Score for matching bases (default: 2)
            mismatch_penalty: Penalty for mismatches (default: -1)
            gap_open_penalty: Penalty for opening a gap (default: -3)
            gap_extend_penalty: Penalty for extending a gap (default: -1)

        Returns:
            Dictionary containing:
                - pos1: Start position in sequence1 (0-based)
                - pos2: Start position in sequence2 (0-based)
                - end1: End position in sequence1 (exclusive)
                - end2: End position in sequence2 (exclusive)
                - length1: Length in sequence1 (bases consumed)
                - length2: Length in sequence2 (bases consumed)
                - seq1_aligned: Aligned sequence from seq1 with gaps ('-')
                - seq2_aligned: Aligned sequence from seq2 with gaps ('-')
                - alignment_length: Total alignment length (including gaps)
                - matches: Number of matching positions
                - mismatches: Number of mismatching positions
                - gaps: Number of gap positions
                - identity: Identity fraction (matches / alignment_length)
                - alignment_string: Visual alignment ('|' match, '.' mismatch, ' ' gap)
                - orientation: 'forward' or 'reverse_complement'
                - seed_pos1: Original seed start in sequence1
                - seed_pos2: Original seed start in sequence2
                - seed_length: Original seed length
            Or None if extension doesn't meet identity threshold

        Example:
            >>> finder = ShortAlignmentFinder(min_length=9, max_mismatches=0)
            >>> seq1 = "ACGTACGTACGTACGT"
            >>> seq2 = "ACGTACGTNNNACGTACGT"  # Has 3bp insertion
            >>> # First find ungapped seed
            >>> seeds = finder.find_alignments_between(seq1, seq2)
            >>> # Then extend with gaps
            >>> gapped = finder.extend_alignment_with_gaps(seq1, seq2, seeds[0])
            >>> print(f"Identity: {gapped['identity']:.2%}")
            >>> print(f"Pattern: {gapped['matches']}/{gapped['alignment_length']}")
        """
        sequence1 = sequence1.upper()
        sequence2 = sequence2.upper()

        # Extract seed information
        seed_pos1 = seed_match['pos1']
        seed_pos2 = seed_match['pos2']
        seed_length = seed_match['length']
        orientation = seed_match.get('orientation', 'forward')

        # For reverse complement, we need to work with the RC of sequence2
        if orientation == 'reverse_complement':
            # Get the region from sequence2 and reverse complement it
            # We'll work in the forward orientation internally
            working_seq2 = reverse_complement(sequence2)
            # Adjust pos2 to the RC coordinate system
            working_pos2 = len(sequence2) - seed_pos2 - seed_length
        else:
            working_seq2 = sequence2
            working_pos2 = seed_pos2

        # Extend left from seed
        left_result = self._extend_gapped_left(
            sequence1, working_seq2,
            seed_pos1, working_pos2,
            max_extension,
            match_score, mismatch_penalty,
            gap_open_penalty, gap_extend_penalty
        )

        # Extend right from seed
        right_result = self._extend_gapped_right(
            sequence1, working_seq2,
            seed_pos1 + seed_length, working_pos2 + seed_length,
            max_extension,
            match_score, mismatch_penalty,
            gap_open_penalty, gap_extend_penalty
        )

        # Get seed sequences
        seed_seq1 = sequence1[seed_pos1:seed_pos1 + seed_length]
        if orientation == 'forward':
            seed_seq2 = sequence2[seed_pos2:seed_pos2 + seed_length]
        else:
            seed_seq2 = reverse_complement(sequence2[seed_pos2:seed_pos2 + seed_length])

        # Combine: left extension + seed + right extension
        full_seq1 = left_result['seq1'] + seed_seq1 + right_result['seq1']
        full_seq2 = left_result['seq2'] + seed_seq2 + right_result['seq2']

        # Calculate final positions
        final_pos1 = seed_pos1 - left_result['consumed1']
        final_end1 = seed_pos1 + seed_length + right_result['consumed1']

        if orientation == 'forward':
            final_pos2 = seed_pos2 - left_result['consumed2']
            final_end2 = seed_pos2 + seed_length + right_result['consumed2']
        else:
            # Map back to original sequence2 coordinates
            final_pos2 = len(sequence2) - (working_pos2 + right_result['consumed2'])
            final_end2 = len(sequence2) - (working_pos2 - left_result['consumed2'])

        # Calculate alignment statistics
        matches = 0
        mismatches = 0
        gaps = 0
        alignment_string = ""

        for b1, b2 in zip(full_seq1, full_seq2):
            if b1 == '-' or b2 == '-':
                gaps += 1
                alignment_string += ' '
            elif b1 == b2:
                matches += 1
                alignment_string += '|'
            else:
                mismatches += 1
                alignment_string += '.'

        alignment_length = len(full_seq1)
        identity = matches / alignment_length if alignment_length > 0 else 0.0

        # Check minimum identity threshold
        if identity < min_identity:
            return None

        # Convert seq2 back to original orientation if needed
        if orientation == 'reverse_complement':
            display_seq2 = reverse_complement(full_seq2)
        else:
            display_seq2 = full_seq2

        return {
            'pos1': final_pos1,
            'pos2': final_pos2,
            'end1': final_end1,
            'end2': final_end2,
            'length1': final_end1 - final_pos1,
            'length2': final_end2 - final_pos2,
            'seq1_aligned': full_seq1,
            'seq2_aligned': display_seq2,
            'alignment_length': alignment_length,
            'matches': matches,
            'mismatches': mismatches,
            'gaps': gaps,
            'identity': identity,
            'alignment_string': alignment_string,
            'orientation': orientation,
            'seed_pos1': seed_pos1,
            'seed_pos2': seed_pos2,
            'seed_length': seed_length
        }

    def _extend_gapped_left(self, seq1: str, seq2: str,
                           pos1: int, pos2: int,
                           max_extension: int,
                           match_score: int,
                           mismatch_penalty: int,
                           gap_open: int,
                           gap_extend: int) -> Dict:
        """
        Extend alignment leftward with gaps using dynamic programming.

        Returns:
            Dictionary with 'seq1', 'seq2', 'consumed1', 'consumed2'
        """
        # Extract regions to align (reverse them for left extension)
        extend1 = seq1[max(0, pos1 - max_extension):pos1][::-1]
        extend2 = seq2[max(0, pos2 - max_extension):pos2][::-1]

        if not extend1 or not extend2:
            return {'seq1': '', 'seq2': '', 'consumed1': 0, 'consumed2': 0}

        # Perform semi-global alignment
        aligned1, aligned2 = self._align_sequences(
            extend1, extend2,
            match_score, mismatch_penalty,
            gap_open, gap_extend
        )

        # Reverse back to original orientation
        aligned1_fwd = aligned1[::-1]
        aligned2_fwd = aligned2[::-1]

        # Count consumed bases (non-gap characters)
        consumed1 = sum(1 for c in aligned1 if c != '-')
        consumed2 = sum(1 for c in aligned2 if c != '-')

        return {
            'seq1': aligned1_fwd,
            'seq2': aligned2_fwd,
            'consumed1': consumed1,
            'consumed2': consumed2
        }

    def _extend_gapped_right(self, seq1: str, seq2: str,
                            pos1: int, pos2: int,
                            max_extension: int,
                            match_score: int,
                            mismatch_penalty: int,
                            gap_open: int,
                            gap_extend: int) -> Dict:
        """
        Extend alignment rightward with gaps using dynamic programming.

        Returns:
            Dictionary with 'seq1', 'seq2', 'consumed1', 'consumed2'
        """
        # Extract regions to align
        extend1 = seq1[pos1:min(len(seq1), pos1 + max_extension)]
        extend2 = seq2[pos2:min(len(seq2), pos2 + max_extension)]

        if not extend1 or not extend2:
            return {'seq1': '', 'seq2': '', 'consumed1': 0, 'consumed2': 0}

        # Perform semi-global alignment
        aligned1, aligned2 = self._align_sequences(
            extend1, extend2,
            match_score, mismatch_penalty,
            gap_open, gap_extend
        )

        # Count consumed bases (non-gap characters)
        consumed1 = sum(1 for c in aligned1 if c != '-')
        consumed2 = sum(1 for c in aligned2 if c != '-')

        return {
            'seq1': aligned1,
            'seq2': aligned2,
            'consumed1': consumed1,
            'consumed2': consumed2
        }

    def _align_sequences(self, seq1: str, seq2: str,
                        match_score: int,
                        mismatch_penalty: int,
                        gap_open: int,
                        gap_extend: int) -> Tuple[str, str]:
        """
        Align two sequences using semi-global alignment with affine gap penalties.

        Uses a simplified Needleman-Wunsch algorithm with affine gap penalties.

        Parameters:
            seq1, seq2: Sequences to align
            match_score: Score for match
            mismatch_penalty: Penalty for mismatch
            gap_open: Gap opening penalty
            gap_extend: Gap extension penalty

        Returns:
            Tuple of (aligned_seq1, aligned_seq2) with gaps as '-'
        """
        len1 = len(seq1)
        len2 = len(seq2)

        # Initialize scoring matrices (main, gap in seq1, gap in seq2)
        M = [[0] * (len2 + 1) for _ in range(len1 + 1)]  # Match/mismatch
        Ix = [[float('-inf')] * (len2 + 1) for _ in range(len1 + 1)]  # Gap in seq2
        Iy = [[float('-inf')] * (len2 + 1) for _ in range(len1 + 1)]  # Gap in seq1

        # Initialize first row and column
        for i in range(len1 + 1):
            M[i][0] = 0  # Semi-global: no penalty for starting gaps
            Ix[i][0] = float('-inf')
            Iy[i][0] = gap_open + i * gap_extend if i > 0 else 0

        for j in range(len2 + 1):
            M[0][j] = 0  # Semi-global: no penalty for starting gaps
            Ix[0][j] = gap_open + j * gap_extend if j > 0 else 0
            Iy[0][j] = float('-inf')

        # Fill matrices
        for i in range(1, len1 + 1):
            for j in range(1, len2 + 1):
                # Calculate match/mismatch score
                if seq1[i-1] == seq2[j-1]:
                    match = match_score
                else:
                    match = mismatch_penalty

                # M[i][j]: align seq1[i-1] with seq2[j-1]
                M[i][j] = max(
                    M[i-1][j-1] + match,
                    Ix[i-1][j-1] + match,
                    Iy[i-1][j-1] + match
                )

                # Ix[i][j]: gap in seq2 (seq1[i-1] aligns with gap)
                Ix[i][j] = max(
                    M[i-1][j] + gap_open + gap_extend,
                    Ix[i-1][j] + gap_extend
                )

                # Iy[i][j]: gap in seq1 (seq2[j-1] aligns with gap)
                Iy[i][j] = max(
                    M[i][j-1] + gap_open + gap_extend,
                    Iy[i][j-1] + gap_extend
                )

        # Traceback to get alignment
        aligned1 = []
        aligned2 = []
        i, j = len1, len2

        # Determine which matrix we're in
        current = 'M'  # Start from match matrix
        max_score = M[i][j]
        if Ix[i][j] > max_score:
            current = 'Ix'
            max_score = Ix[i][j]
        if Iy[i][j] > max_score:
            current = 'Iy'

        while i > 0 or j > 0:
            if current == 'M':
                if i == 0 or j == 0:
                    break
                aligned1.append(seq1[i-1])
                aligned2.append(seq2[j-1])

                # Find where we came from
                if seq1[i-1] == seq2[j-1]:
                    match = match_score
                else:
                    match = mismatch_penalty

                prev_score = M[i][j] - match
                if abs(M[i-1][j-1] - prev_score) < 0.01:
                    current = 'M'
                elif abs(Ix[i-1][j-1] - prev_score) < 0.01:
                    current = 'Ix'
                else:
                    current = 'Iy'

                i -= 1
                j -= 1

            elif current == 'Ix':
                # Gap in seq2
                if i == 0:
                    break
                aligned1.append(seq1[i-1])
                aligned2.append('-')

                # Check where we came from
                if abs(Ix[i-1][j] - (Ix[i][j] - gap_extend)) < 0.01:
                    current = 'Ix'
                else:
                    current = 'M'

                i -= 1

            elif current == 'Iy':
                # Gap in seq1
                if j == 0:
                    break
                aligned1.append('-')
                aligned2.append(seq2[j-1])

                # Check where we came from
                if abs(Iy[i][j-1] - (Iy[i][j] - gap_extend)) < 0.01:
                    current = 'Iy'
                else:
                    current = 'M'

                j -= 1

        # Handle remaining positions (semi-global alignment)
        while i > 0:
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        while j > 0:
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1

        # Reverse to get correct order
        return ''.join(reversed(aligned1)), ''.join(reversed(aligned2))


def find_short_alignments(sequence: str,
                          min_length: int = 9,
                          max_mismatches: int = 0,
                          min_distance: int = 0,
                          max_distance: Optional[int] = None,
                          check_forward: bool = True,
                          check_revcomp: bool = True) -> List[Dict]:
    """
    Convenience function to find short alignments without creating a class instance.

    Parameters:
        sequence: DNA sequence to search
        min_length: Minimum alignment length (default: 9bp)
        max_mismatches: Maximum mismatches allowed (default: 0)
        min_distance: Minimum distance between positions (default: 0)
        max_distance: Maximum distance between positions (default: None)
        check_forward: Find direct repeats (default: True)
        check_revcomp: Find inverted repeats (default: True)

    Returns:
        List of alignment dictionaries

    Example:
        >>> results = find_short_alignments("ACGTACGTNNNNACGTACGT", min_length=8)
        >>> print(f"Found {len(results)} alignments")
    """
    finder = ShortAlignmentFinder(
        min_length=min_length,
        max_mismatches=max_mismatches,
        min_distance=min_distance,
        max_distance=max_distance,
        check_forward=check_forward,
        check_revcomp=check_revcomp
    )
    return finder.find_alignments(sequence)


def find_alignments_between_sequences(sequence1: str,
                                      sequence2: str,
                                      min_length: int = 9,
                                      max_mismatches: int = 0,
                                      check_forward: bool = True,
                                      check_revcomp: bool = True) -> List[Dict]:
    """
    Convenience function to find alignments between two different sequences.

    Useful for finding matches between:
    - Flanking region and non-coding region of transposon
    - Upstream region and downstream region
    - Any two distinct genomic regions

    Parameters:
        sequence1: First DNA sequence (e.g., flanking region)
        sequence2: Second DNA sequence (e.g., non-coding region)
        min_length: Minimum alignment length (default: 9bp)
        max_mismatches: Maximum mismatches allowed (default: 0)
        check_forward: Find direct matches (default: True)
        check_revcomp: Find reverse complement matches (default: True)

    Returns:
        List of alignment dictionaries, each containing:
            - pos1: Position in sequence1 (0-based)
            - pos2: Position in sequence2 (0-based)
            - length: Length of the alignment
            - seq1: Sequence from sequence1
            - seq2: Sequence from sequence2
            - mismatches: Number of mismatches
            - mismatch_positions: List of mismatch positions
            - orientation: 'forward' or 'reverse_complement'

    Example:
        >>> seq1 = "ATGCACGTACGTACGT"  # Flanking region
        >>> seq2 = "NNNNACGTACGTACGTNNNN"  # Non-coding region
        >>> results = find_alignments_between_sequences(seq1, seq2, min_length=9)
        >>> print(f"Found {len(results)} alignments between the two sequences")
    """
    # Note: min_distance and max_distance are not applicable for between-sequence searches
    finder = ShortAlignmentFinder(
        min_length=min_length,
        max_mismatches=max_mismatches,
        min_distance=0,  # Not used in between-sequence search
        max_distance=None,  # Not used in between-sequence search
        check_forward=check_forward,
        check_revcomp=check_revcomp
    )
    return finder.find_alignments_between(sequence1, sequence2)


def extend_alignment_with_gaps(sequence1: str,
                               sequence2: str,
                               seed_match: Dict,
                               min_identity: float = 0.80,
                               max_extension: int = 50,
                               match_score: int = 2,
                               mismatch_penalty: int = -1,
                               gap_open_penalty: int = -3,
                               gap_extend_penalty: int = -1) -> Optional[Dict]:
    """
    Convenience function to extend an alignment with gaps allowed.

    Takes an existing ungapped alignment and extends it bidirectionally with
    gaps allowed. Useful for finding CAST homing spacer-like patterns
    (e.g., 17/20 or 27/32 matches with gaps/indels).

    Parameters:
        sequence1: First DNA sequence
        sequence2: Second DNA sequence
        seed_match: Existing alignment dictionary from find_alignments or find_alignments_between
        min_identity: Minimum identity fraction required (default: 0.80)
        max_extension: Maximum bases to extend in each direction (default: 50)
        match_score: Score for matching bases (default: 2)
        mismatch_penalty: Penalty for mismatches (default: -1)
        gap_open_penalty: Penalty for opening a gap (default: -3)
        gap_extend_penalty: Penalty for extending a gap (default: -1)

    Returns:
        Dictionary with alignment details including gaps, or None if below identity threshold

    Example:
        >>> # First find ungapped seeds
        >>> seeds = find_alignments_between_sequences(seq1, seq2, min_length=9)
        >>> # Then extend best seed with gaps
        >>> if seeds:
        ...     gapped = extend_alignment_with_gaps(seq1, seq2, seeds[0], min_identity=0.85)
        ...     if gapped:
        ...         print(f"Found {gapped['matches']}/{gapped['alignment_length']} match")
        ...         print(f"Identity: {gapped['identity']:.1%}")
    """
    finder = ShortAlignmentFinder(min_length=9)  # Default params, not used for extension
    return finder.extend_alignment_with_gaps(
        sequence1=sequence1,
        sequence2=sequence2,
        seed_match=seed_match,
        min_identity=min_identity,
        max_extension=max_extension,
        match_score=match_score,
        mismatch_penalty=mismatch_penalty,
        gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty
    )
