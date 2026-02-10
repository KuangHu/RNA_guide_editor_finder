"""Tests for alignment_filter module - statistical significance features."""

import math
import sys
import os
import pytest

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from modules.alignment_filter import (
    calculate_ungapped_evalue,
    calculate_gapped_evalue,
    compute_lambda_for_scoring,
    evalue_to_pvalue,
    count_gap_opens,
    Filters,
    FilterEngine,
    FilterPipeline,
    FilterLevel,
    create_default_pipeline,
    create_strict_pipeline,
    create_relaxed_pipeline,
    create_legacy_pipeline,
)
from modules.short_alignment_finder import (
    find_alignments_between_sequences,
    extend_alignment_with_gaps,
)


# ============================================
# Unit tests: calculate_ungapped_evalue
# ============================================

class TestCalculateUngappedEvalue:
    """Tests for the ungapped E-value calculation."""

    def test_8bp_0mm(self):
        """8bp perfect match with G=500 should give E ~ 0.0076."""
        ev = calculate_ungapped_evalue(8, 0, 500)
        assert abs(ev - 0.00763) < 0.001

    def test_9bp_0mm(self):
        """9bp perfect match with G=500 should give E ~ 0.0019."""
        ev = calculate_ungapped_evalue(9, 0, 500)
        assert abs(ev - 0.00191) < 0.001

    def test_9bp_1mm(self):
        """9bp with 1 mismatch should give E ~ 0.052."""
        ev = calculate_ungapped_evalue(9, 1, 500)
        assert abs(ev - 0.0515) < 0.005

    def test_12bp_0mm(self):
        """12bp perfect match should give E ~ 2.98e-5."""
        ev = calculate_ungapped_evalue(12, 0, 500)
        assert abs(ev - 2.98e-5) < 1e-6

    def test_15bp_0mm(self):
        """15bp perfect match should give E ~ 4.66e-7."""
        ev = calculate_ungapped_evalue(15, 0, 500)
        assert abs(ev - 4.66e-7) < 1e-8

    def test_17bp_1mm(self):
        """17bp with 1 mismatch should give E ~ 1.48e-6."""
        ev = calculate_ungapped_evalue(17, 1, 500)
        assert abs(ev - 1.484e-6) < 1e-7

    def test_20bp_2mm(self):
        """20bp with 2 mismatches should give E ~ 7.8e-7."""
        ev = calculate_ungapped_evalue(20, 2, 500)
        assert abs(ev - 7.8e-7) < 2e-7

    def test_evalue_decreases_with_length(self):
        """Longer alignments should have smaller E-values (more significant)."""
        ev_short = calculate_ungapped_evalue(10, 0, 500)
        ev_long = calculate_ungapped_evalue(20, 0, 500)
        assert ev_long < ev_short

    def test_evalue_increases_with_mismatches(self):
        """More mismatches should give larger E-values (less significant)."""
        ev_perfect = calculate_ungapped_evalue(15, 0, 500)
        ev_1mm = calculate_ungapped_evalue(15, 1, 500)
        ev_2mm = calculate_ungapped_evalue(15, 2, 500)
        assert ev_perfect < ev_1mm < ev_2mm

    def test_evalue_scales_with_search_space(self):
        """E-value should scale linearly with search space."""
        ev_small = calculate_ungapped_evalue(10, 0, 100)
        ev_large = calculate_ungapped_evalue(10, 0, 1000)
        assert abs(ev_large / ev_small - 10.0) < 0.01

    def test_invalid_mismatches_more_than_length(self):
        """Mismatches > length should return infinity."""
        ev = calculate_ungapped_evalue(5, 6, 500)
        assert ev == float('inf')

    def test_negative_mismatches(self):
        """Negative mismatches should return infinity."""
        ev = calculate_ungapped_evalue(5, -1, 500)
        assert ev == float('inf')

    def test_zero_length(self):
        """Zero-length alignment should return infinity."""
        ev = calculate_ungapped_evalue(0, 0, 500)
        assert ev == float('inf')


# ============================================
# Unit tests: compute_lambda_for_scoring
# ============================================

class TestComputeLambda:
    """Tests for the Newton-Raphson lambda solver."""

    def test_default_scoring_lambda(self):
        """For match=2, mismatch=-1, lambda ~ 0.2645 (satisfies the equation)."""
        lam = compute_lambda_for_scoring(2, -1)
        assert abs(lam - 0.2645) < 1e-3

    def test_lambda_is_positive(self):
        """Lambda should always be positive."""
        lam = compute_lambda_for_scoring(1, -1)
        assert lam > 0

    def test_lambda_satisfies_equation(self):
        """Lambda should satisfy 0.25*exp(s*lam) + 0.75*exp(r*lam) = 1."""
        s, r = 2, -1
        lam = compute_lambda_for_scoring(s, r)
        result = 0.25 * math.exp(s * lam) + 0.75 * math.exp(r * lam)
        assert abs(result - 1.0) < 1e-8


# ============================================
# Unit tests: evalue_to_pvalue
# ============================================

class TestEvalueToPvalue:
    """Tests for E-value to P-value conversion."""

    def test_zero_evalue(self):
        """E=0 should give P=0."""
        assert evalue_to_pvalue(0) == 0.0

    def test_small_evalue(self):
        """For small E, P ~ E."""
        p = evalue_to_pvalue(0.001)
        assert abs(p - 0.001) < 1e-5

    def test_large_evalue(self):
        """Large E should give P close to 1."""
        p = evalue_to_pvalue(10)
        assert p > 0.999

    def test_very_large_evalue(self):
        """Very large E should give P = 1.0 without overflow."""
        p = evalue_to_pvalue(1000)
        assert p == 1.0

    def test_pvalue_between_0_and_1(self):
        """P-value should always be in [0, 1]."""
        for ev in [0, 0.001, 0.01, 0.1, 1, 10, 100]:
            p = evalue_to_pvalue(ev)
            assert 0 <= p <= 1


# ============================================
# Unit tests: count_gap_opens
# ============================================

class TestCountGapOpens:
    """Tests for gap-open counting."""

    def test_no_gaps(self):
        assert count_gap_opens("ACGTACGT") == 0

    def test_single_gap(self):
        assert count_gap_opens("ACG ACGT") == 1

    def test_multiple_gaps(self):
        assert count_gap_opens("AC GT AC") == 2

    def test_consecutive_gaps_count_as_one(self):
        assert count_gap_opens("AC   GT") == 1

    def test_empty_string(self):
        assert count_gap_opens("") == 0

    def test_all_gaps(self):
        assert count_gap_opens("   ") == 1


# ============================================
# Unit tests: calculate_gapped_evalue
# ============================================

class TestCalculateGappedEvalue:
    """Tests for the Karlin-Altschul gapped E-value calculation."""

    def test_positive_score(self):
        """Positive score should give a finite E-value."""
        ev = calculate_gapped_evalue(20, 100, 500, 0.1, math.log(2))
        assert 0 < ev < float('inf')

    def test_higher_score_lower_evalue(self):
        """Higher scores should give lower E-values."""
        ev_low = calculate_gapped_evalue(10, 100, 500, 0.1, math.log(2))
        ev_high = calculate_gapped_evalue(20, 100, 500, 0.1, math.log(2))
        assert ev_high < ev_low

    def test_zero_score(self):
        """Score of 0 should return infinity."""
        ev = calculate_gapped_evalue(0, 100, 500, 0.1, math.log(2))
        assert ev == float('inf')

    def test_negative_score(self):
        """Negative score should return infinity."""
        ev = calculate_gapped_evalue(-5, 100, 500, 0.1, math.log(2))
        assert ev == float('inf')

    def test_evalue_scales_with_search_space(self):
        """E-value should scale with target_len."""
        ev1 = calculate_gapped_evalue(15, 100, 500, 0.1, math.log(2))
        ev2 = calculate_gapped_evalue(15, 100, 5000, 0.1, math.log(2))
        assert abs(ev2 / ev1 - 10.0) < 0.01


# ============================================
# Integration tests: statistical_significance_filter
# ============================================

class TestStatisticalSignificanceFilter:
    """Integration tests for the filter function."""

    def _make_alignment(self, length, mismatches=0, gaps=0, alignment_string=''):
        """Helper to create a mock alignment dict."""
        seq = 'A' * length
        return seq, {
            'aligned_sequence': seq,
            'mismatches': mismatches,
            'gaps': gaps,
            'alignment_string': alignment_string,
        }

    def test_9bp_0mm_passes_as_low(self):
        """9bp perfect match should pass with 'low' confidence."""
        seq, alignment = self._make_alignment(9, mismatches=0)
        params = {'search_space': 500, 'e_value_reject': 0.01,
                  'e_value_low': 1e-3, 'e_value_high': 1e-6}
        passed, reason, metrics = Filters.statistical_significance_filter(
            seq, params, alignment=alignment)
        assert passed is True
        assert metrics['confidence'] == 'low'

    def test_15bp_0mm_passes_as_high(self):
        """15bp perfect match should pass with 'high' confidence."""
        seq, alignment = self._make_alignment(15, mismatches=0)
        params = {'search_space': 500, 'e_value_reject': 0.01,
                  'e_value_low': 1e-3, 'e_value_high': 1e-6}
        passed, reason, metrics = Filters.statistical_significance_filter(
            seq, params, alignment=alignment)
        assert passed is True
        assert metrics['confidence'] == 'high'

    def test_9bp_1mm_rejected(self):
        """9bp with 1 mismatch should be rejected."""
        seq, alignment = self._make_alignment(9, mismatches=1)
        params = {'search_space': 500, 'e_value_reject': 0.01,
                  'e_value_low': 1e-3, 'e_value_high': 1e-6}
        passed, reason, metrics = Filters.statistical_significance_filter(
            seq, params, alignment=alignment)
        assert passed is False
        assert metrics['confidence'] == 'rejected'
        assert 'E_VALUE_TOO_HIGH' in reason

    def test_8bp_0mm_passes_as_low(self):
        """8bp perfect match: E ~ 0.0076, within (1e-3, 0.01] -> low."""
        seq, alignment = self._make_alignment(8, mismatches=0)
        params = {'search_space': 500, 'e_value_reject': 0.01,
                  'e_value_low': 1e-3, 'e_value_high': 1e-6}
        passed, reason, metrics = Filters.statistical_significance_filter(
            seq, params, alignment=alignment)
        assert passed is True
        assert metrics['confidence'] == 'low'

    def test_no_alignment_defaults(self):
        """Without alignment dict, should assume 0 mismatches and 0 gaps."""
        seq = 'ACGTACGTACGT'  # 12bp
        params = {'search_space': 500, 'e_value_reject': 0.01,
                  'e_value_low': 1e-3, 'e_value_high': 1e-6}
        passed, reason, metrics = Filters.statistical_significance_filter(
            seq, params, alignment=None)
        assert passed is True

    def test_transposon_data_overrides_search_space(self):
        """search_space from transposon_data should override params default."""
        seq, alignment = self._make_alignment(9, mismatches=0)
        params = {'search_space': 500, 'e_value_reject': 0.01,
                  'e_value_low': 1e-3, 'e_value_high': 1e-6}

        # With much larger search space, E-value goes up
        transposon_data = {'search_space': 50000}
        passed_large, _, metrics_large = Filters.statistical_significance_filter(
            seq, params, alignment=alignment, transposon_data=transposon_data)

        # Compare with default search space
        passed_default, _, metrics_default = Filters.statistical_significance_filter(
            seq, params, alignment=alignment)

        assert metrics_large['e_value'] > metrics_default['e_value']

    def test_gapped_alignment(self):
        """Gapped alignment should use Karlin-Altschul path."""
        seq = 'ACGTACGTACGTACGT'  # 16bp
        alignment = {
            'aligned_sequence': seq,
            'mismatches': 1,
            'gaps': 2,
            'alignment_string': 'ACGT ACGTACG ACGT',
        }
        params = {
            'search_space': 500,
            'e_value_reject': 0.01,
            'e_value_low': 1e-3,
            'e_value_high': 1e-6,
            'match_score': 2,
            'mismatch_penalty': -1,
            'gap_open_penalty': -3,
            'gap_extend_penalty': -1,
            'karlin_K': 0.1,
            'karlin_lambda': None,
        }
        passed, reason, metrics = Filters.statistical_significance_filter(
            seq, params, alignment=alignment)
        assert metrics['is_gapped'] is True
        assert 'raw_score' in metrics
        assert 'gap_opens' in metrics
        assert metrics['e_value'] >= 0

    def test_evalue_in_metrics(self):
        """E-value and P-value should always appear in metrics."""
        seq, alignment = self._make_alignment(12, mismatches=0)
        params = {'search_space': 500, 'e_value_reject': 0.01,
                  'e_value_low': 1e-3, 'e_value_high': 1e-6}
        _, _, metrics = Filters.statistical_significance_filter(
            seq, params, alignment=alignment)
        assert 'e_value' in metrics
        assert 'p_value' in metrics
        assert 'confidence' in metrics

    def test_search_space_from_sequence_lengths(self):
        """Search space should be computed from query_length and target_length."""
        seq = 'A' * 10  # 10bp alignment
        alignment = {
            'aligned_sequence': seq,
            'mismatches': 0,
            'gaps': 0,
            'query_length': 50,   # flanking region = 50bp
            'target_length': 200, # non-coding region = 200bp
        }
        params = {'e_value_reject': 0.01, 'e_value_low': 1e-3, 'e_value_high': 1e-6}
        _, _, metrics = Filters.statistical_significance_filter(
            seq, params, alignment=alignment)
        # G = (50-10+1) * (200-10+1) = 41 * 191 = 7831
        expected_search_space = 41 * 191
        assert metrics['search_space'] == expected_search_space

    def test_sequence_lengths_override_params_search_space(self):
        """query_length/target_length in alignment should take priority over params."""
        seq = 'A' * 10
        alignment = {
            'aligned_sequence': seq,
            'mismatches': 0,
            'gaps': 0,
            'query_length': 50,
            'target_length': 200,
        }
        # Even with search_space=500 in params, sequence lengths should win
        params = {'search_space': 500, 'e_value_reject': 0.01,
                  'e_value_low': 1e-3, 'e_value_high': 1e-6}
        _, _, metrics = Filters.statistical_significance_filter(
            seq, params, alignment=alignment)
        assert metrics['search_space'] == 41 * 191  # Not 500

    def test_fallback_to_params_search_space(self):
        """Without query_length/target_length, should fall back to params."""
        seq = 'A' * 10
        alignment = {
            'aligned_sequence': seq,
            'mismatches': 0,
            'gaps': 0,
        }
        params = {'search_space': 500, 'e_value_reject': 0.01,
                  'e_value_low': 1e-3, 'e_value_high': 1e-6}
        _, _, metrics = Filters.statistical_significance_filter(
            seq, params, alignment=alignment)
        assert metrics['search_space'] == 500

    def test_larger_sequences_give_higher_evalue(self):
        """Longer flanking/non-coding regions should increase E-value."""
        seq = 'A' * 12
        alignment_small = {
            'aligned_sequence': seq, 'mismatches': 0, 'gaps': 0,
            'query_length': 30, 'target_length': 50,
        }
        alignment_large = {
            'aligned_sequence': seq, 'mismatches': 0, 'gaps': 0,
            'query_length': 300, 'target_length': 5000,
        }
        params = {'e_value_reject': 1.0, 'e_value_low': 0.5, 'e_value_high': 0.01}
        _, _, metrics_small = Filters.statistical_significance_filter(
            seq, params, alignment=alignment_small)
        _, _, metrics_large = Filters.statistical_significance_filter(
            seq, params, alignment=alignment_large)
        assert metrics_large['e_value'] > metrics_small['e_value']


# ============================================
# Pipeline tests
# ============================================

class TestPipelines:
    """Tests for pipeline presets."""

    def test_default_pipeline_uses_statistical_significance(self):
        """Default pipeline should include statistical_significance filter."""
        pipeline = create_default_pipeline()
        filter_names = [f.name for f in pipeline.filters]
        assert 'statistical_significance' in filter_names
        assert 'length_hard_cutoff' not in filter_names
        assert 'length_confidence' not in filter_names

    def test_legacy_pipeline_uses_length_filters(self):
        """Legacy pipeline should use original length-based filters."""
        pipeline = create_legacy_pipeline()
        filter_names = [f.name for f in pipeline.filters]
        assert 'length_hard_cutoff' in filter_names
        assert 'length_confidence' in filter_names
        assert 'statistical_significance' not in filter_names

    def test_strict_pipeline_has_tighter_thresholds(self):
        """Strict pipeline should have smaller E-value thresholds."""
        pipeline = create_strict_pipeline()
        for f in pipeline.filters:
            if f.name == 'statistical_significance':
                assert f.params['e_value_reject'] == 1e-3
                assert f.params['e_value_low'] == 1e-6
                assert f.params['e_value_high'] == 1e-9
                break
        else:
            pytest.fail("statistical_significance not found in strict pipeline")

    def test_relaxed_pipeline_has_looser_thresholds(self):
        """Relaxed pipeline should have larger E-value thresholds."""
        pipeline = create_relaxed_pipeline()
        for f in pipeline.filters:
            if f.name == 'statistical_significance':
                assert f.params['e_value_reject'] == 0.05
                assert f.params['e_value_low'] == 0.01
                assert f.params['e_value_high'] == 1e-3
                break
        else:
            pytest.fail("statistical_significance not found in relaxed pipeline")

    def test_default_pipeline_engine_runs(self):
        """Default pipeline should work end-to-end through FilterEngine."""
        pipeline = create_default_pipeline()
        engine = FilterEngine(pipeline)

        data = {
            "test_transposon": {
                "start": 100,
                "end": 200,
                "alignments": [
                    {
                        "aligned_sequence": "ACGTACGTACGTACGT",  # 16bp, 0mm
                        "non_coding_start": 110,
                        "non_coding_end": 126,
                        "mismatches": 0,
                        "gaps": 0,
                    },
                    {
                        "aligned_sequence": "ACGTACGT",  # 8bp, 0mm
                        "non_coding_start": 130,
                        "non_coding_end": 138,
                        "mismatches": 0,
                        "gaps": 0,
                    },
                ],
            }
        }

        results = engine.filter_all_alignments(data)
        assert results['total_alignments'] == 2
        # Both should pass (8bp E~0.0076 < 0.01 threshold)
        assert results['passed'] == 2

    def test_legacy_pipeline_engine_runs(self):
        """Legacy pipeline should work end-to-end through FilterEngine."""
        pipeline = create_legacy_pipeline()
        engine = FilterEngine(pipeline)

        data = {
            "test_transposon": {
                "start": 100,
                "end": 200,
                "alignments": [
                    {
                        "aligned_sequence": "ACGTACGTACGTACGT",
                        "non_coding_start": 110,
                        "non_coding_end": 126,
                    },
                ],
            }
        }

        results = engine.filter_all_alignments(data)
        assert results['total_alignments'] == 1
        assert results['passed'] == 1

    def test_statistical_significance_registered_in_engine(self):
        """statistical_significance should be in FilterEngine.FILTER_FUNCTIONS."""
        assert 'statistical_significance' in FilterEngine.FILTER_FUNCTIONS


# ============================================
# End-to-end: two sequences -> alignment -> E-value
# ============================================

def _ungapped_hit_to_filter_alignment(hit, query_length, target_length):
    """Convert an ungapped alignment finder hit to filter-compatible dict."""
    return {
        'aligned_sequence': hit['seq1'],
        'mismatches': hit['mismatches'],
        'gaps': 0,
        'query_length': query_length,
        'target_length': target_length,
    }


def _gapped_hit_to_filter_alignment(hit, query_length, target_length):
    """Convert a gapped alignment finder hit to filter-compatible dict."""
    return {
        'aligned_sequence': hit['seq1_aligned'],
        'mismatches': hit['mismatches'],
        'gaps': hit['gaps'],
        'alignment_string': hit['alignment_string'],
        'query_length': query_length,
        'target_length': target_length,
    }


class TestEndToEnd:
    """End-to-end tests: raw sequences -> alignment finder -> E-value filter."""

    # --- Ungapped cases ---

    def test_perfect_16bp_match(self):
        """16bp perfect match between flanking and non-coding -> high confidence."""
        flanking = 'ACGTACGTACGTACGT'                    # 16bp
        noncoding = 'NNNNNACGTACGTACGTACGTNNNNN'         # 26bp, contains the 16bp

        hits = find_alignments_between_sequences(flanking, noncoding,
                                                  min_length=9, max_mismatches=0)
        assert len(hits) >= 1

        best = hits[0]
        assert best['length'] == 16
        assert best['mismatches'] == 0

        alignment = _ungapped_hit_to_filter_alignment(best, len(flanking), len(noncoding))
        params = {'e_value_reject': 0.01, 'e_value_low': 1e-3, 'e_value_high': 1e-6}
        passed, reason, metrics = Filters.statistical_significance_filter(
            alignment['aligned_sequence'], params, alignment=alignment)

        # G = (16-16+1)*(26-16+1) = 1*11 = 11 -> E = 11/4^16 ~ 2.6e-9
        assert passed is True
        assert metrics['confidence'] == 'high'
        assert metrics['e_value'] < 1e-6
        assert metrics['search_space'] == 1 * 11

    def test_short_9bp_match_in_large_sequences(self):
        """9bp match in large sequences -> bigger search space -> higher E-value."""
        # Embed a 9bp motif in two longer sequences
        motif = 'GCTAGCTAG'  # 9bp
        flanking = 'A' * 50 + motif + 'T' * 41       # 100bp total
        noncoding = 'C' * 200 + motif + 'G' * 291    # 500bp total

        hits = find_alignments_between_sequences(flanking, noncoding,
                                                  min_length=9, max_mismatches=0)
        # Find the hit that matches our motif
        motif_hits = [h for h in hits if h['seq1'] == motif]
        assert len(motif_hits) >= 1

        best = motif_hits[0]
        alignment = _ungapped_hit_to_filter_alignment(best, len(flanking), len(noncoding))
        params = {'e_value_reject': 0.01, 'e_value_low': 1e-3, 'e_value_high': 1e-6}
        passed, reason, metrics = Filters.statistical_significance_filter(
            alignment['aligned_sequence'], params, alignment=alignment)

        # G = (100-9+1)*(500-9+1) = 92*492 = 45264
        assert metrics['search_space'] == 92 * 492
        # E = 45264 / 4^9 ~ 0.173 -> rejected
        assert passed is False
        assert metrics['confidence'] == 'rejected'

    def test_same_9bp_in_small_sequences_passes(self):
        """Same 9bp motif but in short sequences -> small search space -> passes."""
        motif = 'GCTAGCTAG'  # 9bp
        flanking = motif + 'AA'     # 11bp
        noncoding = 'CC' + motif    # 11bp

        hits = find_alignments_between_sequences(flanking, noncoding,
                                                  min_length=9, max_mismatches=0)
        motif_hits = [h for h in hits if h['seq1'] == motif]
        assert len(motif_hits) >= 1

        best = motif_hits[0]
        alignment = _ungapped_hit_to_filter_alignment(best, len(flanking), len(noncoding))
        params = {'e_value_reject': 0.01, 'e_value_low': 1e-3, 'e_value_high': 1e-6}
        passed, reason, metrics = Filters.statistical_significance_filter(
            alignment['aligned_sequence'], params, alignment=alignment)

        # G = (11-9+1)*(11-9+1) = 3*3 = 9
        assert metrics['search_space'] == 9
        # E = 9 / 4^9 ~ 3.4e-5 -> low confidence, passes
        assert passed is True
        assert metrics['confidence'] == 'low'

    def test_match_with_1_mismatch(self):
        """12bp match with 1 mismatch between flanking and non-coding."""
        flanking = 'TTTACGTACGTACTTT'              # 16bp, motif at [3:15]
        noncoding = 'NNNACGTACGCACNNN'              # 16bp, 1mm at pos 10 (T->C)

        hits = find_alignments_between_sequences(flanking, noncoding,
                                                  min_length=9, max_mismatches=1)
        # Should find a ~12bp hit with 1 mismatch
        long_hits = [h for h in hits if h['length'] >= 10 and h['orientation'] == 'forward']
        assert len(long_hits) >= 1

        best = max(long_hits, key=lambda h: h['length'])
        alignment = _ungapped_hit_to_filter_alignment(best, len(flanking), len(noncoding))
        params = {'e_value_reject': 0.01, 'e_value_low': 1e-3, 'e_value_high': 1e-6}
        passed, reason, metrics = Filters.statistical_significance_filter(
            alignment['aligned_sequence'], params, alignment=alignment)

        assert metrics['mismatches'] == 1
        assert metrics['alignment_length'] == best['length']
        # With 1mm the E-value is higher than 0mm of same length
        ev_0mm = calculate_ungapped_evalue(best['length'], 0, metrics['search_space'])
        assert metrics['e_value'] > ev_0mm

    # --- Gapped case ---

    def test_gapped_alignment_end_to_end(self):
        """Gapped alignment: flanking has insertion relative to non-coding."""
        #                   seed region (9bp perfect)
        flanking = 'ACGTACGTANNNACGTACGTA'          # 21bp, 3bp insertion in middle
        noncoding = 'CCCACGTACGTAACGTACGTACCC'      # 24bp, continuous

        # Find ungapped seed first
        seeds = find_alignments_between_sequences(flanking, noncoding,
                                                   min_length=9, max_mismatches=0)
        assert len(seeds) >= 1

        # Extend with gaps
        gapped = extend_alignment_with_gaps(flanking, noncoding, seeds[0],
                                             min_identity=0.70)
        if gapped is None:
            pytest.skip("Gapped extension didn't meet identity threshold for this seed")

        alignment = _gapped_hit_to_filter_alignment(gapped, len(flanking), len(noncoding))
        params = {
            'e_value_reject': 0.01, 'e_value_low': 1e-3, 'e_value_high': 1e-6,
            'match_score': 2, 'mismatch_penalty': -1,
            'gap_open_penalty': -3, 'gap_extend_penalty': -1,
            'karlin_K': 0.1, 'karlin_lambda': None,
        }
        passed, reason, metrics = Filters.statistical_significance_filter(
            alignment['aligned_sequence'], params, alignment=alignment)

        assert metrics['is_gapped'] is True
        assert 'raw_score' in metrics
        assert 'e_value' in metrics
        assert metrics['e_value'] >= 0

    # --- Full pipeline end-to-end ---

    def test_full_pipeline_from_sequences(self):
        """Full pipeline: sequences -> finder -> FilterEngine -> results."""
        flanking = 'GCTAGCTAGCTAGCTA'                # 16bp
        noncoding = 'AAAAGCTAGCTAGCTAGCTAAAAA'       # 24bp

        hits = find_alignments_between_sequences(flanking, noncoding,
                                                  min_length=9, max_mismatches=0)
        assert len(hits) >= 1
        best = hits[0]

        # Build the data structure FilterEngine expects
        # Use start/end far from the alignment so boundary filter doesn't trip
        data = {
            "test_transposon": {
                "start": -1000,
                "end": len(noncoding) + 1000,
                "alignments": [
                    {
                        "aligned_sequence": best['seq1'],
                        "mismatches": best['mismatches'],
                        "gaps": 0,
                        "query_length": len(flanking),
                        "target_length": len(noncoding),
                        "non_coding_start": best['pos2'],
                        "non_coding_end": best['pos2'] + best['length'],
                    }
                ],
            }
        }

        pipeline = create_default_pipeline()
        engine = FilterEngine(pipeline)
        results = engine.filter_all_alignments(data)

        assert results['total_alignments'] == 1
        assert results['passed'] == 1

        detail = results['details'][0]['filter_result']
        sig_metrics = detail['filter_details']['statistical_significance']['metrics']
        assert 'e_value' in sig_metrics
        assert sig_metrics['search_space'] == (len(flanking) - best['length'] + 1) * (len(noncoding) - best['length'] + 1)
        print(f"\n  Alignment: {best['length']}bp, {best['mismatches']}mm")
        print(f"  Search space: {sig_metrics['search_space']}")
        print(f"  E-value: {sig_metrics['e_value']:.2e}")
        print(f"  Confidence: {sig_metrics['confidence']}")
