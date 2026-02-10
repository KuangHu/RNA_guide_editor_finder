"""
Alignment Filter Module for RNA Guide Editor Finder

This module provides a flexible, configurable filtering pipeline for
RNA guide alignments with the following capabilities:

1. Composition filters: AT/GC content, entropy, homopolymers, dinucleotide repeats
2. Length filters: Hard cutoffs and confidence-based length filtering
3. Topology filters: Boundary artifact detection

Usage:
    from modules.alignment_filter import FilterEngine, create_default_pipeline

    pipeline = create_default_pipeline()
    engine = FilterEngine(pipeline)
    results = engine.filter_all_alignments(data)
"""

import math
from collections import Counter
from typing import Dict, List, Tuple
from dataclasses import dataclass, field
from enum import Enum


# ============================================
# Configuration Classes
# ============================================

class FilterLevel(Enum):
    """Filter level categories"""
    COMPOSITION = 1  # Sequence composition
    LENGTH = 2       # Length and probability
    TOPOLOGY = 3     # Topological position


@dataclass
class FilterConfig:
    """Configuration for a single filter"""
    name: str
    enabled: bool = True
    level: FilterLevel = FilterLevel.COMPOSITION
    params: Dict = field(default_factory=dict)
    description: str = ""


@dataclass
class FilterPipeline:
    """Filter pipeline configuration"""
    filters: List[FilterConfig] = field(default_factory=list)

    def add_filter(self, name: str, level: FilterLevel, enabled: bool = True,
                   params: Dict = None, description: str = ""):
        """Add a filter to the pipeline"""
        self.filters.append(FilterConfig(
            name=name,
            enabled=enabled,
            level=level,
            params=params or {},
            description=description
        ))

    def enable_filter(self, name: str):
        """Enable a specific filter"""
        for f in self.filters:
            if f.name == name:
                f.enabled = True
                break

    def disable_filter(self, name: str):
        """Disable a specific filter"""
        for f in self.filters:
            if f.name == name:
                f.enabled = False
                break

    def update_params(self, name: str, params: Dict):
        """Update filter parameters"""
        for f in self.filters:
            if f.name == name:
                f.params.update(params)
                break

    def get_enabled_filters(self) -> List[FilterConfig]:
        """Get all enabled filters"""
        return [f for f in self.filters if f.enabled]


# ============================================
# Basic Calculation Functions
# ============================================

def calculate_at_content(sequence: str) -> float:
    """Calculate AT content percentage of a sequence"""
    sequence = sequence.upper()
    at_count = sequence.count('A') + sequence.count('T')
    return (at_count / len(sequence)) * 100 if len(sequence) > 0 else 0


def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content percentage of a sequence"""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0


def calculate_shannon_entropy(sequence: str) -> float:
    """Calculate Shannon entropy (complexity) of a sequence"""
    sequence = sequence.upper()
    if len(sequence) == 0:
        return 0

    freq = Counter(sequence)
    entropy = 0
    for count in freq.values():
        p = count / len(sequence)
        if p > 0:
            entropy -= p * math.log2(p)

    return entropy


def has_homopolymer(sequence: str, min_length: int = 5) -> bool:
    """Check if sequence has a homopolymer run"""
    sequence = sequence.upper()
    for base in ['A', 'T', 'G', 'C']:
        if base * min_length in sequence:
            return True
    return False


def get_longest_homopolymer(sequence: str) -> int:
    """Get the length of the longest homopolymer run"""
    sequence = sequence.upper()
    max_length = 0
    for base in ['A', 'T', 'G', 'C']:
        current = 0
        for char in sequence:
            if char == base:
                current += 1
                max_length = max(max_length, current)
            else:
                current = 0
    return max_length


def has_dinucleotide_repeat(sequence: str, threshold: float = 0.8) -> bool:
    """Check if sequence has excessive dinucleotide repeats"""
    sequence = sequence.upper()
    if len(sequence) < 6:
        return False

    dinuc_patterns = ['AT', 'TA', 'GT', 'TG', 'AC', 'CA', 'AG', 'GA', 'CT', 'TC', 'GC', 'CG']

    for pattern in dinuc_patterns:
        matches = sum(1 for i in range(len(sequence) - 1)
                     if sequence[i:i+2] == pattern)

        if matches / (len(sequence) - 1) >= threshold:
            return True

    return False


def calculate_dinucleotide_ratio(sequence: str) -> float:
    """Calculate the maximum dinucleotide repeat ratio"""
    sequence = sequence.upper()
    if len(sequence) < 2:
        return 0.0

    dinuc_patterns = ['AT', 'TA', 'GT', 'TG', 'AC', 'CA', 'AG', 'GA', 'CT', 'TC', 'GC', 'CG']
    max_ratio = 0.0

    for pattern in dinuc_patterns:
        matches = sum(1 for i in range(len(sequence) - 1)
                     if sequence[i:i+2] == pattern)
        ratio = matches / (len(sequence) - 1)
        max_ratio = max(max_ratio, ratio)

    return max_ratio


def calculate_ungapped_evalue(alignment_length: int, mismatches: int,
                              search_space: int) -> float:
    """Calculate E-value for an ungapped alignment using exact combinatorial probability.

    Formula (computed in log-space to avoid overflow):
        E = G * C(n,k) * 3^k / 4^n

    where G = search_space, n = alignment_length, k = mismatches.
    """
    n = alignment_length
    k = mismatches
    if k > n or k < 0 or n <= 0:
        return float('inf')

    # log(C(n,k)) + k*log(3) - n*log(4) + log(G)
    log_comb = math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    log_e = math.log(search_space) + log_comb + k * math.log(3) - n * math.log(4)
    return math.exp(log_e)


def calculate_gapped_evalue(score: float, query_len: int, target_len: int,
                            K: float, lam: float) -> float:
    """Calculate E-value for a gapped alignment using Karlin-Altschul statistics.

    Formula: E = K * m * n * exp(-lambda * S)
    """
    if score <= 0:
        return float('inf')
    log_e = math.log(K) + math.log(query_len) + math.log(target_len) - lam * score
    return math.exp(log_e)


def compute_lambda_for_scoring(match_score: int = 2, mismatch_penalty: int = -1,
                               max_iter: int = 100, tol: float = 1e-10) -> float:
    """Compute Karlin-Altschul lambda using Newton-Raphson.

    Solves: 0.25*exp(match*lam) + 0.75*exp(mismatch*lam) - 1 = 0
    assuming equal base frequencies (p=0.25 each), so match prob = 0.25,
    mismatch prob = 0.75.

    For match=2, mismatch=-1: lambda ~ ln(2) ~ 0.693
    """
    s = match_score
    r = mismatch_penalty
    lam = 0.5  # initial guess

    for _ in range(max_iter):
        f = 0.25 * math.exp(s * lam) + 0.75 * math.exp(r * lam) - 1.0
        fp = 0.25 * s * math.exp(s * lam) + 0.75 * r * math.exp(r * lam)
        if abs(fp) < 1e-15:
            break
        lam_new = lam - f / fp
        if lam_new <= 0:
            lam_new = lam / 2.0
        if abs(lam_new - lam) < tol:
            lam = lam_new
            break
        lam = lam_new

    return lam


def evalue_to_pvalue(evalue: float) -> float:
    """Convert E-value to P-value: P = 1 - exp(-E)"""
    if evalue >= 700:
        return 1.0
    return 1.0 - math.exp(-evalue)


def count_gap_opens(alignment_string: str) -> int:
    """Count gap-open events (non-gap -> gap transitions) in an alignment string.

    Gaps are represented as spaces ' ' in the alignment string.
    """
    if not alignment_string:
        return 0
    opens = 0
    in_gap = False
    for ch in alignment_string:
        if ch == ' ':
            if not in_gap:
                opens += 1
                in_gap = True
        else:
            in_gap = False
    return opens


def is_boundary_artifact(alignment: Dict, transposon_start: int, transposon_end: int,
                        boundary_distance: int = 10) -> bool:
    """Check if alignment is a boundary annotation artifact"""
    non_coding_start = alignment['non_coding_start']
    non_coding_end = alignment['non_coding_end']

    distance_to_start = abs(non_coding_start - transposon_start)
    distance_to_end = abs(non_coding_end - transposon_end)

    return min(distance_to_start, distance_to_end) < boundary_distance


# ============================================
# Modular Filter Functions
# ============================================

class Filters:
    """Collection of all available filter functions"""

    @staticmethod
    def at_content_filter(sequence: str, params: Dict) -> Tuple[bool, str, Dict]:
        """AT content filter"""
        threshold = params.get('max_at_percent', 75.0)
        at_content = calculate_at_content(sequence)

        passed = at_content <= threshold
        reason = "" if passed else f"AT_CONTENT_TOO_HIGH: {at_content:.1f}% > {threshold}%"
        metrics = {'at_content': round(at_content, 2)}

        return passed, reason, metrics

    @staticmethod
    def gc_content_filter(sequence: str, params: Dict) -> Tuple[bool, str, Dict]:
        """GC content filter"""
        min_gc = params.get('min_gc_percent', 0.0)
        max_gc = params.get('max_gc_percent', 100.0)
        gc_content = calculate_gc_content(sequence)

        passed = min_gc <= gc_content <= max_gc
        reason = "" if passed else f"GC_CONTENT_OUT_OF_RANGE: {gc_content:.1f}% not in [{min_gc}, {max_gc}]"
        metrics = {'gc_content': round(gc_content, 2)}

        return passed, reason, metrics

    @staticmethod
    def entropy_filter(sequence: str, params: Dict) -> Tuple[bool, str, Dict]:
        """Sequence complexity (entropy) filter"""
        min_entropy = params.get('min_entropy', 1.5)
        entropy = calculate_shannon_entropy(sequence)

        passed = entropy >= min_entropy
        reason = "" if passed else f"LOW_COMPLEXITY: entropy={entropy:.2f} < {min_entropy}"
        metrics = {'entropy': round(entropy, 2)}

        return passed, reason, metrics

    @staticmethod
    def homopolymer_filter(sequence: str, params: Dict) -> Tuple[bool, str, Dict]:
        """Homopolymer filter"""
        max_homopolymer = params.get('max_homopolymer_length', 5)
        longest = get_longest_homopolymer(sequence)

        passed = longest <= max_homopolymer
        reason = "" if passed else f"HOMOPOLYMER_TOO_LONG: {longest}bp > {max_homopolymer}bp"
        metrics = {'longest_homopolymer': longest}

        return passed, reason, metrics

    @staticmethod
    def dinucleotide_filter(sequence: str, params: Dict) -> Tuple[bool, str, Dict]:
        """Dinucleotide repeat filter"""
        max_ratio = params.get('max_dinucleotide_ratio', 0.8)
        ratio = calculate_dinucleotide_ratio(sequence)

        passed = ratio <= max_ratio
        reason = "" if passed else f"DINUCLEOTIDE_REPEAT: {ratio:.1%} > {max_ratio:.1%}"
        metrics = {'dinucleotide_ratio': round(ratio, 3)}

        return passed, reason, metrics

    @staticmethod
    def length_hard_cutoff_filter(sequence: str, params: Dict) -> Tuple[bool, str, Dict]:
        """Length hard cutoff filter"""
        min_length = params.get('min_length', 9)
        length = len(sequence)

        passed = length >= min_length
        reason = "" if passed else f"TOO_SHORT: {length}bp < {min_length}bp"
        metrics = {'length': length}

        return passed, reason, metrics

    @staticmethod
    def length_confidence_filter(sequence: str, params: Dict) -> Tuple[bool, str, Dict]:
        """Length confidence filter"""
        low_conf_min = params.get('low_confidence_min', 9)
        high_conf_min = params.get('high_confidence_min', 14)
        length = len(sequence)

        if length < low_conf_min:
            confidence = 'rejected'
            reason = f"LENGTH_REJECTED: {length}bp < {low_conf_min}bp"
        elif length < high_conf_min:
            confidence = 'low'
            reason = f"LOW_CONFIDENCE_LENGTH: {length}bp in doubt zone [{low_conf_min}, {high_conf_min})"
        else:
            confidence = 'high'
            reason = ""

        passed = confidence != 'rejected'
        metrics = {'length': length, 'confidence': confidence}

        return passed, reason, metrics

    @staticmethod
    def boundary_artifact_filter(sequence: str, params: Dict,
                                 alignment: Dict = None,
                                 transposon_data: Dict = None) -> Tuple[bool, str, Dict]:
        """Boundary artifact filter"""
        if alignment is None or transposon_data is None:
            return True, "", {}

        boundary_distance = params.get('boundary_distance', 10)
        is_artifact = is_boundary_artifact(
            alignment,
            transposon_data['start'],
            transposon_data['end'],
            boundary_distance
        )

        passed = not is_artifact
        reason = "" if passed else "BOUNDARY_ARTIFACT: too close to transposon boundary"

        non_coding_start = alignment['non_coding_start']
        transposon_start = transposon_data['start']
        transposon_end = transposon_data['end']

        metrics = {
            'distance_to_start': abs(non_coding_start - transposon_start),
            'distance_to_end': abs(alignment['non_coding_end'] - transposon_end)
        }

        return passed, reason, metrics

    @staticmethod
    def statistical_significance_filter(sequence: str, params: Dict,
                                        alignment: Dict = None,
                                        transposon_data: Dict = None) -> Tuple[bool, str, Dict]:
        """Statistical significance filter using E-value calculations.

        Determines alignment confidence based on how unlikely the match
        is to occur by chance, accounting for mismatches, gaps, and
        search space size.

        The search space is computed from the full sequence lengths:
            G = (query_length - n + 1) * (target_length - n + 1)
        where query_length = len(flanking region), target_length = len(non-coding region),
        and n = alignment_length. These lengths are read from the alignment dict
        ('query_length' and 'target_length' fields). Falls back to params['search_space']
        if not provided.
        """
        e_value_reject = params.get('e_value_reject', 0.01)
        e_value_low = params.get('e_value_low', 1e-3)
        e_value_high = params.get('e_value_high', 1e-6)
        match_score = params.get('match_score', 2)
        mismatch_penalty = params.get('mismatch_penalty', -1)
        gap_open_penalty = params.get('gap_open_penalty', -3)
        gap_extend_penalty = params.get('gap_extend_penalty', -1)
        karlin_K = params.get('karlin_K', 0.1)
        karlin_lambda = params.get('karlin_lambda', None)

        # Get alignment info
        mismatches = 0
        gaps = 0
        alignment_string = ''
        alignment_length = len(sequence)
        query_length = None
        target_length = None

        if alignment is not None:
            mismatches = alignment.get('mismatches', 0)
            gaps = alignment.get('gaps', 0)
            alignment_string = alignment.get('alignment_string', '')
            query_length = alignment.get('query_length', None)
            target_length = alignment.get('target_length', None)

        # Compute search space from full sequence lengths if available
        if query_length is not None and target_length is not None:
            n = alignment_length
            search_space = max(1, query_length - n + 1) * max(1, target_length - n + 1)
        else:
            # Fallback to explicit search_space param
            search_space = params.get('search_space', 500)
            if transposon_data is not None:
                search_space = transposon_data.get('search_space', search_space)

        is_gapped = gaps > 0

        if is_gapped:
            # Gapped path: compute raw score, then Karlin-Altschul E-value
            gap_opens = count_gap_opens(alignment_string)
            gap_extensions = gaps - gap_opens
            matches = alignment_length - mismatches - gaps
            raw_score = (matches * match_score +
                         mismatches * mismatch_penalty +
                         gap_opens * gap_open_penalty +
                         gap_extensions * gap_extend_penalty)

            if karlin_lambda is None:
                karlin_lambda = compute_lambda_for_scoring(match_score, mismatch_penalty)

            # For gapped Karlin-Altschul, use full sequence lengths directly
            ka_query = query_length if query_length is not None else alignment_length
            ka_target = target_length if target_length is not None else search_space

            evalue = calculate_gapped_evalue(raw_score, ka_query, ka_target,
                                             karlin_K, karlin_lambda)
            metrics = {
                'e_value': evalue,
                'p_value': evalue_to_pvalue(evalue),
                'raw_score': raw_score,
                'gap_opens': gap_opens,
                'search_space': search_space,
                'is_gapped': True,
            }
        else:
            # Ungapped path: exact combinatorial E-value
            evalue = calculate_ungapped_evalue(alignment_length, mismatches, search_space)
            metrics = {
                'e_value': evalue,
                'p_value': evalue_to_pvalue(evalue),
                'alignment_length': alignment_length,
                'mismatches': mismatches,
                'search_space': search_space,
                'is_gapped': False,
            }

        # Assign confidence based on E-value thresholds
        if evalue > e_value_reject:
            confidence = 'rejected'
            reason = f"E_VALUE_TOO_HIGH: E={evalue:.2e} > {e_value_reject:.2e}"
            passed = False
        elif evalue > e_value_low:
            confidence = 'low'
            reason = f"LOW_CONFIDENCE_EVALUE: E={evalue:.2e} in ({e_value_low:.2e}, {e_value_reject:.2e}]"
            passed = True
        elif evalue > e_value_high:
            confidence = 'low'
            reason = ""
            passed = True
        else:
            confidence = 'high'
            reason = ""
            passed = True

        metrics['confidence'] = confidence
        return passed, reason, metrics


# ============================================
# Main Filter Engine
# ============================================

class FilterEngine:
    """Flexible filtering engine"""

    # Map filter names to functions
    FILTER_FUNCTIONS = {
        'at_content': Filters.at_content_filter,
        'gc_content': Filters.gc_content_filter,
        'entropy': Filters.entropy_filter,
        'homopolymer': Filters.homopolymer_filter,
        'dinucleotide': Filters.dinucleotide_filter,
        'length_hard_cutoff': Filters.length_hard_cutoff_filter,
        'length_confidence': Filters.length_confidence_filter,
        'boundary_artifact': Filters.boundary_artifact_filter,
        'statistical_significance': Filters.statistical_significance_filter,
    }

    def __init__(self, pipeline: FilterPipeline):
        self.pipeline = pipeline

    def filter_alignment(self, alignment: Dict, transposon_data: Dict = None) -> Dict:
        """
        Apply filter pipeline to a single alignment.

        Returns:
        {
            'pass': bool,
            'confidence': str,
            'filter_reasons': List[str],
            'failed_filters': List[str],
            'passed_filters': List[str],
            'metrics': Dict,
            'filter_details': Dict
        }
        """
        sequence = alignment['aligned_sequence']

        result = {
            'pass': True,
            'confidence': 'high',
            'filter_reasons': [],
            'failed_filters': [],
            'passed_filters': [],
            'metrics': {},
            'filter_details': {}
        }

        # Group filters by level
        enabled_filters = self.pipeline.get_enabled_filters()
        filters_by_level = {level: [] for level in FilterLevel}
        for f in enabled_filters:
            filters_by_level[f.level].append(f)

        # Execute filters by level
        for level in [FilterLevel.COMPOSITION, FilterLevel.LENGTH, FilterLevel.TOPOLOGY]:
            for filter_config in filters_by_level[level]:
                filter_func = self.FILTER_FUNCTIONS.get(filter_config.name)

                if filter_func is None:
                    continue

                # Execute filter
                try:
                    # Some filters need additional context
                    if filter_config.name in ['boundary_artifact', 'statistical_significance']:
                        passed, reason, metrics = filter_func(
                            sequence,
                            filter_config.params,
                            alignment=alignment,
                            transposon_data=transposon_data
                        )
                    else:
                        passed, reason, metrics = filter_func(
                            sequence,
                            filter_config.params
                        )

                    # Record results
                    result['metrics'].update(metrics)
                    result['filter_details'][filter_config.name] = {
                        'passed': passed,
                        'reason': reason,
                        'metrics': metrics
                    }

                    if passed:
                        result['passed_filters'].append(filter_config.name)

                        # Special handling for confidence filters
                        if filter_config.name in ('length_confidence', 'statistical_significance') and 'confidence' in metrics:
                            result['confidence'] = metrics['confidence']
                    else:
                        result['pass'] = False
                        result['failed_filters'].append(filter_config.name)
                        if reason:
                            result['filter_reasons'].append(reason)

                except Exception as e:
                    result['filter_details'][filter_config.name] = {
                        'error': str(e)
                    }

        # Add basic sequence info
        result['metrics']['sequence'] = sequence
        result['metrics']['length'] = len(sequence)

        return result

    def filter_all_alignments(self, data: Dict) -> Dict:
        """
        Filter all alignments.

        Returns statistics and detailed results.
        """
        results = {
            'total_alignments': 0,
            'passed': 0,
            'failed': 0,
            'by_confidence': {'high': 0, 'low': 0, 'rejected': 0},
            'failed_by_filter': {},
            'details': []
        }

        # Initialize filter counts
        for filter_config in self.pipeline.filters:
            if filter_config.enabled:
                results['failed_by_filter'][filter_config.name] = 0

        for key, transposon_data in data.items():
            alignments = transposon_data.get('alignments', [])
            results['total_alignments'] += len(alignments)

            for alignment in alignments:
                filter_result = self.filter_alignment(alignment, transposon_data)

                # Statistics
                if filter_result['pass']:
                    results['passed'] += 1
                    confidence = filter_result.get('confidence', 'high')
                    results['by_confidence'][confidence] = results['by_confidence'].get(confidence, 0) + 1
                else:
                    results['failed'] += 1
                    results['by_confidence']['rejected'] += 1

                    # Count failures per filter
                    for failed_filter in filter_result['failed_filters']:
                        results['failed_by_filter'][failed_filter] += 1

                # Save detailed results
                results['details'].append({
                    'transposon_id': key,
                    'alignment': alignment,
                    'filter_result': filter_result
                })

        return results


# ============================================
# Preset Configurations
# ============================================

def create_default_pipeline() -> FilterPipeline:
    """Create default filter pipeline"""
    pipeline = FilterPipeline()

    # Level 1: Composition filters
    pipeline.add_filter(
        'at_content',
        FilterLevel.COMPOSITION,
        enabled=True,
        params={'max_at_percent': 75.0},
        description="Filter sequences with AT content > 75%"
    )

    pipeline.add_filter(
        'entropy',
        FilterLevel.COMPOSITION,
        enabled=True,
        params={'min_entropy': 1.5},
        description="Filter low-complexity sequences (Shannon entropy < 1.5)"
    )

    pipeline.add_filter(
        'homopolymer',
        FilterLevel.COMPOSITION,
        enabled=True,
        params={'max_homopolymer_length': 5},
        description="Filter sequences with homopolymers > 5bp"
    )

    pipeline.add_filter(
        'dinucleotide',
        FilterLevel.COMPOSITION,
        enabled=True,
        params={'max_dinucleotide_ratio': 0.8},
        description="Filter sequences with excessive dinucleotide repeats"
    )

    # Level 2: Statistical significance filter
    pipeline.add_filter(
        'statistical_significance',
        FilterLevel.LENGTH,
        enabled=True,
        params={
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
        },
        description="E-value based statistical significance filter"
    )

    # Level 3: Topology filters
    pipeline.add_filter(
        'boundary_artifact',
        FilterLevel.TOPOLOGY,
        enabled=True,
        params={'boundary_distance': 10},
        description="Filter boundary annotation artifacts"
    )

    return pipeline


def create_strict_pipeline() -> FilterPipeline:
    """Create strict filter pipeline"""
    pipeline = create_default_pipeline()

    # Stricter parameters
    pipeline.update_params('at_content', {'max_at_percent': 70.0})
    pipeline.update_params('entropy', {'min_entropy': 1.8})
    pipeline.update_params('homopolymer', {'max_homopolymer_length': 4})
    pipeline.update_params('statistical_significance', {
        'e_value_reject': 1e-3,
        'e_value_low': 1e-6,
        'e_value_high': 1e-9,
    })

    return pipeline


def create_relaxed_pipeline() -> FilterPipeline:
    """Create relaxed filter pipeline"""
    pipeline = create_default_pipeline()

    # More relaxed parameters
    pipeline.update_params('at_content', {'max_at_percent': 80.0})
    pipeline.update_params('entropy', {'min_entropy': 1.2})
    pipeline.update_params('homopolymer', {'max_homopolymer_length': 6})
    pipeline.update_params('statistical_significance', {
        'e_value_reject': 0.05,
        'e_value_low': 0.01,
        'e_value_high': 1e-3,
    })

    return pipeline


def create_legacy_pipeline() -> FilterPipeline:
    """Create legacy filter pipeline with original length-based filters.

    Provides backwards compatibility with the pre-statistical-significance
    filtering approach using hard length cutoffs and confidence zones.
    """
    pipeline = FilterPipeline()

    # Level 1: Composition filters
    pipeline.add_filter(
        'at_content',
        FilterLevel.COMPOSITION,
        enabled=True,
        params={'max_at_percent': 75.0},
        description="Filter sequences with AT content > 75%"
    )

    pipeline.add_filter(
        'entropy',
        FilterLevel.COMPOSITION,
        enabled=True,
        params={'min_entropy': 1.5},
        description="Filter low-complexity sequences (Shannon entropy < 1.5)"
    )

    pipeline.add_filter(
        'homopolymer',
        FilterLevel.COMPOSITION,
        enabled=True,
        params={'max_homopolymer_length': 5},
        description="Filter sequences with homopolymers > 5bp"
    )

    pipeline.add_filter(
        'dinucleotide',
        FilterLevel.COMPOSITION,
        enabled=True,
        params={'max_dinucleotide_ratio': 0.8},
        description="Filter sequences with excessive dinucleotide repeats"
    )

    # Level 2: Length filters (original)
    pipeline.add_filter(
        'length_hard_cutoff',
        FilterLevel.LENGTH,
        enabled=True,
        params={'min_length': 9},
        description="Hard cutoff: reject sequences < 9bp"
    )

    pipeline.add_filter(
        'length_confidence',
        FilterLevel.LENGTH,
        enabled=True,
        params={'low_confidence_min': 9, 'high_confidence_min': 14},
        description="Assign confidence levels based on length"
    )

    # Level 3: Topology filters
    pipeline.add_filter(
        'boundary_artifact',
        FilterLevel.TOPOLOGY,
        enabled=True,
        params={'boundary_distance': 10},
        description="Filter boundary annotation artifacts"
    )

    return pipeline


# ============================================
# Example Usage
# ============================================

if __name__ == "__main__":

    # Example 1: Using default configuration
    print("=" * 60)
    print("Example 1: Using default configuration")
    print("=" * 60)

    pipeline = create_default_pipeline()
    engine = FilterEngine(pipeline)

    sample_data = {
        "GCA_000016765.1_ASM1676v1_genomic.fna_CP000699.1_IS_13": {
            "start": 3348491,
            "end": 3350016,
            "alignments": [
                {
                    "sequence_source": "reverse",
                    "non_coding_region_index": 1,
                    "non_coding_start": 3349940,
                    "non_coding_end": 3350016,
                    "sequence_position": 10,
                    "non_coding_position": 0,
                    "length": 9,
                    "aligned_sequence": "CGGATCATC",
                    "alignment_type": "forward"
                }
            ]
        }
    }

    results = engine.filter_all_alignments(sample_data)
    print(f"Total: {results['total_alignments']}")
    print(f"Passed: {results['passed']}")
    print(f"Failed: {results['failed']}")
    print(f"By confidence: {results['by_confidence']}")
    print(f"Failed by filter: {results['failed_by_filter']}")

    # Example 2: Statistical significance E-values
    print("\n" + "=" * 60)
    print("Example 2: Statistical significance - E-value comparison")
    print("=" * 60)

    test_cases = [
        ("8bp, 0mm", 8, 0),
        ("9bp, 0mm", 9, 0),
        ("9bp, 1mm", 9, 1),
        ("12bp, 0mm", 12, 0),
        ("15bp, 0mm", 15, 0),
        ("17bp, 1mm", 17, 1),
        ("20bp, 2mm", 20, 2),
    ]

    print(f"{'Alignment':<12} {'E-value (G=500)':<20} {'Confidence'}")
    print("-" * 50)
    for label, length, mm in test_cases:
        ev = calculate_ungapped_evalue(length, mm, 500)
        # Determine confidence using default thresholds
        if ev > 0.01:
            conf = "rejected"
        elif ev > 1e-6:
            conf = "low"
        else:
            conf = "high"
        print(f"{label:<12} {ev:<20.2e} {conf}")

    # Example 3: Dynamic parameter adjustment
    print("\n" + "=" * 60)
    print("Example 3: Dynamic parameter adjustment")
    print("=" * 60)

    adjustable_pipeline = create_default_pipeline()

    # Disable some filters
    adjustable_pipeline.disable_filter('dinucleotide')
    adjustable_pipeline.disable_filter('boundary_artifact')

    # Adjust E-value thresholds
    adjustable_pipeline.update_params('at_content', {'max_at_percent': 85.0})
    adjustable_pipeline.update_params('statistical_significance', {'e_value_reject': 0.05})

    adjustable_engine = FilterEngine(adjustable_pipeline)
    adjustable_results = adjustable_engine.filter_all_alignments(sample_data)

    print(f"Enabled filters: {[f.name for f in adjustable_pipeline.get_enabled_filters()]}")
    print(f"Total: {adjustable_results['total_alignments']}")
    print(f"Passed: {adjustable_results['passed']}")
    print(f"Failed: {adjustable_results['failed']}")

    # Example 4: Detailed filter results
    print("\n" + "=" * 60)
    print("Example 4: Detailed filter results")
    print("=" * 60)

    for detail in results['details']:
        print(f"\nTransposon: {detail['transposon_id']}")
        print(f"Sequence: {detail['alignment']['aligned_sequence']}")
        print(f"Pass: {detail['filter_result']['pass']}")
        print(f"Confidence: {detail['filter_result']['confidence']}")
        print(f"Failed filters: {detail['filter_result']['failed_filters']}")
        print(f"Metrics: {detail['filter_result']['metrics']}")
        print("Filter details:")
        for fname, fdetail in detail['filter_result']['filter_details'].items():
            print(f"  {fname}: {fdetail}")
