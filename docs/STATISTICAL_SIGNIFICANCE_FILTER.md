# Statistical Significance Filter

E-value based filtering for RNA guide alignments. Replaces the old length-based
hard cutoff + confidence zone approach with a probabilistic model that accounts
for alignment length, mismatches, gaps, and search space size.

## Why E-values instead of length cutoffs

The old filter used arbitrary thresholds: <9bp rejected, 9-13bp low, 14+ high.
This can't distinguish between:

| Alignment | Old filter | E-value (G=500) | Actual significance |
|-----------|-----------|-----------------|---------------------|
| 9bp, 0 mismatches | low | 1.91e-03 | Likely real |
| 9bp, 1 mismatch | low | 5.15e-02 | Probably noise |

The E-value answers: "How unlikely is it that a random sequence would produce
this alignment?" A 27x difference in significance that the old filter missed.

## Full pipeline: flanking region to E-value

```
STAGE 1                    STAGE 2                   STAGE 3              STAGE 4
Region Extraction    ->    Alignment Finding    ->    Field Mapping   ->   Filter Pipeline
(region_extractor)         (short_alignment_finder)   (glue code)          (alignment_filter)
```

### Stage 1: Region Extraction

`RegionExtractor` extracts flanking and coding regions from a genome:

```python
from modules.region_extractor import RegionExtractor

extractor = RegionExtractor("genome.fna")
regions = extractor.extract_protein_regions(
    contig_id="CP000699.1",
    start=3348491,          # transposase gene start (1-based)
    end=3350016,            # transposase gene end
    strand=1,
    upstream_length=350,
    downstream_length=250
)

flanking = regions['upstream_sequence']    # up to 350bp
noncoding = ...                            # internal non-coding region of transposon
```

**Output:**
- `upstream_sequence` - flanking upstream of the transposase gene
- `downstream_sequence` - flanking downstream
- `coding_sequence` - the transposase gene itself

The non-coding region is the internal part of the transposon between coding
genes, extracted separately from annotation data.

### Stage 2: Alignment Finding

`find_alignments_between_sequences` searches for short matches between the
flanking region and the non-coding region:

```python
from modules.short_alignment_finder import (
    find_alignments_between_sequences,
    extend_alignment_with_gaps,
)

# Ungapped alignment (exact or near-exact matches)
hits = find_alignments_between_sequences(
    flanking,              # sequence A
    noncoding,             # sequence B
    min_length=9,
    max_mismatches=1
)
```

**Algorithm:**
1. For each position `i` in flanking, `j` in noncoding:
   compare `flanking[i:i+9]` vs `noncoding[j:j+9]`
2. If mismatches <= limit: seed match found
3. Extend seed bidirectionally while staying within mismatch limit
4. Consolidate overlapping hits, keep longest

**Ungapped hit output:**
```python
{
    'pos1': 24,                    # position in flanking
    'pos2': 16,                    # position in noncoding
    'length': 15,                  # alignment length
    'seq1': 'TTTAAACCCGGGAAA',     # matched sequence from flanking
    'seq2': 'TTTAAACCCGGGAAA',     # matched sequence from noncoding
    'mismatches': 0,
    'mismatch_positions': [],
    'orientation': 'forward'       # or 'reverse_complement'
}
```

**Optional: gapped extension** for alignments with insertions/deletions:

```python
gapped = extend_alignment_with_gaps(flanking, noncoding, hits[0],
                                     min_identity=0.80)
```

**Gapped hit output** (additional fields):
```python
{
    ...
    'seq1_aligned': 'ACGT---ACGT',     # with '-' for gaps
    'seq2_aligned': 'ACGTACGACGT',
    'alignment_length': 14,             # total including gaps
    'matches': 11,
    'mismatches': 0,
    'gaps': 3,
    'identity': 0.786,
    'alignment_string': '|||| |||||||',  # '|'=match, '.'=mismatch, ' '=gap
}
```

### Stage 3: Field Mapping

The alignment finder and filter use different field names. This mapping
connects them:

| Alignment finder output | Filter input field | Notes |
|---|---|---|
| `seq1` (ungapped) | `aligned_sequence` | The matched sequence text |
| `seq1_aligned` (gapped) | `aligned_sequence` | With '-' gap characters |
| `mismatches` | `mismatches` | Same field name |
| *(not in output)* | `gaps` | Set to 0 for ungapped hits |
| *(not in output)* | `alignment_string` | Only needed for gapped hits |
| *(not in output)* | `query_length` | `len(flanking)` - full sequence A length |
| *(not in output)* | `target_length` | `len(noncoding)` - full sequence B length |

**Ungapped hit conversion:**

```python
def ungapped_hit_to_filter_alignment(hit, query_length, target_length):
    return {
        'aligned_sequence': hit['seq1'],
        'mismatches': hit['mismatches'],
        'gaps': 0,
        'query_length': query_length,
        'target_length': target_length,
    }
```

**Gapped hit conversion:**

```python
def gapped_hit_to_filter_alignment(hit, query_length, target_length):
    return {
        'aligned_sequence': hit['seq1_aligned'],
        'mismatches': hit['mismatches'],
        'gaps': hit['gaps'],
        'alignment_string': hit['alignment_string'],
        'query_length': query_length,
        'target_length': target_length,
    }
```

### Stage 4: Statistical Significance Filter

The filter computes an E-value and assigns confidence.

#### Search space computation

When `query_length` and `target_length` are provided in the alignment dict:

```
G = (query_length - n + 1) x (target_length - n + 1)
```

where `n` = alignment length. This is the number of positions in both sequences
where this alignment could start by chance.

Example: flanking=350bp, noncoding=200bp, alignment=15bp:
```
G = (350 - 15 + 1) x (200 - 15 + 1) = 336 x 186 = 62,496
```

Falls back to `params['search_space']` (default 500) if sequence lengths
are not provided.

#### Ungapped E-value (gaps == 0)

Exact combinatorial probability computed in log-space:

```
E = G x C(n,k) x 3^k / 4^n
```

where `n` = alignment length, `k` = mismatches, `G` = search space.

- `1/4^n` = probability of a specific n-bp sequence occurring by chance
- `C(n,k) x 3^k` = number of ways to place k mismatches (each has 3 wrong bases)
- `G` = number of positions to test

#### Gapped E-value (gaps > 0)

Uses Karlin-Altschul statistics:

```
E = K x m x n x exp(-lambda x S)
```

where:
- `S` = raw alignment score = `matches*2 + mismatches*(-1) + gap_opens*(-3) + gap_extensions*(-1)`
- `lambda` = solved by Newton-Raphson from `0.25*exp(2*lam) + 0.75*exp(-lam) = 1`
  (for match=2, mismatch=-1: lambda ~ 0.2645)
- `K` = 0.1 (default Karlin-Altschul constant)
- `m`, `n` = full query and target sequence lengths

Gap opens are counted from the `alignment_string` field (transitions from
non-gap to gap character `' '`).

#### Confidence assignment

| E-value range | Confidence | pass |
|---|---|---|
| E > 0.01 | rejected | False |
| 1e-3 < E <= 0.01 | low | True |
| 1e-6 < E <= 1e-3 | low | True |
| E <= 1e-6 | high | True |

#### Filter output

```python
passed, reason, metrics = Filters.statistical_significance_filter(
    sequence, params, alignment=alignment
)

# metrics contains:
{
    'e_value': 5.82e-05,
    'p_value': 5.82e-05,
    'search_space': 62496,
    'alignment_length': 15,
    'mismatches': 0,
    'is_gapped': False,
    'confidence': 'low',
}
```

## Complete end-to-end example

```python
from modules.short_alignment_finder import find_alignments_between_sequences
from modules.alignment_filter import (
    FilterEngine, create_default_pipeline, Filters
)

# --- Input sequences ---
flanking = "ATGCTAGCTAGCTTTAAACCCGGGAAAGCTAGCTAGCTA"   # 39bp
noncoding = "NNNNNNNNTTTAAACCCGGGAAANNNNNNNNN"         # 32bp

# --- Stage 2: Find alignments ---
hits = find_alignments_between_sequences(
    flanking, noncoding, min_length=9, max_mismatches=1
)

# --- Stage 3: Map fields ---
best = hits[0]
alignment = {
    'aligned_sequence': best['seq1'],
    'mismatches': best['mismatches'],
    'gaps': 0,
    'query_length': len(flanking),      # 39
    'target_length': len(noncoding),    # 32
    'non_coding_start': best['pos2'] + 100,   # offset for boundary filter
    'non_coding_end': best['pos2'] + best['length'] + 100,
}

# --- Stage 4a: Direct filter call ---
params = {
    'e_value_reject': 0.01,
    'e_value_low': 1e-3,
    'e_value_high': 1e-6,
}
passed, reason, metrics = Filters.statistical_significance_filter(
    alignment['aligned_sequence'], params, alignment=alignment
)
print(f"E-value: {metrics['e_value']:.2e}")
print(f"Search space: {metrics['search_space']}")
print(f"Confidence: {metrics['confidence']}")

# --- Stage 4b: Or use full pipeline ---
data = {
    "transposon_1": {
        "start": -1000,
        "end": 2000,
        "alignments": [alignment],
    }
}
pipeline = create_default_pipeline()
engine = FilterEngine(pipeline)
results = engine.filter_all_alignments(data)
```

## Pipeline presets

| Pipeline | E reject | E low | E high | Use case |
|---|---|---|---|---|
| `create_default_pipeline()` | 0.01 | 1e-3 | 1e-6 | General analysis |
| `create_strict_pipeline()` | 1e-3 | 1e-6 | 1e-9 | High-confidence only |
| `create_relaxed_pipeline()` | 0.05 | 0.01 | 1e-3 | Exploratory / short sequences |
| `create_legacy_pipeline()` | *(n/a)* | *(n/a)* | *(n/a)* | Old length-based filters |

The legacy pipeline uses the original `length_hard_cutoff` + `length_confidence`
filters for backwards compatibility.

## Reference E-values

Ungapped, default search space G=500:

| Alignment | E-value | Confidence |
|-----------|---------|------------|
| 8bp, 0mm | 7.63e-03 | low |
| 9bp, 0mm | 1.91e-03 | low |
| 9bp, 1mm | 5.15e-02 | **rejected** |
| 12bp, 0mm | 2.98e-05 | low |
| 15bp, 0mm | 4.66e-07 | **high** |
| 17bp, 1mm | 1.48e-06 | low |
| 20bp, 2mm | 7.78e-07 | high |

The same alignment can be significant or not depending on search space:

| Alignment | Flanking | Non-coding | Search space G | E-value | Confidence |
|-----------|----------|------------|----------------|---------|------------|
| 9bp, 0mm | 11bp | 11bp | 9 | 3.4e-05 | low |
| 9bp, 0mm | 100bp | 500bp | 45,264 | 0.173 | **rejected** |

## Functions added to alignment_filter.py

**Calculation functions** (in "Basic Calculation Functions" section):

- `calculate_ungapped_evalue(alignment_length, mismatches, search_space)` - exact combinatorial E-value
- `calculate_gapped_evalue(score, query_len, target_len, K, lam)` - Karlin-Altschul E-value
- `compute_lambda_for_scoring(match_score, mismatch_penalty)` - Newton-Raphson lambda solver
- `evalue_to_pvalue(evalue)` - E-value to P-value conversion
- `count_gap_opens(alignment_string)` - count gap-open events in alignment string

**Filter function** (in `Filters` class):

- `Filters.statistical_significance_filter(sequence, params, alignment, transposon_data)` - the main filter

**Filter params dict:**

```python
{
    'e_value_reject': 0.01,       # E > this -> rejected
    'e_value_low': 1e-3,          # E > this -> low confidence
    'e_value_high': 1e-6,         # E <= this -> high confidence
    'match_score': 2,             # scoring for gapped path
    'mismatch_penalty': -1,
    'gap_open_penalty': -3,
    'gap_extend_penalty': -1,
    'karlin_K': 0.1,
    'karlin_lambda': None,        # auto-computed if None
    'search_space': 500,          # fallback if query/target_length not in alignment
}
```

## Tests

Run with:
```
python -m pytest tests/test_alignment_filter.py -v
```

Test classes:
- `TestCalculateUngappedEvalue` - E-value formula correctness against known values
- `TestComputeLambda` - Newton-Raphson solver convergence
- `TestEvalueToPvalue` - conversion edge cases
- `TestCountGapOpens` - gap counting logic
- `TestCalculateGappedEvalue` - Karlin-Altschul formula
- `TestStatisticalSignificanceFilter` - filter confidence assignment, search space computation
- `TestPipelines` - preset pipeline configuration
- `TestEndToEnd` - full flow from raw sequences through alignment finder to E-value
