# Finding Alignments Between Two Sequences

## Overview

The **Between-Sequences Alignment Finder** extends the Short Alignment Finder to identify alignments between **two different DNA sequences**, rather than within a single sequence.

This is particularly useful for:
- **Transposon analysis**: Finding shared sequences between flanking regions and non-coding regions
- **Regulatory element discovery**: Identifying motifs shared between upstream and downstream regions
- **Comparative genomics**: Finding conserved sequences between different genomic regions
- **Guide RNA identification**: Detecting guide RNA pairs in IS110 transposons across different regions

## Key Differences from Single-Sequence Search

| Feature | Single Sequence | Between Sequences |
|---------|----------------|-------------------|
| Input | One sequence | Two sequences |
| Distance constraints | Yes (min_distance, max_distance) | No (not applicable) |
| Position reference | pos1 and pos2 in same sequence | pos1 in seq1, pos2 in seq2 |
| Self-matching | Need to avoid | Not an issue |
| Use case | Find repeats within a region | Find shared motifs between regions |

## Quick Start

### Basic Example

```python
from modules.short_alignment_finder import find_alignments_between_sequences

# Two different sequences
flanking_region = "ATGCATGC" + "TTTAAACCCGGG" + "GCTAGCTA"
noncoding_region = "NNNNNNNN" + "TTTAAACCCGGG" + "NNNNNNNN"

# Find alignments between them
results = find_alignments_between_sequences(
    flanking_region,
    noncoding_region,
    min_length=9
)

# Output:
# [{
#     'pos1': 8,          # Position in flanking_region
#     'pos2': 8,          # Position in noncoding_region
#     'length': 12,
#     'seq1': 'TTTAAACCCGGG',
#     'seq2': 'TTTAAACCCGGG',
#     'mismatches': 0,
#     'mismatch_positions': [],
#     'orientation': 'forward'
# }]
```

### IS110 Transposon Use Case

```python
from modules.short_alignment_finder import ShortAlignmentFinder

# Flanking region (e.g., 100bp upstream of transposon)
flanking = extract_flanking_region(transposon_id)

# Non-coding region (internal region of transposon)
noncoding = extract_noncoding_region(transposon_id)

# Find guide RNA candidates
finder = ShortAlignmentFinder(
    min_length=9,          # IS110 guide RNAs typically 9-15bp
    max_mismatches=1,      # Allow biological variation
    check_forward=True,    # Direct matches
    check_revcomp=True     # Inverted repeats
)

results = finder.find_alignments_between(flanking, noncoding)

# Filter for high-quality candidates
guide_rna_candidates = [
    r for r in results
    if r['length'] >= 12 and r['mismatches'] == 0
]
```

## API Reference

### Method: `find_alignments_between()`

```python
finder = ShortAlignmentFinder(min_length=9, max_mismatches=1)
results = finder.find_alignments_between(sequence1, sequence2)
```

**Parameters:**
- `sequence1` (str): First DNA sequence (e.g., flanking region)
- `sequence2` (str): Second DNA sequence (e.g., non-coding region)

**Returns:**
- List of alignment dictionaries (see Output Format below)

**Note:** The `min_distance` and `max_distance` parameters set in the constructor are **ignored** for between-sequence searches (they don't make sense when comparing different sequences).

### Convenience Function

```python
from modules.short_alignment_finder import find_alignments_between_sequences

results = find_alignments_between_sequences(
    sequence1,
    sequence2,
    min_length=9,
    max_mismatches=0,
    check_forward=True,
    check_revcomp=True
)
```

**Parameters:**
- `sequence1` (str): First DNA sequence
- `sequence2` (str): Second DNA sequence
- `min_length` (int): Minimum alignment length (default: 9bp)
- `max_mismatches` (int): Maximum mismatches allowed (default: 0)
- `check_forward` (bool): Find direct matches (default: True)
- `check_revcomp` (bool): Find reverse complement matches (default: True)

## Output Format

Each alignment is returned as a dictionary:

```python
{
    'pos1': 12,                     # Position in sequence1 (0-based)
    'pos2': 8,                      # Position in sequence2 (0-based)
    'length': 15,                   # Length of alignment (bp)
    'seq1': 'TTTAAACCCGGGAAA',     # Sequence from sequence1
    'seq2': 'TTTAAACCCGGGAAA',     # Sequence from sequence2 (as it appears)
    'mismatches': 0,                # Number of mismatches
    'mismatch_positions': [],       # List of mismatch positions
    'orientation': 'forward'        # 'forward' or 'reverse_complement'
}
```

### Key Differences from Single-Sequence Output

- **No `distance` field**: Distance doesn't make sense between different sequences
- **pos1 and pos2**: Refer to positions in different sequences (not the same sequence)
- **Interpretation**: pos1 is in sequence1, pos2 is in sequence2

## Examples

### Example 1: Transposon Flanking vs Non-coding

```python
from modules.short_alignment_finder import ShortAlignmentFinder

# Extract regions from your transposon data
flanking = get_flanking_region(transposon_id, upstream_bp=100)
noncoding = get_noncoding_region(transposon_id)

# Find shared sequences
finder = ShortAlignmentFinder(min_length=9, max_mismatches=1)
alignments = finder.find_alignments_between(flanking, noncoding)

# Analyze results
for alignment in alignments:
    if alignment['length'] >= 12 and alignment['mismatches'] == 0:
        print(f"High-quality guide RNA candidate:")
        print(f"  Flanking position: {alignment['pos1']}")
        print(f"  Noncoding position: {alignment['pos2']}")
        print(f"  Sequence: {alignment['seq1']}")
```

### Example 2: Finding Inverted Repeats Between Regions

```python
# Look for palindromic sequences shared between regions
finder = ShortAlignmentFinder(
    min_length=8,
    max_mismatches=0,
    check_forward=False,      # Only inverted repeats
    check_revcomp=True
)

results = finder.find_alignments_between(upstream_region, downstream_region)

for r in results:
    print(f"Inverted repeat found:")
    print(f"  Upstream: {r['seq1']}")
    print(f"  Downstream: {r['seq2']} (reverse complement)")
```

### Example 3: Batch Processing Multiple Transposons

```python
from modules.short_alignment_finder import find_alignments_between_sequences

# Process multiple transposons
transposons = load_transposon_data()

results_summary = []

for transposon in transposons:
    alignments = find_alignments_between_sequences(
        transposon['flanking'],
        transposon['noncoding'],
        min_length=9,
        max_mismatches=1
    )

    # Get best match
    if alignments:
        best = max(alignments, key=lambda x: (x['length'], -x['mismatches']))
        results_summary.append({
            'transposon_id': transposon['id'],
            'best_match_length': best['length'],
            'best_match_mismatches': best['mismatches'],
            'best_match_seq': best['seq1']
        })

# Convert to DataFrame for analysis
import pandas as pd
df = pd.DataFrame(results_summary)
print(df)
```

### Example 4: With Mismatch Tolerance

```python
# Allow mismatches for biological variation
finder = ShortAlignmentFinder(min_length=9, max_mismatches=1)
results = finder.find_alignments_between(seq1, seq2)

# Filter by quality
high_quality = [r for r in results if r['mismatches'] == 0 and r['length'] >= 12]
medium_quality = [r for r in results if r['mismatches'] == 1 and r['length'] >= 12]

print(f"Perfect matches: {len(high_quality)}")
print(f"Matches with 1 mismatch: {len(medium_quality)}")
```

## Common Use Cases

### 1. IS110 Guide RNA Discovery

**Scenario**: Find guide RNA pairs between transposon flanking and non-coding regions

```python
finder = ShortAlignmentFinder(
    min_length=9,          # IS110 typical length
    max_mismatches=1,      # Allow minor variation
    check_forward=True,
    check_revcomp=True
)

results = finder.find_alignments_between(flanking, noncoding)

# Filter for guide RNA candidates
candidates = [
    r for r in results
    if 9 <= r['length'] <= 15  # Typical guide RNA length
    and r['mismatches'] <= 1
]
```

### 2. Regulatory Element Identification

**Scenario**: Find conserved regulatory motifs between upstream and downstream regions

```python
finder = ShortAlignmentFinder(
    min_length=6,          # Regulatory motifs can be shorter
    max_mismatches=0,      # Require perfect conservation
    check_forward=True,
    check_revcomp=False    # Only forward matches
)

conserved_motifs = finder.find_alignments_between(upstream, downstream)
```

### 3. Palindrome Detection Across Regions

**Scenario**: Find palindromic sequences that span different regions

```python
finder = ShortAlignmentFinder(
    min_length=8,
    max_mismatches=0,
    check_forward=False,   # Ignore direct matches
    check_revcomp=True     # Only palindromes
)

palindromes = finder.find_alignments_between(region1, region2)
```

## Integration with Other Modules

### With RegionExtractor

```python
from modules.region_extractor import RegionExtractor
from modules.short_alignment_finder import find_alignments_between_sequences

# Extract regions
extractor = RegionExtractor(fna_folder='/path/to/genomes/')

# Get flanking region
flanking_seq = extractor.extract_region(
    contig_id='CONTIG_123',
    start=1000,
    end=1100,
    strand=1
)

# Get non-coding region
noncoding_seq = extractor.extract_region(
    contig_id='CONTIG_123',
    start=2000,
    end=2100,
    strand=1
)

# Find alignments
alignments = find_alignments_between_sequences(flanking_seq, noncoding_seq, min_length=9)
```

### Export to DataFrame

```python
import pandas as pd

# Get results
results = find_alignments_between_sequences(seq1, seq2, min_length=9)

# Convert to DataFrame
df = pd.DataFrame(results)

# Add metadata
df['seq1_length'] = len(seq1)
df['seq2_length'] = len(seq2)
df['identity'] = 1 - (df['mismatches'] / df['length'])

# Analyze
print(df[['length', 'mismatches', 'orientation', 'identity']].describe())

# Export
df.to_csv('between_sequences_alignments.csv', index=False)
```

## Algorithm Details

### How It Works

1. **Seed Finding**
   - For each position `i` in sequence1, extract k-mer of `min_length`
   - Search for this k-mer in all positions of sequence2
   - Record seeds with ≤ `max_mismatches`

2. **Extension**
   - For each seed, extend left and right
   - Continue extending while total mismatches ≤ `max_mismatches`
   - Stop at sequence boundaries

3. **Consolidation**
   - Sort matches by length (descending) and mismatches (ascending)
   - Remove overlapping matches, keeping longest/best ones

### Performance

**Time Complexity**: O(n × m × k) where:
- n = length of sequence1
- m = length of sequence2
- k = average match length

**Typical Performance**:
- 100bp × 100bp: <0.1 seconds
- 500bp × 500bp: ~1-2 seconds
- 1kb × 1kb: ~5-10 seconds

**Optimization Tips**:
1. Use higher `min_length` to reduce seed matches
2. Set stricter `max_mismatches` for faster processing
3. Pre-filter sequences to remove low-complexity regions (N-rich)

## Troubleshooting

### Issue: Too Many Results

**Problem**: Finding hundreds of alignments, mostly low-quality

**Solutions**:
1. Increase `min_length` (e.g., from 9 to 12)
2. Decrease `max_mismatches` (e.g., from 1 to 0)
3. Post-filter by length and quality:
   ```python
   results = [r for r in results if r['length'] >= 12 and r['mismatches'] == 0]
   ```

### Issue: No Results Found

**Problem**: No alignments detected

**Solutions**:
1. Decrease `min_length` (e.g., from 15 to 9)
2. Increase `max_mismatches` (e.g., from 0 to 1)
3. Check both orientations:
   ```python
   finder = ShortAlignmentFinder(check_forward=True, check_revcomp=True)
   ```
4. Verify sequence format (should be A/T/G/C)

### Issue: Slow Performance

**Problem**: Takes too long with large sequences

**Solutions**:
1. Break sequences into smaller chunks
2. Use higher `min_length` to reduce seed matches
3. Pre-filter sequences:
   ```python
   # Remove N-rich regions
   seq1_filtered = remove_low_complexity(seq1)
   seq2_filtered = remove_low_complexity(seq2)
   ```

## Differences from Single-Sequence Mode

| Aspect | Single Sequence | Between Sequences |
|--------|----------------|-------------------|
| Distance constraints | Applied | **Ignored** |
| Self-matching | Needs handling | **Not applicable** |
| Position interpretation | Both in same sequence | **pos1 in seq1, pos2 in seq2** |
| Typical use | Find repeats | **Find shared motifs** |
| Output field | Includes 'distance' | **No 'distance' field** |

## Demo Script

Run the comprehensive demo to see all features:

```bash
python3 examples/demo_between_sequences_alignment.py
```

The demo includes:
1. Basic usage
2. Mismatch tolerance
3. Inverted repeats
4. IS110 transposon analysis
5. Batch processing
6. Convenience function usage
7. JSON export

## Testing

Run the test suite:

```bash
python3 -m pytest tests/test_short_alignment_finder.py::TestBetweenSequences -v
```

Tests cover:
- Basic between-sequences alignment
- Reverse complement matches
- Mismatch handling
- Auto-extension
- Transposon use cases
- Edge cases

## Comparison Table

### When to Use Each Mode

| Scenario | Mode | Reason |
|----------|------|--------|
| Find repeats within a transposon upstream region | **Single sequence** | Looking for patterns within one region |
| Find guide RNAs shared between flanking and non-coding | **Between sequences** | Comparing two distinct regions |
| Detect palindromic sequences in a gene | **Single sequence** | Within-region analysis |
| Compare upstream vs downstream regulatory elements | **Between sequences** | Cross-region comparison |
| Find tandem repeats | **Single sequence** | Within-region, distance matters |
| Identify conserved motifs across genomic features | **Between sequences** | Comparing separate features |

## Best Practices

1. **Choose appropriate parameters**
   - Start with `min_length=9, max_mismatches=0`
   - Adjust based on results

2. **Filter results**
   ```python
   # Keep only high-quality matches
   filtered = [r for r in results if r['length'] >= 12 and r['mismatches'] <= 1]
   ```

3. **Check both orientations**
   - IS110 guide RNAs can be in either orientation
   - Use `check_forward=True, check_revcomp=True`

4. **Validate biologically**
   - Cross-reference with known guide RNA databases
   - Check spacer distances in original genomic context

5. **Document your analysis**
   ```python
   import json

   analysis_record = {
       'parameters': {'min_length': 9, 'max_mismatches': 1},
       'sequences': {'seq1_len': len(seq1), 'seq2_len': len(seq2)},
       'results': results
   }

   with open('analysis.json', 'w') as f:
       json.dump(analysis_record, f, indent=2)
   ```

## Related Documentation

- [Short Alignment Finder (single sequence) README](SHORT_ALIGNMENT_FINDER_README.md)
- [Parser Quick Reference](PARSER_QUICK_REFERENCE.md)
- [Complete Workflow Guide](COMPLETE_WORKFLOW.md)
- [Module Organization](MODULE_ORGANIZATION.md)

## Author

Kuang Hu
Date: 2026-01-26

## Support

For issues, feature requests, or questions:
- Check existing documentation
- Run demo scripts for examples
- Review test cases for usage patterns
