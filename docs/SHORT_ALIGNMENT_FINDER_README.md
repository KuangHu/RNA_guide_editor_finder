# Short Alignment Finder Module

## Overview

The **Short Alignment Finder** is a module for identifying short DNA alignments (repeats) within a single DNA sequence. It finds both **direct repeats** (forward-forward matches) and **inverted repeats** (forward-reverse complement matches) with configurable parameters.

This is particularly useful for:
- Finding guide RNA pairs in IS110 transposon upstream regions
- Identifying potential regulatory elements
- Analyzing repeat structure in genomic sequences
- Detecting palindromic sequences

## Key Features

✅ **Configurable Parameters**
- Minimum alignment length (e.g., 9bp for IS110)
- Mismatch tolerance (0, 1, or more mismatches)
- Distance constraints (min/max distance between repeats)
- Separate control for forward and reverse complement matches

✅ **Auto-Extension**
- Automatically extends seed matches to their maximum length
- If a 9bp match is part of an 11bp match, returns the 11bp match
- Consolidates overlapping matches

✅ **Detailed Output**
- Both sequences (even if they have mismatches)
- Mismatch positions for quality assessment
- Orientation (forward or reverse_complement)
- Full coordinate information

✅ **Performance**
- Efficient seed-and-extend algorithm
- De-duplication to avoid redundant results

## Installation

The module is part of the RNA Guide Editor Finder toolkit. No additional installation required.

```python
from modules.short_alignment_finder import ShortAlignmentFinder, find_short_alignments
```

## Quick Start

### Example 1: Basic Usage

```python
from modules.short_alignment_finder import find_short_alignments

# DNA sequence with a direct repeat
sequence = "ACGTACGTAC" + "NNNN" + "ACGTACGTAC"

# Find repeats with minimum 9bp length
results = find_short_alignments(sequence, min_length=9)

# Output format:
# [
#     {
#         'pos1': 0,
#         'pos2': 14,
#         'length': 10,
#         'distance': 14,
#         'seq1': 'ACGTACGTAC',
#         'seq2': 'ACGTACGTAC',
#         'mismatches': 0,
#         'mismatch_positions': [],
#         'orientation': 'forward'
#     }
# ]
```

### Example 2: IS110 Guide RNA Finding

```python
from modules.short_alignment_finder import ShortAlignmentFinder

# IS110 upstream region
upstream_sequence = "ATGC...TGCA"  # Your sequence here

# IS110-specific parameters
finder = ShortAlignmentFinder(
    min_length=9,          # IS110 guide RNAs typically 9-15bp
    max_mismatches=1,      # Allow 1 mismatch for biological variation
    min_distance=20,       # Repeats separated by spacers
    max_distance=200,      # Within reasonable genomic window
    check_forward=True,    # Find direct repeats
    check_revcomp=True     # Find inverted repeats
)

results = finder.find_alignments(upstream_sequence)

# Analyze results
for result in results:
    print(f"Found {result['length']}bp repeat at positions {result['pos1']} and {result['pos2']}")
    print(f"  Orientation: {result['orientation']}")
    print(f"  Sequence 1: {result['seq1']}")
    print(f"  Sequence 2: {result['seq2']}")
    if result['mismatches'] > 0:
        print(f"  Mismatches: {result['mismatches']} at positions {result['mismatch_positions']}")
```

## API Reference

### ShortAlignmentFinder Class

```python
finder = ShortAlignmentFinder(
    min_length=9,                # Minimum alignment length to report
    max_mismatches=0,            # Maximum number of mismatches allowed
    min_distance=0,              # Minimum distance between positions (pos2 - pos1)
    max_distance=None,           # Maximum distance (None = no limit)
    check_forward=True,          # Find direct repeats
    check_revcomp=True           # Find inverted repeats
)
```

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_length` | int | 9 | Minimum alignment length to report (bp) |
| `max_mismatches` | int | 0 | Maximum mismatches allowed |
| `min_distance` | int | 0 | Minimum distance between alignment starts |
| `max_distance` | int or None | None | Maximum distance (None = no limit) |
| `check_forward` | bool | True | Find direct repeats (forward-forward) |
| `check_revcomp` | bool | True | Find inverted repeats (forward-RC) |

#### Methods

**`find_alignments(sequence: str) -> List[Dict]`**

Find all short alignments in the sequence.

**Parameters:**
- `sequence`: DNA sequence to search (string, A/T/G/C, case-insensitive)

**Returns:**
- List of alignment dictionaries (see Output Format below)

### Convenience Function

```python
results = find_short_alignments(
    sequence,
    min_length=9,
    max_mismatches=0,
    min_distance=0,
    max_distance=None,
    check_forward=True,
    check_revcomp=True
)
```

Same parameters as `ShortAlignmentFinder`, but as a one-line function call.

## Output Format

Each alignment is returned as a dictionary:

```python
{
    'pos1': 153,                    # Start position of first occurrence (0-based)
    'pos2': 241,                    # Start position of second occurrence (0-based)
    'length': 10,                   # Length of the alignment (bp)
    'distance': 88,                 # Distance between positions (pos2 - pos1)
    'seq1': 'ACACCCCTTG',          # Sequence at pos1
    'seq2': 'ACACCCCTTG',          # Sequence at pos2 (as it appears in input)
    'mismatches': 0,                # Number of mismatches
    'mismatch_positions': [],       # List of mismatch positions (0-based within alignment)
    'orientation': 'forward'        # 'forward' or 'reverse_complement'
}
```

### Interpretation

- **pos1, pos2**: 0-based coordinates in the input sequence
- **distance**: Measured start-to-start (pos2 - pos1)
- **seq2**: Always shown as it appears in the input sequence
  - For `'forward'` orientation: seq2 matches seq1 directly
  - For `'reverse_complement'` orientation: seq2 is the reverse complement of seq1
- **mismatch_positions**: Positions within the alignment (0 to length-1) where bases differ

## Distance Constraints

Distance is measured **start-to-start** (pos2 - pos1), not gap distance.

```
Example sequence:
Position:  0         10        20
           ACGTACGTAC|NNNN|ACGTACGTAC

pos1 = 0
pos2 = 14
distance = pos2 - pos1 = 14

Gap distance (end-to-start) = pos2 - (pos1 + length) = 14 - 10 = 4
```

### Common Distance Settings

| Setting | Use Case |
|---------|----------|
| `min_distance=0, max_distance=None` | Find all repeats (default) |
| `min_distance=20, max_distance=200` | IS110 guide RNAs (typical spacer range) |
| `min_distance=10` | Exclude very close/overlapping repeats |
| `max_distance=50` | Only nearby repeats |

## Understanding Orientations

### Direct Repeats (forward)

Both occurrences in the same orientation:

```
Sequence: ACGTACGT----ACGTACGT
          ^^^^^^^^    ^^^^^^^^
          pos1=0      pos2=12
          orientation='forward'
```

### Inverted Repeats (reverse_complement)

Second occurrence is the reverse complement of the first:

```
Sequence: ACGTACGT----ACGTACGT
          ^^^^^^^^    ^^^^^^^^
          pos1=0      pos2=12 (shown as appears in sequence)

Comparison:
  pos1: ACGTACGT
  pos2: ACGTACGT (reverse complement of ACGTACGT = ACGTACGT, palindrome!)

  orientation='reverse_complement'
```

## Mismatch Handling

### Recording Mismatches

When `max_mismatches > 0`, alignments with mismatches are reported:

```python
# Sequence with 1 mismatch at position 4
sequence = "ACGTACGTAC" + "NNNN" + "ACGTGCGTAC"
#           ACGTACGTAC              ACGTGCGTAC
#               ^                       ^ (A->G mismatch)

finder = ShortAlignmentFinder(min_length=9, max_mismatches=1)
results = finder.find_alignments(sequence)

# Output:
# {
#     'length': 10,
#     'seq1': 'ACGTACGTAC',
#     'seq2': 'ACGTGCGTAC',
#     'mismatches': 1,
#     'mismatch_positions': [4],  # Position within the alignment
#     ...
# }
```

### Extension with Mismatches

The algorithm extends matches **even when encountering mismatches**, as long as the total stays within `max_mismatches`.

Example with `max_mismatches=1`:
1. Find 9bp seed match with 0 mismatches
2. Extend right: next base is a mismatch, total = 1 (OK, continue)
3. Extend right: next base matches (total still = 1, continue)
4. Extend right: next base is a mismatch, total = 2 (STOP, limit exceeded)

## Auto-Extension Feature

The algorithm automatically extends matches to their maximum length and consolidates overlapping matches.

### Example

```python
# A 15bp perfect repeat
sequence = "ACGTACGTACGTACG" + "NNNN" + "ACGTACGTACGTACG"

finder = ShortAlignmentFinder(min_length=9)
results = finder.find_alignments(sequence)

# Even though min_length=9, the algorithm finds and reports the full 15bp match
# It does NOT report multiple overlapping 9bp, 10bp, 11bp, etc. matches
# Result: 1 alignment of length 15bp
```

### How It Works

1. **Seed Finding**: Find all possible 9bp matches (min_length)
2. **Extension**: For each seed, extend left and right as far as possible
3. **Consolidation**: If multiple extended matches overlap, keep the longest one

This ensures you get biologically meaningful results without redundant shorter matches.

## Performance Considerations

### Computational Complexity

- **Seed finding**: O(n²) where n = sequence length
- **Extension**: O(m) where m = match length
- **Consolidation**: O(k²) where k = number of matches (typically small)

### Optimization Tips

1. **Use distance constraints**: Setting `max_distance` significantly reduces search space
2. **Appropriate min_length**: Higher values reduce number of seed matches
3. **Avoid N-rich sequences**: Ambiguous bases (N) match everything, creating many false positives

### Typical Performance

- 1kb sequence: <1 second
- 10kb sequence: ~1-2 seconds
- 100kb sequence: ~30-60 seconds (consider breaking into chunks)

## Common Use Cases

### 1. IS110 Guide RNA Identification

```python
finder = ShortAlignmentFinder(
    min_length=9,
    max_mismatches=1,
    min_distance=20,
    max_distance=200,
    check_forward=True,
    check_revcomp=True
)
```

### 2. Perfect Palindrome Detection

```python
# Only look for inverted repeats with no mismatches
finder = ShortAlignmentFinder(
    min_length=8,
    max_mismatches=0,
    check_forward=False,      # Disable direct repeats
    check_revcomp=True        # Only inverted repeats
)
```

### 3. Nearby Tandem Repeats

```python
# Short-range repeats only
finder = ShortAlignmentFinder(
    min_length=6,
    max_mismatches=0,
    min_distance=5,
    max_distance=50
)
```

### 4. Distant Repeats Only

```python
# Long-range repeats (e.g., transposon ends)
finder = ShortAlignmentFinder(
    min_length=15,
    max_mismatches=2,
    min_distance=1000,      # At least 1kb apart
    max_distance=None       # No upper limit
)
```

## Integration with Other Modules

### Export to DataFrame

```python
import pandas as pd
from modules.short_alignment_finder import find_short_alignments

# Find alignments
results = find_short_alignments(sequence, min_length=9)

# Convert to DataFrame
df = pd.DataFrame(results)

# Analyze
print(df[['pos1', 'pos2', 'length', 'distance', 'orientation', 'mismatches']])

# Filter
high_quality = df[(df['length'] >= 12) & (df['mismatches'] == 0)]

# Export
df.to_csv('alignments.csv', index=False)
```

### Batch Processing

```python
from utils.parsers import parse_fasta

# Process multiple sequences
for record in parse_fasta('sequences.fasta'):
    results = find_short_alignments(record.sequence, min_length=9)
    print(f"{record.id}: Found {len(results)} alignments")
    for r in results:
        print(f"  - {r['length']}bp at positions {r['pos1']}, {r['pos2']}")
```

## Troubleshooting

### Issue: Too Many Results

**Problem**: Thousands of alignments found

**Solutions**:
1. Increase `min_length` (e.g., from 9 to 12)
2. Set stricter distance constraints (`min_distance`, `max_distance`)
3. Reduce `max_mismatches` (e.g., from 1 to 0)
4. Filter N-rich regions before analysis

### Issue: No Results Found

**Problem**: No alignments detected

**Solutions**:
1. Decrease `min_length` (e.g., from 15 to 9)
2. Increase `max_mismatches` (e.g., from 0 to 1)
3. Check both orientations (`check_forward=True, check_revcomp=True`)
4. Verify sequence format (should be A/T/G/C)

### Issue: Unexpected Short Matches

**Problem**: Getting many 9bp matches instead of longer ones

**Possible causes**:
- Sequence has many mismatches beyond the auto-extension limit
- True biological variation (not all repeats are perfect)

**Solutions**:
- Check `mismatch_positions` to understand where differences occur
- Increase `max_mismatches` to allow more variation during extension
- Filter results by length in downstream analysis

## Examples

See `examples/demo_short_alignment_finder.py` for comprehensive examples including:
1. Basic usage with perfect direct repeats
2. Finding repeats with mismatches
3. Inverted repeat detection
4. Distance constraints
5. Auto-extension demonstration
6. JSON output
7. Complete IS110 workflow

Run the demo:
```bash
python3 examples/demo_short_alignment_finder.py
```

## Testing

Run the test suite:
```bash
python3 -m pytest tests/test_short_alignment_finder.py -v
```

Tests cover:
- Direct and inverted repeat finding
- Auto-extension and consolidation
- Mismatch tolerance
- Distance constraints
- Edge cases and parameter validation

## Author

Kuang Hu
Date: 2026-01-26

## Related Modules

- **sequence_search.py**: Find sequences in genomic databases using MMseqs2
- **alignment_filter.py**: Filter alignments by composition, length, topology
- **region_extractor.py**: Extract genomic regions with flanking sequences
- **utils/parsers.py**: Parse bioinformatics file formats (includes `reverse_complement()`)
