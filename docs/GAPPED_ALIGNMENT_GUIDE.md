# Gapped Alignment Extension Guide

## Overview

The gapped alignment extension feature allows you to extend existing ungapped alignments with gaps (insertions/deletions) allowed. This is particularly useful for finding **CAST homing spacer patterns** where sequences have ~75-85% identity with mismatches and indels.

## Key Features

- **Seed-based extension**: Uses existing ungapped alignments as anchor points
- **Bidirectional extension**: Extends both 5' and 3' from seed position
- **Affine gap penalties**: Separate penalties for gap opening and extension
- **Identity filtering**: Configurable minimum identity threshold
- **Gap support**: Handles insertions, deletions, and mismatches
- **Visual output**: Provides alignment strings for visualization

## When to Use

Use gapped alignment extension when you need to:
- Find CAST homing spacer-like patterns (17/20, 27/32 matches)
- Detect alignments with indels (not just substitutions)
- Extend local alignments to include surrounding mismatched regions
- Compare spacer sequences to att sites with gaps

## Basic Workflow

### 1. Find Ungapped Seeds

First, find short ungapped (or nearly ungapped) alignments to use as seeds:

```python
from modules.short_alignment_finder import find_alignments_between_sequences

spacer = "AAACGAACCACCTTCATTCGAGTACAAGAGCA"
attsite = "AAGCGTACAACGTTTATCCGAGTTCAAGAGCA"

# Find seeds with up to 2 mismatches
seeds = find_alignments_between_sequences(
    spacer,
    attsite,
    min_length=8,
    max_mismatches=2
)
```

### 2. Extend with Gaps

Extend the best seed with gaps allowed:

```python
from modules.short_alignment_finder import extend_alignment_with_gaps

# Extend best seed
gapped = extend_alignment_with_gaps(
    spacer,
    attsite,
    seeds[0],
    min_identity=0.75,  # 75% minimum identity
    max_extension=40
)

if gapped:
    print(f"Identity: {gapped['identity']:.1%}")
    print(f"Pattern: {gapped['matches']}/{gapped['alignment_length']}")
    print(f"Alignment: {gapped['seq1_aligned']}")
    print(f"           {gapped['alignment_string']}")
    print(f"           {gapped['seq2_aligned']}")
```

## Parameters

### Required Parameters

- `sequence1`: First DNA sequence (str)
- `sequence2`: Second DNA sequence (str)
- `seed_match`: Seed alignment dictionary from `find_alignments` or `find_alignments_between_sequences`

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_identity` | 0.80 | Minimum identity (matches/total) required (0.0-1.0) |
| `max_extension` | 50 | Maximum bases to extend in each direction |
| `match_score` | 2 | Score for matching bases |
| `mismatch_penalty` | -1 | Penalty for mismatches |
| `gap_open_penalty` | -3 | Penalty for opening a gap |
| `gap_extend_penalty` | -1 | Penalty for extending a gap |

## Return Value

Returns a dictionary with:

```python
{
    'pos1': 0,                      # Start in sequence1
    'pos2': 0,                      # Start in sequence2
    'end1': 32,                     # End in sequence1
    'end2': 32,                     # End in sequence2
    'length1': 32,                  # Length in sequence1
    'length2': 32,                  # Length in sequence2
    'seq1_aligned': 'AAACGAA...',  # Aligned seq1 with gaps
    'seq2_aligned': 'AAGCGTA...',  # Aligned seq2 with gaps
    'alignment_length': 32,         # Total alignment length
    'matches': 25,                  # Number of matches
    'mismatches': 7,                # Number of mismatches
    'gaps': 0,                      # Number of gaps
    'identity': 0.781,              # Identity fraction
    'alignment_string': '||.||...', # Visual alignment
    'orientation': 'forward',       # 'forward' or 'reverse_complement'
    'seed_pos1': 0,                 # Original seed start
    'seed_pos2': 0,                 # Original seed start
    'seed_length': 8                # Original seed length
}
```

Returns `None` if extension doesn't meet `min_identity` threshold.

## Real CAST Examples

Here are results from real CAST homing spacer alignments:

### Example 1: 78% Identity

```
Spacer:  AAACGAACCACCTTCATTCGAGTACAAGAGCA
         ||.||.||.||.||.||.|||||.||||||||
AttSite: AAGCGTACAACGTTTATCCGAGTTCAAGAGCA

Identity: 78.1% (25/32)
Matches: 25, Mismatches: 7, Gaps: 0
```

### Example 2: 81% Identity

```
Spacer:  AAATTTACATCTCCAGGCTCAGCCAAAAAGC
         |||||.||.||.|||||.|||||.|||||.|
AttSite: AAATTAACGTCGCCAGGTTCAGCTAAAAAAC

Identity: 80.6% (25/31)
Matches: 25, Mismatches: 6, Gaps: 0
```

### Example 3: 75% Identity

```
Spacer:  AGGACAGGAAGAAAACACCCAAGTTGGG
         |||||.|||||..|.||.|||||..|||
AttSite: AGGACCGGAAGGTAGCAGCCAAGGCGGG

Identity: 75.0% (21/28)
Matches: 21, Mismatches: 7, Gaps: 0
```

### Example 4: 76% Identity

```
Spacer:  AAAAACAATCGATTCGATGTTTAAGTTAA
         |||||.||.||.||.|||||.||.||.||
AttSite: AAAAAGAACCGTTTTGATGTCTATGTAAA

Identity: 75.9% (22/29)
Matches: 22, Mismatches: 7, Gaps: 0
```

## Tips and Best Practices

### Seed Selection

1. **Use flexible seed parameters**: Allow 1-3 mismatches in seeds for better coverage
   ```python
   seeds = find_alignments_between_sequences(
       seq1, seq2,
       min_length=8,
       max_mismatches=2  # Allow some mismatches
   )
   ```

2. **Try multiple seeds**: Extend several seeds and pick the best result
   ```python
   best_result = None
   best_identity = 0

   for seed in seeds[:5]:  # Try top 5 seeds
       result = extend_alignment_with_gaps(seq1, seq2, seed)
       if result and result['identity'] > best_identity:
           best_result = result
           best_identity = result['identity']
   ```

### Identity Thresholds

Choose identity threshold based on your requirements:

- **≥85%**: Very stringent, for highly similar sequences
- **75-85%**: Typical CAST homing spacer range (recommended)
- **70-75%**: More permissive, may include false positives
- **<70%**: Very permissive, likely to include many false matches

### Gap Penalties

Adjust gap penalties to control gap vs. mismatch preference:

```python
# Prefer gaps over mismatches (lenient gap penalty)
extend_alignment_with_gaps(
    seq1, seq2, seed,
    gap_open_penalty=-1,   # Less penalty
    gap_extend_penalty=-1
)

# Prefer mismatches over gaps (strict gap penalty)
extend_alignment_with_gaps(
    seq1, seq2, seed,
    gap_open_penalty=-5,   # More penalty
    gap_extend_penalty=-2
)
```

### Performance Tips

1. **Limit max_extension**: Use realistic extension limits
   ```python
   max_extension=40  # For ~30bp sequences
   ```

2. **Filter seeds first**: Only extend high-quality seeds
   ```python
   # Sort seeds by length and mismatches
   sorted_seeds = sorted(seeds, key=lambda x: (-x['length'], x['mismatches']))
   best_seeds = sorted_seeds[:3]
   ```

3. **Use appropriate min_length**: Longer seeds = faster
   ```python
   min_length=9  # 9-11bp works well for most cases
   ```

## Working with RNA Sequences

If you have RNA sequences (with U instead of T):

```python
# Convert U to T
seq1_dna = seq1.replace('U', 'T')
seq2_dna = seq2.replace('U', 'T')

# Then proceed normally
seeds = find_alignments_between_sequences(seq1_dna, seq2_dna, ...)
```

## Complete Example

```python
#!/usr/bin/env python3
from modules.short_alignment_finder import (
    find_alignments_between_sequences,
    extend_alignment_with_gaps
)

# Real CAST homing spacer example
spacer = "AAACGAACCACCTTCATTCGAGTACAAGAGCA"
attsite = "AAGCGTACAACGTTTATCCGAGTTCAAGAGCA"

# Step 1: Find seeds
print("Finding seeds...")
seeds = find_alignments_between_sequences(
    spacer, attsite,
    min_length=8,
    max_mismatches=2
)
print(f"Found {len(seeds)} seeds")

# Step 2: Extend best seeds
print("\nExtending seeds...")
best_result = None
best_identity = 0

for i, seed in enumerate(seeds[:3], 1):
    result = extend_alignment_with_gaps(
        spacer, attsite, seed,
        min_identity=0.70,
        max_extension=40
    )

    if result and result['identity'] > best_identity:
        best_result = result
        best_identity = result['identity']
        print(f"  Seed {i}: {result['matches']}/{result['alignment_length']} "
              f"({result['identity']:.1%})")

# Step 3: Display best result
if best_result:
    print(f"\nBest alignment:")
    print(f"  {best_result['seq1_aligned']}")
    print(f"  {best_result['alignment_string']}")
    print(f"  {best_result['seq2_aligned']}")
    print(f"\nStatistics:")
    print(f"  Identity: {best_result['identity']:.1%}")
    print(f"  Matches: {best_result['matches']}")
    print(f"  Mismatches: {best_result['mismatches']}")
    print(f"  Gaps: {best_result['gaps']}")
else:
    print("No alignment met threshold")
```

## See Also

- `SHORT_ALIGNMENT_FINDER_README.md` - Ungapped alignment finding
- `BETWEEN_SEQUENCES_GUIDE.md` - Between-sequence alignment basics
- `examples/demo_gapped_alignment.py` - Interactive examples
- `tests/test_gapped_alignment.py` - Test suite

## References

- CAST homing spacer patterns typically show 70-85% identity
- Affine gap penalties: gap_cost = open + (n-1) × extend
- Semi-global alignment: no penalty for terminal gaps
