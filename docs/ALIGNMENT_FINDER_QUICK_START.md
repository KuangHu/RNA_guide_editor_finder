# Short Alignment Finder - Quick Start Guide

## Two Modes of Operation

The Short Alignment Finder has **two modes** for different use cases:

| Mode | Use Case | Example |
|------|----------|---------|
| **Single Sequence** | Find repeats within one sequence | Find direct/inverted repeats in a transposon upstream region |
| **Between Sequences** | Find shared motifs between two sequences | Find guide RNAs shared between flanking and non-coding regions |

---

## Mode 1: Single Sequence (Find Repeats)

### When to Use
- Finding repeats within a transposon upstream region
- Detecting palindromic sequences
- Identifying tandem repeats
- Any within-sequence repeat analysis

### Quick Example

```python
from modules.short_alignment_finder import find_short_alignments

# Single sequence with repeats
sequence = "ATGC" + "ACGTACGTAC" + "NNNN" + "ACGTACGTAC" + "GCTA"

# Find repeats
results = find_short_alignments(
    sequence,
    min_length=9,
    max_mismatches=0,
    min_distance=10,      # Repeats must be at least 10bp apart
    max_distance=200,     # And within 200bp
    check_forward=True,   # Direct repeats
    check_revcomp=True    # Inverted repeats
)

# Output includes 'distance' field
print(results[0])
# {
#     'pos1': 4,
#     'pos2': 18,
#     'distance': 14,          # pos2 - pos1
#     'length': 10,
#     'seq1': 'ACGTACGTAC',
#     'seq2': 'ACGTACGTAC',
#     'mismatches': 0,
#     'mismatch_positions': [],
#     'orientation': 'forward'
# }
```

### Key Features
- ✅ Distance constraints (min_distance, max_distance)
- ✅ Find repeats at specific spacing
- ✅ Useful for structural analysis within a region
- ✅ Output includes 'distance' field

---

## Mode 2: Between Sequences (Find Shared Motifs)

### When to Use
- Finding guide RNAs between flanking and non-coding regions
- Comparing upstream vs downstream regions
- Identifying conserved motifs across genomic features
- Any cross-region comparison

### Quick Example

```python
from modules.short_alignment_finder import find_alignments_between_sequences

# Two different sequences
flanking_region = "ATGC" + "TTTAAACCCGGG" + "GCTA"
noncoding_region = "NNNN" + "TTTAAACCCGGG" + "NNNN"

# Find shared sequences
results = find_alignments_between_sequences(
    flanking_region,
    noncoding_region,
    min_length=9,
    max_mismatches=0,
    check_forward=True,   # Direct matches
    check_revcomp=True    # Reverse complement matches
)

# Output: pos1 in seq1, pos2 in seq2 (no distance field)
print(results[0])
# {
#     'pos1': 4,            # Position in flanking_region
#     'pos2': 4,            # Position in noncoding_region
#     'length': 12,
#     'seq1': 'TTTAAACCCGGG',
#     'seq2': 'TTTAAACCCGGG',
#     'mismatches': 0,
#     'mismatch_positions': [],
#     'orientation': 'forward'
# }
```

### Key Features
- ✅ Compare two different sequences
- ✅ No distance constraints (not applicable)
- ✅ pos1 and pos2 refer to different sequences
- ✅ No 'distance' field in output

---

## Side-by-Side Comparison

| Feature | Single Sequence | Between Sequences |
|---------|----------------|-------------------|
| **Function** | `find_short_alignments(seq)` | `find_alignments_between_sequences(seq1, seq2)` |
| **Input** | One sequence | Two sequences |
| **min_distance** | ✅ Applies | ❌ Ignored |
| **max_distance** | ✅ Applies | ❌ Ignored |
| **Output: pos1** | Position in sequence | Position in sequence1 |
| **Output: pos2** | Position in sequence (later) | Position in sequence2 |
| **Output: distance** | ✅ Included | ❌ Not included |
| **Use case** | Find repeats | Find shared motifs |

---

## Complete Examples

### IS110 Analysis: Finding Repeats Within Upstream Region

```python
from modules.short_alignment_finder import ShortAlignmentFinder

# Upstream region of transposon
upstream_region = "ATGC...TGCA"  # Your sequence

# Find repeats within this region
finder = ShortAlignmentFinder(
    min_length=9,
    max_mismatches=1,
    min_distance=20,      # Repeats separated by spacers
    max_distance=200,     # Within reasonable window
    check_forward=True,
    check_revcomp=True
)

results = finder.find_alignments(upstream_region)

print(f"Found {len(results)} repeats within upstream region")
for r in results:
    print(f"  {r['length']}bp at positions {r['pos1']} and {r['pos2']}, distance={r['distance']}")
```

### IS110 Analysis: Finding Guide RNAs Between Regions

```python
from modules.short_alignment_finder import ShortAlignmentFinder

# Two different regions
flanking_region = extract_flanking(transposon_id)
noncoding_region = extract_noncoding(transposon_id)

# Find shared sequences between regions
finder = ShortAlignmentFinder(
    min_length=9,
    max_mismatches=1,
    check_forward=True,
    check_revcomp=True
)

results = finder.find_alignments_between(flanking_region, noncoding_region)

print(f"Found {len(results)} shared sequences between flanking and non-coding")
for r in results:
    print(f"  {r['length']}bp: flanking pos {r['pos1']}, noncoding pos {r['pos2']}")
```

---

## Parameter Guide

### Common Parameters (Both Modes)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_length` | int | 9 | Minimum alignment length (bp) |
| `max_mismatches` | int | 0 | Maximum mismatches allowed |
| `check_forward` | bool | True | Find direct matches |
| `check_revcomp` | bool | True | Find reverse complement matches |

### Single-Sequence Only Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_distance` | int | 0 | Minimum distance between repeats (pos2 - pos1) |
| `max_distance` | int | None | Maximum distance (None = no limit) |

**Note**: These are **ignored** in between-sequences mode.

---

## Choosing the Right Mode

### Use Single Sequence When:
- ✅ Analyzing one genomic region
- ✅ Looking for repeats at specific spacing
- ✅ Distance between matches matters
- ✅ Finding tandem repeats or palindromes within a region

**Example**: "Find all 9-15bp repeats in the transposon upstream region that are 20-200bp apart"

### Use Between Sequences When:
- ✅ Comparing two different genomic regions
- ✅ Finding shared motifs across features
- ✅ Distance doesn't make sense (different sequences)
- ✅ Looking for conserved elements

**Example**: "Find guide RNA sequences that appear in both the flanking region and the non-coding region"

---

## Output Format

### Single Sequence Output

```python
{
    'pos1': 10,                    # First occurrence (0-based)
    'pos2': 45,                    # Second occurrence (0-based)
    'distance': 35,                # pos2 - pos1
    'length': 12,
    'seq1': 'ACGTACGTACGT',
    'seq2': 'ACGTACGTACGT',
    'mismatches': 0,
    'mismatch_positions': [],
    'orientation': 'forward'
}
```

### Between Sequences Output

```python
{
    'pos1': 10,                    # Position in sequence1 (0-based)
    'pos2': 25,                    # Position in sequence2 (0-based)
    # No 'distance' field
    'length': 12,
    'seq1': 'ACGTACGTACGT',
    'seq2': 'ACGTACGTACGT',
    'mismatches': 0,
    'mismatch_positions': [],
    'orientation': 'forward'
}
```

---

## Common Patterns

### Pattern 1: Strict Perfect Matches

```python
# Single sequence
results = find_short_alignments(seq, min_length=12, max_mismatches=0)

# Between sequences
results = find_alignments_between_sequences(seq1, seq2, min_length=12, max_mismatches=0)
```

### Pattern 2: Allow Biological Variation

```python
# Single sequence
results = find_short_alignments(seq, min_length=9, max_mismatches=1)

# Between sequences
results = find_alignments_between_sequences(seq1, seq2, min_length=9, max_mismatches=1)
```

### Pattern 3: Only Inverted Repeats

```python
# Single sequence
results = find_short_alignments(seq, min_length=8, check_forward=False, check_revcomp=True)

# Between sequences
results = find_alignments_between_sequences(seq1, seq2, min_length=8,
                                           check_forward=False, check_revcomp=True)
```

### Pattern 4: Distance-Constrained (Single Sequence Only)

```python
# Find repeats 50-150bp apart
results = find_short_alignments(seq, min_length=10, min_distance=50, max_distance=150)
```

---

## Class vs Function Usage

### Using the Class (More Control)

```python
from modules.short_alignment_finder import ShortAlignmentFinder

# Create finder with parameters
finder = ShortAlignmentFinder(min_length=9, max_mismatches=1)

# Use for single sequence
results1 = finder.find_alignments(sequence)

# Reuse for between sequences
results2 = finder.find_alignments_between(seq1, seq2)
```

### Using Functions (Quick & Simple)

```python
from modules.short_alignment_finder import (
    find_short_alignments,
    find_alignments_between_sequences
)

# One-line calls
results1 = find_short_alignments(sequence, min_length=9)
results2 = find_alignments_between_sequences(seq1, seq2, min_length=9)
```

---

## Integration Example: Complete IS110 Workflow

```python
from modules.short_alignment_finder import ShortAlignmentFinder
from modules.region_extractor import RegionExtractor

# Setup
extractor = RegionExtractor(fna_folder='/path/to/genomes/')
finder = ShortAlignmentFinder(min_length=9, max_mismatches=1)

# Extract regions
upstream = extractor.extract_region('CONTIG_123', 1000, 1200, strand=1)
flanking = extractor.extract_region('CONTIG_123', 1200, 1300, strand=1)
noncoding = extractor.extract_region('CONTIG_123', 2000, 2100, strand=1)

# Analysis 1: Find repeats within upstream region
upstream_repeats = finder.find_alignments(upstream)
print(f"Found {len(upstream_repeats)} repeats in upstream region")

# Analysis 2: Find guide RNAs between flanking and non-coding
guide_rnas = finder.find_alignments_between(flanking, noncoding)
print(f"Found {len(guide_rnas)} potential guide RNA pairs")

# Filter for high-quality candidates
candidates = [r for r in guide_rnas if r['length'] >= 12 and r['mismatches'] == 0]
print(f"High-quality candidates: {len(candidates)}")
```

---

## Demo Scripts

### Single Sequence Demo
```bash
python3 examples/demo_short_alignment_finder.py
```

### Between Sequences Demo
```bash
python3 examples/demo_between_sequences_alignment.py
```

---

## Further Reading

- **[SHORT_ALIGNMENT_FINDER_README.md](SHORT_ALIGNMENT_FINDER_README.md)** - Complete single-sequence documentation
- **[BETWEEN_SEQUENCES_GUIDE.md](BETWEEN_SEQUENCES_GUIDE.md)** - Complete between-sequences documentation
- **[COMPLETE_WORKFLOW.md](COMPLETE_WORKFLOW.md)** - Full IS110 analysis pipeline

---

## Quick Decision Tree

```
Are you comparing two different sequences?
├─ YES → Use find_alignments_between_sequences()
│        (flanking vs non-coding, upstream vs downstream, etc.)
│
└─ NO  → Use find_short_alignments()
         (find repeats within a single region)

Do you care about distance between matches?
├─ YES → Use find_short_alignments() with min_distance/max_distance
│        (e.g., "repeats must be 20-200bp apart")
│
└─ NO  → Either mode works, choose based on whether you have
         one sequence or two sequences
```

---

## Author

Kuang Hu
Date: 2026-01-26
