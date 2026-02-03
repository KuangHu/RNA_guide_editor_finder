# Parser Module Quick Reference

## Import Everything You Need

```python
from utils.parsers import (
    # Data classes
    DiamondBlastHit,
    HmmHit,
    # Parsers
    parse_prodigal_faa,
    parse_diamond_blastp_with_positions,
    parse_hmm_tblout_with_positions,
    # Sequence retrieval
    retrieve_sequences_from_hits,
    retrieve_sequence_from_fasta,
    reverse_complement
)
```

## Typical Workflow

### 1. Diamond BLASTP Workflow

```python
# Step 1: Parse Prodigal output for positions
position_map = parse_prodigal_faa("genome.faa")

# Step 2: Parse Diamond BLASTP with positions
hits = parse_diamond_blastp_with_positions(
    "diamond_results.m8",
    position_map,
    file_basename="genome_name"
)

# Step 3: Retrieve sequences
hits = retrieve_sequences_from_hits(hits, ["genome.fna"])

# Step 4: Use the data
for hit in hits:
    print(f"{hit.file_basename}: {hit.contig_id}:{hit.start}-{hit.end}")
    print(f"  Protein: {hit.protein_id}")
    print(f"  Subject: {hit.subject_id}")
    print(f"  E-value: {hit.evalue}, Identity: {hit.pident}%")
    if hasattr(hit, 'genomic_sequence'):
        print(f"  Sequence: {hit.genomic_sequence}")
```

### 2. HMMER Workflow

```python
# Same as Diamond, just use HMMER parser
position_map = parse_prodigal_faa("genome.faa")

hits = parse_hmm_tblout_with_positions(
    "hmmer_results.tbl",
    position_map,
    file_basename="genome_name"
)

hits = retrieve_sequences_from_hits(hits, ["genome.fna"])

for hit in hits:
    print(f"HMM: {hit.query_name}, Score: {hit.score}")
```

## Data Available in Each Hit

### DiamondBlastHit
- `hit.file_basename` - Which file this came from
- `hit.contig_id` - Which contig/chromosome
- `hit.protein_id` - Protein identifier
- `hit.start` - Start position (1-based)
- `hit.end` - End position (1-based)
- `hit.strand` - 1 (forward) or -1 (reverse)
- `hit.subject_id` - What it matched
- `hit.pident` - Percent identity
- `hit.evalue` - E-value
- `hit.bitscore` - Bit score
- `hit.genomic_sequence` - DNA sequence (after retrieval)
- `hit.length` - Length of region
- `hit.is_forward` - True if forward strand

### HmmHit
- `hit.file_basename` - Which file this came from
- `hit.contig_id` - Which contig/chromosome
- `hit.protein_id` - Protein identifier
- `hit.start` - Start position (1-based)
- `hit.end` - End position (1-based)
- `hit.strand` - 1 (forward) or -1 (reverse)
- `hit.query_name` - HMM profile name
- `hit.evalue` - E-value
- `hit.score` - Score
- `hit.genomic_sequence` - DNA sequence (after retrieval)
- `hit.length` - Length of region
- `hit.is_forward` - True if forward strand

## Filtering Examples

```python
# Filter by E-value
good_hits = [h for h in hits if h.evalue < 1e-50]

# Filter by identity (Diamond only)
high_identity = [h for h in hits if h.pident > 90]

# Filter by score (HMMER only)
high_score = [h for h in hits if h.score > 100]

# Filter by strand
forward_only = [h for h in hits if h.is_forward]

# Filter by length
long_enough = [h for h in hits if h.length >= 100]
```

## Memory-Efficient Processing

For very large files, use iterators:

```python
from utils.parsers import iter_diamond_blastp_with_positions

position_map = parse_prodigal_faa("genome.faa")

# Process one hit at a time
for hit in iter_diamond_blastp_with_positions("results.m8", position_map):
    if hit.evalue < 1e-50:
        # Process this hit
        sequence = retrieve_sequence_from_fasta(
            "genome.fna",
            hit.contig_id,
            hit.start,
            hit.end,
            hit.strand
        )
```

## Testing

```bash
# Run all tests
conda activate opfi
python -m pytest tests/test_parsers.py -v

# Run specific test
python -m pytest tests/test_parsers.py::TestDiamondBlastpParsers -v
```

## Demo

```bash
python examples/demo_parsers.py
```
