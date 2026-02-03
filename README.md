# RNA Guide Editor Finder

A Python bioinformatics toolkit for identifying and analyzing RNA guide and editor sequences in IS110 transposons across bacterial genomes.

## Features

### 1. Sequence Search Module (`modules/sequence_search.py`)
- Search DNA sequences (50-200bp) against GTDB genomic databases using MMseqs2
- Retrieve matching sequences from .fna files
- Deduplicate results based on sequence similarity
- Support for genomic context around hits

### 2. Alignment Filter Module (`modules/alignment_filter.py`)
- Configurable multi-level filtering pipeline:
  - **Composition filters**: AT/GC content, Shannon entropy, homopolymers, dinucleotide repeats
  - **Length filters**: Hard cutoffs and confidence-based filtering
  - **Topology filters**: Boundary artifact detection
- Preset configurations: default, strict, relaxed
- Detailed metrics and filtering results

### 3. Parsers Module (`utils/parsers.py`)
- Comprehensive parsers for bioinformatics file formats
- Automatic sequence retrieval from genomic FASTA files
- Full genomic position tracking

## Installation

### Requirements
- Python 3.8+
- MMseqs2 (for sequence search)
- pytest (for testing)

### Setup
```bash
# Activate conda environment
conda activate opfi

# Install dependencies
pip install pytest

# Run tests
python -m pytest tests/ -v
```

## Parsers Module

The parsers module provides comprehensive support for parsing bioinformatics file formats with automatic position and sequence retrieval.

### Supported Formats

#### 1. Prodigal Output
- **FAA files** (protein sequences)
- **GFF files** (gene annotations)

```python
from utils.parsers import parse_prodigal_faa, parse_prodigal_faa_full

# Parse positions only
position_map = parse_prodigal_faa("proteins.faa")
# Returns: {'contig_1_1': (100, 300, 1), ...}
#           protein_id -> (start, end, strand)

# Parse full information including sequences
genes = parse_prodigal_faa_full("proteins.faa")
# Returns: {'contig_1_1': ProdigalGene(...), ...}
```

#### 2. FASTA Files

```python
from utils.parsers import parse_fasta, parse_fasta_records, iter_fasta

# Parse as dictionary
sequences = parse_fasta("genome.fna")
# Returns: {'contig_1': 'ATCG...', ...}

# Parse as records
records = parse_fasta_records("genome.fna")
# Returns: [FastaRecord(header='...', sequence='...'), ...]

# Memory-efficient iteration
for record in iter_fasta("large_genome.fna"):
    print(f"{record.id}: {record.length} bp")
```

#### 3. GFF/GFF3 Files

```python
from utils.parsers import parse_gff, iter_gff

# Parse entire file
features = parse_gff("annotations.gff")
# Returns: [GffFeature(...), ...]

# Iterate for large files
for feature in iter_gff("large_annotations.gff"):
    if feature.feature_type == "CDS":
        print(f"{feature.seqid}: {feature.start}..{feature.end}")
```

#### 4. BED Files

```python
from utils.parsers import parse_bed, iter_bed

# Parse BED file
records = parse_bed("regions.bed")
# Returns: [BedRecord(...), ...]
```

#### 5. Diamond BLASTP Results

```python
from utils.parsers import (
    parse_prodigal_faa,
    parse_diamond_blastp_with_positions,
    retrieve_sequences_from_hits
)

# Step 1: Get protein positions from Prodigal
position_map = parse_prodigal_faa("proteins.faa")

# Step 2: Parse Diamond BLASTP with positions
hits = parse_diamond_blastp_with_positions(
    "diamond_results.m8",
    position_map,
    file_basename="my_genome"
)

# Each hit contains:
# - file_basename: Source file
# - contig_id: Contig identifier
# - protein_id: Protein identifier
# - start, end: Genomic coordinates (1-based)
# - strand: 1 (forward) or -1 (reverse)
# - query_id, subject_id: Alignment IDs
# - pident, evalue, bitscore: Alignment metrics

# Step 3: Retrieve genomic sequences
hits_with_sequences = retrieve_sequences_from_hits(
    hits,
    ["genome.fna", "other_genome.fna"],
    verbose=True
)

# Access sequences
for hit in hits_with_sequences:
    if hasattr(hit, 'genomic_sequence'):
        print(f"{hit.protein_id}: {hit.genomic_sequence[:50]}...")
```

#### 6. HMMER Results

```python
from utils.parsers import (
    parse_prodigal_faa,
    parse_hmm_tblout_with_positions,
    retrieve_sequences_from_hits
)

# Step 1: Get protein positions from Prodigal
position_map = parse_prodigal_faa("proteins.faa")

# Step 2: Parse HMMER tblout with positions
hits = parse_hmm_tblout_with_positions(
    "hmmer_results.tbl",
    position_map,
    file_basename="my_genome"
)

# Each hit contains:
# - file_basename: Source file
# - contig_id: Contig identifier
# - protein_id: Protein identifier
# - start, end: Genomic coordinates (1-based)
# - strand: 1 (forward) or -1 (reverse)
# - query_name: HMM profile name
# - target_name: Target sequence name
# - evalue, score, bias: HMM metrics

# Step 3: Retrieve genomic sequences
hits_with_sequences = retrieve_sequences_from_hits(
    hits,
    ["genome.fna"],
    verbose=True
)
```

### Sequence Retrieval Functions

#### Direct Sequence Retrieval

```python
from utils.parsers import retrieve_sequence_from_fasta

# Retrieve specific sequence region
sequence = retrieve_sequence_from_fasta(
    "genome.fna",
    contig_id="contig_1",
    start=1000,      # 1-based, inclusive
    end=1100,        # 1-based, inclusive
    strand=1         # 1=forward, -1=reverse (auto rev-comp)
)
```

#### Batch Sequence Retrieval

```python
from utils.parsers import retrieve_sequences_from_hits

# Works with both DiamondBlastHit and HmmHit objects
hits_with_sequences = retrieve_sequences_from_hits(
    hits,
    fasta_files=["genome1.fna", "genome2.fna"],
    verbose=True
)

# Sequences are added as 'genomic_sequence' attribute
for hit in hits_with_sequences:
    if hasattr(hit, 'genomic_sequence'):
        seq = hit.genomic_sequence
        print(f"{hit.protein_id}: {len(seq)} bp")
```

#### Reverse Complement

```python
from utils.parsers import reverse_complement

rev_comp = reverse_complement("ATCGATCG")
# Returns: "CGATCGAT"
```

### Data Classes

#### DiamondBlastHit
```python
@dataclass
class DiamondBlastHit:
    file_basename: str      # Source file identifier
    contig_id: str          # Contig/chromosome identifier
    protein_id: str         # Protein identifier
    start: int              # Genomic start (1-based)
    end: int                # Genomic end (1-based)
    strand: int             # 1=forward, -1=reverse
    query_id: str           # Query sequence ID
    subject_id: str         # Subject sequence ID
    pident: float           # Percentage identity
    alignment_length: int   # Alignment length
    evalue: float           # E-value
    bitscore: float         # Bit score
    # ... additional fields

    # Properties
    @property
    def length(self) -> int:
        """Genomic region length"""

    @property
    def is_forward(self) -> bool:
        """True if on forward strand"""
```

#### HmmHit
```python
@dataclass
class HmmHit:
    file_basename: str      # Source file identifier
    contig_id: str          # Contig/chromosome identifier
    protein_id: str         # Protein identifier
    start: int              # Genomic start (1-based)
    end: int                # Genomic end (1-based)
    strand: int             # 1=forward, -1=reverse
    query_name: str         # HMM profile name
    target_name: str        # Target sequence name
    evalue: float           # E-value
    score: float            # Score
    bias: float             # Bias
    # ... additional domain fields

    # Properties
    @property
    def length(self) -> int:
        """Genomic region length"""

    @property
    def is_forward(self) -> bool:
        """True if on forward strand"""
```

## Complete Workflow Examples

### Example 1: Diamond BLASTP Analysis

```python
from utils.parsers import (
    parse_prodigal_faa,
    parse_diamond_blastp_with_positions,
    retrieve_sequences_from_hits
)

# 1. Parse Prodigal output
position_map = parse_prodigal_faa("genome.faa")

# 2. Parse Diamond BLASTP results
hits = parse_diamond_blastp_with_positions(
    "diamond_results.m8",
    position_map,
    file_basename="my_genome"
)

# 3. Filter hits by E-value
good_hits = [h for h in hits if h.evalue < 1e-50]

# 4. Retrieve sequences
hits_with_seq = retrieve_sequences_from_hits(
    good_hits,
    ["genome.fna"],
    verbose=True
)

# 5. Analyze results
for hit in hits_with_seq:
    print(f"{hit.subject_id}: {hit.contig_id}:{hit.start}-{hit.end}")
    print(f"  E-value: {hit.evalue}, Identity: {hit.pident}%")
    if hasattr(hit, 'genomic_sequence'):
        print(f"  Sequence: {hit.genomic_sequence[:50]}...")
```

### Example 2: HMMER Analysis

```python
from utils.parsers import (
    parse_prodigal_faa,
    parse_hmm_tblout_with_positions,
    retrieve_sequences_from_hits
)

# 1. Parse Prodigal output
position_map = parse_prodigal_faa("genome.faa")

# 2. Parse HMMER results
hits = parse_hmm_tblout_with_positions(
    "hmmer_results.tbl",
    position_map,
    file_basename="my_genome"
)

# 3. Filter by score
good_hits = [h for h in hits if h.score > 100]

# 4. Group by HMM profile
from collections import defaultdict
by_profile = defaultdict(list)
for hit in good_hits:
    by_profile[hit.query_name].append(hit)

# 5. Retrieve sequences
hits_with_seq = retrieve_sequences_from_hits(
    good_hits,
    ["genome.fna"],
    verbose=True
)

# 6. Analyze results
for profile, profile_hits in by_profile.items():
    print(f"{profile}: {len(profile_hits)} hits")
```

### Example 3: Memory-Efficient Processing

```python
from utils.parsers import (
    parse_prodigal_faa,
    iter_diamond_blastp_with_positions,
    retrieve_sequence_from_fasta
)

# Parse positions
position_map = parse_prodigal_faa("genome.faa")

# Process hits one at a time (memory efficient)
for hit in iter_diamond_blastp_with_positions("diamond_results.m8", position_map):
    if hit.evalue < 1e-50:
        # Retrieve sequence for this specific hit
        sequence = retrieve_sequence_from_fasta(
            "genome.fna",
            hit.contig_id,
            hit.start,
            hit.end,
            hit.strand
        )

        if sequence:
            print(f"{hit.subject_id}: {len(sequence)} bp")
```

## Testing

Run the comprehensive test suite:

```bash
# Activate environment
conda activate opfi

# Run all tests
python -m pytest tests/ -v

# Run specific test file
python -m pytest tests/test_parsers.py -v

# Run specific test
python -m pytest tests/test_parsers.py::TestDiamondBlastpParsers::test_parse_diamond_blastp_with_positions -v
```

## Demo Scripts

See example usage:

```bash
python examples/demo_parsers.py
```

## Project Structure

```
RNA_guide_editor_finder/
├── modules/
│   ├── __init__.py
│   ├── sequence_search.py      # MMseqs2 search functionality
│   └── alignment_filter.py     # Alignment filtering pipeline
├── utils/
│   ├── __init__.py
│   └── parsers.py             # Comprehensive parsers module
├── tests/
│   ├── test_parsers.py        # Parser unit tests
│   └── test_sequence_search_integration.py
├── examples/
│   └── demo_parsers.py        # Demo scripts
└── README.md                  # This file
```

## Contributing

When adding new parsers:
1. Add parser functions to `utils/parsers.py`
2. Create corresponding data classes if needed
3. Add unit tests to `tests/test_parsers.py`
4. Update `utils/__init__.py` to export new functions
5. Update this README with usage examples

## License

This project is for research use in studying IS110 transposons.

## Contact

For questions or issues, please contact the project maintainers.
