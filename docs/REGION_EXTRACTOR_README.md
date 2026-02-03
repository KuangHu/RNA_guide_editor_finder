# Region Extractor Module

A reusable Python module for extracting coding regions and flanking sequences from genomic data.

## Overview

The `region_extractor` module provides tools to extract:
- **Coding regions** (gene sequences) - the transposase/protein-coding sequence
- **Upstream flanking regions** (non-coding sequences before the gene)
- **Downstream flanking regions** (non-coding sequences after the gene)

All sequences are properly oriented according to gene strand direction, with reverse complement applied for reverse-strand genes.

### What the coordinates mean

```
Genomic DNA:

[upstream region] [====== CODING REGION ======] [downstream region]
      650-999          1000-2000                    2001-2250
     (350 bp)      (transposase gene)               (250 bp)
                   ↑                  ↑
                 start              end
              (parameter)        (parameter)
```

- **start, end** = The transposase gene coordinates (protein coding sequence)
- **upstream** = Non-coding region before the gene (may contain promoters, IRs)
- **downstream** = Non-coding region after the gene (may contain terminators, IRs)

## Key Features

- ✓ Handles both forward and reverse strand orientations
- ✓ Automatic reverse complementation for reverse strand genes
- ✓ 1-based coordinate system (compatible with Prodigal/GFF)
- ✓ Boundary-aware extraction (handles contig edges)
- ✓ Integration with Prodigal FAA file parsing
- ✓ Clean dictionary output format for downstream processing
- ✓ Integrates with existing `utils.parsers` module (no code duplication)

## Quick Start

### Basic Usage

```python
from modules.region_extractor import RegionExtractor

# Initialize with genome file
extractor = RegionExtractor("/path/to/genome.fna")

# Extract regions for a protein
result = extractor.extract_protein_regions(
    contig_id="AAXX01000001.1",
    start=1000,
    end=2000,
    strand=1,  # 1 for forward, -1 for reverse
    upstream_length=350,
    downstream_length=250
)

# Access sequences
coding_seq = result['coding_sequence']
upstream_seq = result['upstream_sequence']
downstream_seq = result['downstream_sequence']

# Access coordinates
print(result['coding_coords'])  # {'start': 1000, 'end': 2000, 'length': 1001}
```

### Convenience Function

For quick one-off extractions:

```python
from modules.region_extractor import create_transposon_dict

transposon = create_transposon_dict(
    genome_fna_path="/path/to/genome.fna",
    contig_id="AAXX01000001.1",
    protein_id="AAXX01000001.1_123",
    start=1000,
    end=2000,
    strand=1,
    upstream_length=350,
    downstream_length=250,
    additional_metadata={'annotation': 'IS110 transposase'}
)
```

### Integration with Prodigal

```python
from modules.region_extractor import (
    RegionExtractor,
    parse_prodigal_faa,
    extract_contig_from_protein_id
)

# Parse Prodigal FAA file
positions = parse_prodigal_faa("/path/to/proteins.faa")

# Initialize extractor
extractor = RegionExtractor("/path/to/genome.fna")

# Process a specific protein
protein_id = "AAXX01000001.1_123"
if protein_id in positions:
    start, end, strand = positions[protein_id]
    contig_id = extract_contig_from_protein_id(protein_id)

    result = extractor.extract_protein_regions(
        contig_id=contig_id,
        start=start,
        end=end,
        strand=strand
    )
```

## Output Format

The `extract_protein_regions` method returns a dictionary with:

```python
{
    'contig_id': 'AAXX01000001.1',
    'strand': 1,
    'coding_sequence': 'ATGCGT...',
    'upstream_sequence': 'TTAACG...',
    'downstream_sequence': 'CGATAA...',
    'coding_coords': {
        'start': 1000,
        'end': 2000,
        'length': 1001
    },
    'upstream_coords': {
        'start': 650,
        'end': 999,
        'length': 350
    },
    'downstream_coords': {
        'start': 2001,
        'end': 2250,
        'length': 250
    }
}
```

## API Reference

### RegionExtractor Class

#### `__init__(genome_fna_path=None)`
Initialize the extractor, optionally loading a genome file.

#### `load_genome(fna_path)`
Load genomic sequences from a FASTA file.

#### `extract_protein_regions(contig_id, start, end, strand, upstream_length=350, downstream_length=250)`
Extract coding and flanking regions for a gene.

**Parameters:**
- `contig_id` (str): Contig/chromosome identifier
- `start` (int): Gene start position (1-based, inclusive)
- `end` (int): Gene end position (1-based, inclusive)
- `strand` (int): 1 for forward, -1 for reverse
- `upstream_length` (int): Length of upstream region to extract
- `downstream_length` (int): Length of downstream region to extract

**Returns:** Dictionary with sequences and coordinates

### Utility Functions

#### `reverse_complement(seq)`
Return reverse complement of a DNA sequence.
*Note: Imported from `utils.parsers` module*

#### `parse_prodigal_faa(faa_path)`
Parse Prodigal FAA file to extract gene positions.
Returns: `{protein_id: (start, end, strand)}`
*Note: Imported from `utils.parsers` module*

#### `extract_contig_from_protein_id(protein_id)`
Extract contig ID from protein ID.
Example: `'AAXX01000001.1_281'` → `'AAXX01000001.1'`

#### `create_transposon_dict(...)`
Convenience function to create a complete transposon dictionary entry.

## Complete Workflow Example

Here's how to process multiple proteins from a CSV file (similar to the GTDB1 workflow):

```python
import csv
import json
from collections import defaultdict
from modules.region_extractor import (
    RegionExtractor,
    parse_prodigal_faa,
    extract_contig_from_protein_id
)

# 1. Read input CSV
proteins_data = []
with open('proteins.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        proteins_data.append(row)

# 2. Group by genome to minimize file I/O
proteins_by_genome = defaultdict(list)
for row in proteins_data:
    genome = row['genome_name']
    proteins_by_genome[genome].append(row)

# 3. Process each genome
all_results = []
for genome_name, protein_rows in proteins_by_genome.items():
    # Load genome and parse Prodigal annotations
    fna_path = f"/path/to/genomes/{genome_name}.fna"
    faa_path = f"/path/to/proteins/{genome_name}.faa"

    extractor = RegionExtractor(fna_path)
    positions = parse_prodigal_faa(faa_path)

    # Process each protein
    for row in protein_rows:
        protein_id = row['protein_id']

        if protein_id not in positions:
            continue

        start, end, strand = positions[protein_id]
        contig_id = extract_contig_from_protein_id(protein_id)

        # Extract regions
        result = extractor.extract_protein_regions(
            contig_id=contig_id,
            start=start,
            end=end,
            strand=strand,
            upstream_length=350,
            downstream_length=250
        )

        # Build complete dictionary
        transposon = {
            'protein_id': protein_id,
            'genome': genome_name,
            'contig_id': contig_id,
            'strand': strand,
            'coordinates': {
                'coding': result['coding_coords'],
                'upstream': result['upstream_coords'],
                'downstream': result['downstream_coords']
            },
            'sequences': {
                'coding': result['coding_sequence'],
                'upstream': result['upstream_sequence'],
                'downstream': result['downstream_sequence']
            },
            # Add any additional metadata from CSV
            'domains': {
                'DEDD': {
                    'evalue': row.get('DEDD_evalue'),
                    'start': row.get('DEDD_start'),
                    'end': row.get('DEDD_end')
                }
            }
        }

        all_results.append(transposon)

# 4. Save results
with open('output.json', 'w') as f:
    json.dump(all_results, f, indent=2)

print(f"Saved {len(all_results)} transposons")
```

## Coordinate System

**Important:** This module uses **1-based, inclusive** coordinates (Prodigal/GFF format):
- `start=100, end=200` includes positions 100, 101, ..., 200 (101 bases total)
- Length = end - start + 1

## Strand Orientation

### Forward Strand (strand=1)
```
Genomic:  [upstream] [====gene====] [downstream]
          ^                          ^
          start-350                  end+250
```

### Reverse Strand (strand=-1)
```
Genomic:  [downstream] [====gene====] [upstream]
          ^                            ^
          start-250                    end+350

Note: Sequences are reverse complemented
```

## Testing

Run the manual test:
```bash
python3 test_manual.py
```

Or use pytest (if installed):
```bash
pytest tests/test_region_extractor.py -v
```

## Files

- `modules/region_extractor.py` - Main module
- `examples/extract_regions_example.py` - Usage examples
- `tests/test_region_extractor.py` - Unit tests
- `test_manual.py` - Simple manual tests

## Integration with Existing Code

To integrate with your existing GTDB1 script:

1. Replace the standalone functions with imports from this module
2. Use `RegionExtractor` class instead of separate `load_genome_sequence` calls
3. Use `parse_prodigal_faa` and `extract_contig_from_protein_id` helper functions

The module is drop-in compatible with your existing workflow!
