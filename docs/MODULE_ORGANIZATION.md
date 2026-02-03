# Module Organization

## Overview

The RNA Guide Editor Finder project is organized into reusable modules that avoid code duplication.

```
RNA_guide_editor_finder/
├── modules/          # High-level feature modules
│   ├── region_extractor.py    # Extract genomic regions
│   ├── sequence_search.py     # Sequence searching
│   └── alignment_filter.py    # Filter alignments
├── utils/            # Low-level utility functions
│   └── parsers.py             # File format parsers
└── tests/            # Unit tests
```

## Module Dependencies

### modules/region_extractor.py

**Purpose**: Extract coding regions and flanking sequences from genomic data.

**Dependencies**:
- `utils.parsers` → imports `parse_prodigal_faa()`, `reverse_complement()`
- `Bio.SeqIO` → for FASTA file reading

**Exports**:
- `RegionExtractor` class
- `extract_contig_from_protein_id()` function
- `create_transposon_dict()` convenience function
- Re-exports: `parse_prodigal_faa()`, `reverse_complement()`

### utils/parsers.py

**Purpose**: Low-level parsers for bioinformatics file formats.

**Key Functions**:
- `parse_prodigal_faa()` - Parse Prodigal FAA files → (start, end, strand)
- `parse_prodigal_faa_full()` - Parse with full gene info + sequences
- `parse_prodigal_gff()` - Parse Prodigal GFF files
- `parse_fasta()` - Parse FASTA files
- `parse_gff()` - Parse GFF/GFF3 files
- `parse_bed()` - Parse BED files
- `parse_diamond_blastp()` - Parse DIAMOND BLAST output
- `parse_hmm_tblout()` - Parse HMMER tblout files
- `reverse_complement()` - DNA reverse complement
- `retrieve_sequence_from_fasta()` - Extract sequences by coordinates

**Data Classes**:
- `ProdigalGene`
- `FastaRecord`
- `GffFeature`
- `BedRecord`
- `DiamondBlastHit`
- `HmmHit`

## Design Principles

### 1. No Code Duplication

Instead of duplicating `parse_prodigal_faa()` and `reverse_complement()` in multiple modules, we:
- Define them once in `utils.parsers`
- Import them where needed
- Re-export from `modules` for convenience

### 2. Layered Architecture

```
User Scripts
     ↓
modules/          (High-level features)
     ↓
utils/            (Low-level utilities)
     ↓
External Libraries (Biopython, etc.)
```

### 3. Clear Responsibilities

**utils/** = Generic, reusable utilities
- Can be used in any project
- No domain-specific logic
- Pure functions when possible

**modules/** = Domain-specific features
- Transposon-specific logic
- Combines multiple utilities
- Provides convenient APIs

## Usage Examples

### Using region_extractor

```python
# Import from modules (recommended)
from modules.region_extractor import RegionExtractor, create_transposon_dict

# Or import everything from modules
from modules import (
    RegionExtractor,
    parse_prodigal_faa,
    reverse_complement
)
```

### Using parsers directly

```python
# For direct access to parsing utilities
from utils.parsers import (
    parse_prodigal_faa,
    parse_fasta,
    parse_gff,
    reverse_complement
)
```

## Why This Organization?

### Benefits

1. **Maintainability**: Fix a bug once, it's fixed everywhere
2. **Testability**: Test utilities independently of higher-level features
3. **Reusability**: Use parsers in multiple modules without duplication
4. **Clarity**: Clear separation between generic utilities and domain logic
5. **Flexibility**: Can use parsers standalone or through higher-level modules

### Example: parse_prodigal_faa

**Before refactoring**:
- Duplicate code in `region_extractor.py` and `parsers.py`
- Bug fixes need to be applied twice
- Tests needed in multiple places

**After refactoring**:
- Single source of truth in `utils.parsers`
- Imported by `region_extractor.py`
- One set of tests
- Fixes propagate automatically

## Best Practices

### When to add to utils/

Add functions to `utils/` when they are:
- Generic file format parsers
- Reusable across projects
- Have no domain-specific logic

### When to add to modules/

Add code to `modules/` when it:
- Combines multiple utilities
- Has transposon/IS110-specific logic
- Provides high-level workflows

### When to create a new module

Create a new module when:
- It has a distinct, cohesive purpose
- It will be >200 lines of code
- Multiple scripts will use it
- It deserves its own documentation

## Summary

The `region_extractor` module provides high-level functionality for extracting transposon regions, while relying on battle-tested utilities from `utils.parsers` for low-level file parsing. This design eliminates code duplication and makes the codebase more maintainable.
