# Understanding Coordinates in region_extractor

## What do start and end mean?

**start** and **end** refer to the **transposase gene's coding region** - the actual protein-coding sequence that gets translated.

## Visual Example

### Forward Strand (strand=1)

```
Genomic DNA (5' to 3'):
Position:  1    ...   650  ...  999 1000 ... 2000 2001 ... 2250  ... 3000

           [    genome    ][upstream][===== GENE =====][downstream][  genome  ]
                              350bp      1001bp CDS        250bp

Parameters:
  start = 1000  ← Start of transposase gene
  end = 2000    ← End of transposase gene
  upstream_length = 350
  downstream_length = 250

Extracted Regions:
  upstream:    positions 650-999   (350 bp before gene)
  coding:      positions 1000-2000 (1001 bp, the transposase gene itself)
  downstream:  positions 2001-2250 (250 bp after gene)
```

### Reverse Strand (strand=-1)

```
Genomic DNA (5' to 3'):
Position:  1    ...   650  ...  999 1000 ... 2000 2001 ... 2250  ... 3000

           [    genome    ][downstream][===== GENE =====][upstream][  genome  ]
                              250bp      1001bp CDS        350bp

Parameters:
  start = 1000  ← Gene start (still 1000, same as genomic position)
  end = 2000    ← Gene end
  strand = -1   ← Reverse strand
  upstream_length = 350
  downstream_length = 250

Extracted Regions (in gene orientation):
  upstream:    positions 2001-2350 (genomically downstream, but gene's upstream)
  coding:      positions 1000-2000 (reverse complemented)
  downstream:  positions 650-999   (genomically upstream, but gene's downstream)

Note: All sequences are reverse complemented to match gene orientation!
```

## Real Example: IS110 Transposon

```
IS110 Transposon Structure:

[IR-L]  [promoter]  [======= transposase gene =======]  [terminator]  [IR-R]
  ↑                  ↑                                ↑                   ↑
  │                  start=X                         end=Y               │
  │                  (protein starts)                (protein ends)      │
  │                                                                       │
  └─── upstream region ───┘                          └─── downstream ────┘
       (may contain IR-L,                                 (may contain IR-R,
        promoter, etc.)                                    terminator, etc.)
```

When you call:
```python
result = extractor.extract_protein_regions(
    contig_id="...",
    start=X,        # ← Transposase gene start
    end=Y,          # ← Transposase gene end
    strand=1,
    upstream_length=350,
    downstream_length=250
)
```

You get:
- **coding_sequence**: The transposase protein-coding sequence (what gets translated)
- **upstream_sequence**: 350 bp before the gene (promoter region, may have IR-L)
- **downstream_sequence**: 250 bp after the gene (terminator region, may have IR-R)

## Common Use Cases

### 1. Extract full transposon with flanking regions
```python
# Get transposase + 500bp on each side
result = extractor.extract_protein_regions(
    contig_id="NC_000001",
    start=10000,      # Transposase starts here
    end=11200,        # Transposase ends here
    strand=1,
    upstream_length=500,
    downstream_length=500
)
```

### 2. Look for promoters
```python
# Get 300bp upstream to search for promoter motifs
result = extractor.extract_protein_regions(
    start=gene_start,
    end=gene_end,
    upstream_length=300,  # Search for promoters here
    downstream_length=50
)

promoter_region = result['upstream_sequence']
# Search for -10, -35 boxes, etc.
```

### 3. Find inverted repeats (IRs)
```python
# IS elements often have IRs flanking the transposase
result = extractor.extract_protein_regions(
    start=gene_start,
    end=gene_end,
    upstream_length=150,    # IR-L might be here
    downstream_length=150   # IR-R might be here
)

left_flank = result['upstream_sequence'][-50:]   # Last 50bp of upstream
right_flank = result['downstream_sequence'][:50]  # First 50bp of downstream
# Check if they're inverted repeats
```

## Summary

| Parameter | Meaning |
|-----------|---------|
| `start` | Where the transposase gene **starts** |
| `end` | Where the transposase gene **ends** |
| `upstream_length` | How many bp before the gene to extract |
| `downstream_length` | How many bp after the gene to extract |

The **coding region** (start to end) is the transposase itself. The **flanking regions** are the non-coding sequences around it.
