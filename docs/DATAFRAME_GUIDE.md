# DataFrame Analysis Guide

Complete guide for converting parser results to pandas DataFrames for analysis.

## Installation

```bash
# Install pandas (if not already installed)
conda activate opfi
pip install pandas

# For Excel export (optional)
pip install openpyxl
```

## Quick Start

```python
from utils.parsers import (
    parse_prodigal_faa,
    parse_diamond_blastp_with_positions,
    diamond_hits_to_dataframe,
    retrieve_sequences_from_hits,
    save_hits_to_csv
)

# 1. Parse Diamond BLASTP results
position_map = parse_prodigal_faa("genome.faa")
hits = parse_diamond_blastp_with_positions(
    "diamond_results.m8",
    position_map,
    file_basename="my_genome"
)

# 2. Retrieve sequences (optional)
hits = retrieve_sequences_from_hits(hits, ["genome.fna"])

# 3. Convert to DataFrame
df = diamond_hits_to_dataframe(hits)

# 4. Analyze!
print(df.head())
print(df.describe())
```

## DataFrame Columns

### Diamond BLASTP DataFrame

| Column | Type | Description |
|--------|------|-------------|
| `file_basename` | str | Source genome file |
| `contig_id` | str | Contig/chromosome identifier |
| `protein_id` | str | Protein identifier |
| `start` | int | Genomic start position (1-based) |
| `end` | int | Genomic end position (1-based) |
| `strand` | int | 1 (forward) or -1 (reverse) |
| `query_id` | str | Query sequence ID |
| `subject_id` | str | Subject sequence ID (what it matched) |
| `pident` | float | Percentage identity |
| `alignment_length` | int | Alignment length |
| `mismatches` | int | Number of mismatches |
| `gap_opens` | int | Number of gap openings |
| `query_start` | int | Query alignment start |
| `query_end` | int | Query alignment end |
| `subject_start` | int | Subject alignment start |
| `subject_end` | int | Subject alignment end |
| `evalue` | float | E-value |
| `bitscore` | float | Bit score |
| `length` | int | Genomic region length |
| `is_forward` | bool | True if forward strand |
| `genomic_sequence` | str | DNA sequence (if retrieved) |

### HMMER DataFrame

| Column | Type | Description |
|--------|------|-------------|
| `file_basename` | str | Source genome file |
| `contig_id` | str | Contig/chromosome identifier |
| `protein_id` | str | Protein identifier |
| `start` | int | Genomic start position (1-based) |
| `end` | int | Genomic end position (1-based) |
| `strand` | int | 1 (forward) or -1 (reverse) |
| `query_name` | str | HMM profile name |
| `target_name` | str | Target sequence name |
| `evalue` | float | E-value |
| `score` | float | Score |
| `bias` | float | Bias |
| `domain_evalue` | float | Best domain E-value |
| `domain_score` | float | Best domain score |
| `domain_number` | int | Number of domains |
| `length` | int | Genomic region length |
| `is_forward` | bool | True if forward strand |
| `genomic_sequence` | str | DNA sequence (if retrieved) |

## Common Analysis Tasks

### 1. Basic Information

```python
# View first few rows
df.head()

# View last few rows
df.tail()

# Get DataFrame info
df.info()

# Summary statistics
df.describe()

# Count rows
len(df)

# Column names
df.columns.tolist()
```

### 2. Filtering

```python
# Filter by E-value
high_confidence = df[df['evalue'] < 1e-50]

# Filter by identity
high_identity = df[df['pident'] > 90]

# Filter by strand
forward_only = df[df['is_forward']]
reverse_only = df[~df['is_forward']]

# Filter by genome
genome1_hits = df[df['file_basename'] == 'genome1']

# Multiple conditions (AND)
filtered = df[
    (df['evalue'] < 1e-50) &
    (df['pident'] > 90) &
    (df['is_forward'])
]

# Multiple conditions (OR)
filtered = df[
    (df['evalue'] < 1e-100) |
    (df['pident'] > 95)
]

# Filter by subject ID
transposase_hits = df[df['subject_id'] == 'IS110_transposase']

# Filter by subject ID pattern
transposase_hits = df[df['subject_id'].str.contains('transposase')]

# Filter by length
long_hits = df[df['length'] >= 200]
```

### 3. Sorting

```python
# Sort by E-value (ascending)
df_sorted = df.sort_values('evalue')

# Sort by identity (descending)
df_sorted = df.sort_values('pident', ascending=False)

# Sort by multiple columns
df_sorted = df.sort_values(['subject_id', 'evalue'])

# Sort and get top 10
top10 = df.sort_values('bitscore', ascending=False).head(10)
```

### 4. Grouping and Aggregation

```python
# Count hits per genome
df.groupby('file_basename').size()

# Average identity per subject
df.groupby('subject_id')['pident'].mean()

# Multiple aggregations
summary = df.groupby('subject_id').agg({
    'protein_id': 'count',
    'pident': ['mean', 'min', 'max'],
    'evalue': 'min',
    'bitscore': 'max'
})

# Group by multiple columns
summary = df.groupby(['file_basename', 'subject_id']).agg({
    'protein_id': 'count',
    'pident': 'mean'
})

# Count by strand
df.groupby('is_forward').size()

# Group and filter
good_subjects = df.groupby('subject_id').filter(
    lambda x: (x['pident'].mean() > 90) and (len(x) >= 3)
)
```

### 5. Statistics

```python
# Mean identity
df['pident'].mean()

# Median E-value
df['evalue'].median()

# Standard deviation of bitscore
df['bitscore'].std()

# Min and max
df['length'].min()
df['length'].max()

# Quantiles
df['pident'].quantile([0.25, 0.5, 0.75])

# Value counts
df['subject_id'].value_counts()

# Correlation between columns
df[['pident', 'bitscore', 'evalue']].corr()
```

### 6. Adding Custom Columns

```python
# Add log E-value
df['log_evalue'] = -df['evalue'].apply(lambda x: np.log10(x) if x > 0 else 300)

# Add identity category
df['identity_category'] = pd.cut(
    df['pident'],
    bins=[0, 80, 90, 95, 100],
    labels=['low', 'medium', 'high', 'very_high']
)

# Add length category
df['length_category'] = df['length'].apply(
    lambda x: 'short' if x < 150 else 'medium' if x < 250 else 'long'
)

# Add genome + contig combined
df['genome_contig'] = df['file_basename'] + ':' + df['contig_id']

# Check if sequence retrieved
df['has_sequence'] = df['genomic_sequence'].notna()
```

### 7. Selection and Projection

```python
# Select specific columns
subset = df[['protein_id', 'subject_id', 'evalue', 'pident']]

# Select columns by pattern
numeric_cols = df.select_dtypes(include=['float', 'int'])

# Rename columns
df_renamed = df.rename(columns={
    'pident': 'percent_identity',
    'evalue': 'e_value'
})

# Drop columns
df_slim = df.drop(columns=['genomic_sequence'])
```

### 8. Finding Best Hits

```python
# Best hit per genome (by bitscore)
best_hits = df.loc[df.groupby('file_basename')['bitscore'].idxmax()]

# Best hit per subject (by E-value)
best_per_subject = df.loc[df.groupby('subject_id')['evalue'].idxmin()]

# Top 5 hits per genome
top_per_genome = df.groupby('file_basename').apply(
    lambda x: x.nsmallest(5, 'evalue')
)
```

### 9. Comparison Analysis

```python
# Compare two genomes
genome1 = df[df['file_basename'] == 'genome1']
genome2 = df[df['file_basename'] == 'genome2']

print(f"Genome1: {len(genome1)} hits, mean identity: {genome1['pident'].mean():.2f}%")
print(f"Genome2: {len(genome2)} hits, mean identity: {genome2['pident'].mean():.2f}%")

# Compare subject distributions
import matplotlib.pyplot as plt

df.groupby(['file_basename', 'subject_id']).size().unstack().plot(kind='bar')
plt.title('Hit Distribution by Genome and Subject')
plt.ylabel('Number of Hits')
plt.show()
```

### 10. Export Results

```python
# Export to CSV
df.to_csv('results.csv', index=False)

# Export to CSV (selected columns only)
df[['protein_id', 'subject_id', 'evalue', 'pident']].to_csv(
    'summary.csv',
    index=False
)

# Export to Excel
df.to_excel('results.xlsx', index=False, sheet_name='Diamond_Hits')

# Export multiple DataFrames to Excel
with pd.ExcelWriter('results.xlsx') as writer:
    df_diamond.to_excel(writer, sheet_name='Diamond', index=False)
    df_hmmer.to_excel(writer, sheet_name='HMMER', index=False)

# Export filtered results
high_confidence = df[df['evalue'] < 1e-50]
high_confidence.to_csv('high_confidence.csv', index=False)

# Use built-in export functions
from utils.parsers import save_hits_to_csv, save_hits_to_excel

save_hits_to_csv(hits, 'diamond_results.csv')
save_hits_to_excel(hits, 'diamond_results.xlsx')
```

## Complete Workflow Examples

### Example 1: Find High-Confidence Transposase Hits

```python
from utils.parsers import (
    parse_prodigal_faa,
    parse_diamond_blastp_with_positions,
    diamond_hits_to_dataframe,
    save_hits_to_csv
)

# Parse results
position_map = parse_prodigal_faa("genome.faa")
hits = parse_diamond_blastp_with_positions(
    "diamond_results.m8",
    position_map,
    file_basename="my_genome"
)

# Convert to DataFrame
df = diamond_hits_to_dataframe(hits)

# Filter for high-confidence transposase hits
transposase = df[
    (df['subject_id'].str.contains('transposase', case=False)) &
    (df['evalue'] < 1e-50) &
    (df['pident'] > 90)
]

print(f"Found {len(transposase)} high-confidence transposase hits")
print(transposase[['protein_id', 'subject_id', 'pident', 'evalue']])

# Export
transposase.to_csv('transposase_hits.csv', index=False)
```

### Example 2: Compare Multiple Genomes

```python
import pandas as pd
from utils.parsers import (
    parse_prodigal_faa,
    parse_diamond_blastp_with_positions,
    diamond_hits_to_dataframe
)

# Parse results from multiple genomes
all_hits = []

for genome_file in ['genome1.faa', 'genome2.faa', 'genome3.faa']:
    position_map = parse_prodigal_faa(genome_file)
    hits = parse_diamond_blastp_with_positions(
        f"{genome_file.replace('.faa', '_diamond.m8')}",
        position_map,
        file_basename=genome_file.replace('.faa', '')
    )
    all_hits.extend(hits)

# Convert to DataFrame
df = diamond_hits_to_dataframe(all_hits)

# Compare genomes
comparison = df.groupby('file_basename').agg({
    'protein_id': 'count',
    'pident': 'mean',
    'evalue': 'min',
    'bitscore': 'max'
})
comparison.columns = ['num_hits', 'mean_identity', 'best_evalue', 'max_bitscore']

print(comparison)

# Find subjects present in all genomes
subjects_per_genome = df.groupby(['file_basename', 'subject_id']).size().unstack(fill_value=0)
common_subjects = subjects_per_genome.loc[:, (subjects_per_genome > 0).all()]

print(f"\nSubjects found in all genomes: {len(common_subjects.columns)}")
print(common_subjects.columns.tolist())
```

### Example 3: Quality Control Report

```python
from utils.parsers import (
    parse_diamond_blastp_with_positions,
    parse_prodigal_faa,
    diamond_hits_to_dataframe
)
import pandas as pd

# Parse and convert
position_map = parse_prodigal_faa("genome.faa")
hits = parse_diamond_blastp_with_positions(
    "diamond_results.m8",
    position_map,
    file_basename="my_genome"
)
df = diamond_hits_to_dataframe(hits)

# Generate QC report
print("=" * 70)
print("QUALITY CONTROL REPORT")
print("=" * 70)

print(f"\nTotal hits: {len(df)}")
print(f"Unique subjects: {df['subject_id'].nunique()}")
print(f"Unique proteins: {df['protein_id'].nunique()}")

print("\nE-value distribution:")
print(f"  < 1e-100: {len(df[df['evalue'] < 1e-100])}")
print(f"  < 1e-50: {len(df[df['evalue'] < 1e-50])}")
print(f"  < 1e-10: {len(df[df['evalue'] < 1e-10])}")

print("\nIdentity distribution:")
print(f"  > 95%: {len(df[df['pident'] > 95])}")
print(f"  > 90%: {len(df[df['pident'] > 90])}")
print(f"  > 80%: {len(df[df['pident'] > 80])}")

print("\nLength distribution:")
print(df['length'].describe())

print("\nTop 10 subjects by hit count:")
print(df['subject_id'].value_counts().head(10))

print("\nStrand distribution:")
print(f"  Forward: {len(df[df['is_forward']])}")
print(f"  Reverse: {len(df[~df['is_forward']])}")
```

## Tips and Best Practices

1. **Memory Management:**
   - For large datasets, exclude `genomic_sequence` column when exporting
   - Use iterators for very large files (see `iter_diamond_blastp_with_positions`)

2. **Performance:**
   - Filter early to reduce DataFrame size
   - Use vectorized operations instead of loops
   - Convert to appropriate dtypes (e.g., int8 for strand)

3. **Data Quality:**
   - Check for missing values: `df.isna().sum()`
   - Check for duplicates: `df.duplicated().sum()`
   - Validate ranges: `df['pident'].between(0, 100).all()`

4. **Visualization:**
   - Use matplotlib or seaborn for plots
   - Common plots: histograms, scatter plots, bar charts

5. **Export Format:**
   - CSV: Best for large datasets, universal compatibility
   - Excel: Good for small-medium datasets, nice formatting
   - Parquet: Best for very large datasets (requires pyarrow)

## Demo Script

Run the comprehensive demo:

```bash
conda activate opfi
python examples/demo_dataframe_analysis.py
```

This shows examples of:
- Converting hits to DataFrames
- Filtering and sorting
- Grouping and aggregation
- Cross-genome comparisons
- Export to CSV/Excel
