# Complete Workflow: From Raw Files to DataFrame Analysis

```
┌────────────────────────────────────────────────────────────────────────┐
│                         INPUT FILES                                     │
└────────────────────────────────────────────────────────────────────────┘

  genome.faa              diamond_results.m8         genome.fna
  (Prodigal)              (Diamond BLASTP)           (FASTA sequences)
      │                          │                          │
      │                          │                          │
      ▼                          ▼                          │
┌─────────────┐          ┌──────────────┐                  │
│Parse        │          │Parse Diamond │                  │
│Prodigal FAA │          │BLASTP        │                  │
│             │          │              │                  │
│Returns:     │          │Needs:        │                  │
│position_map │─────────▶│position_map  │                  │
│             │          │              │                  │
│             │          │Returns:      │                  │
│             │          │List[Diamond  │                  │
│             │          │BlastHit]     │                  │
└─────────────┘          └──────────────┘                  │
                                │                          │
                                │ Each hit has:            │
                                │ - file_basename          │
                                │ - contig_id             │
                                │ - protein_id            │
                                │ - start, end            │
                                │ - strand                │
                                │ - subject_id            │
                                │ - pident, evalue        │
                                │                          │
                                ▼                          │
                         ┌──────────────┐                  │
                         │Retrieve      │◀─────────────────┘
                         │Sequences     │
                         │              │
                         │Adds:         │
                         │genomic_      │
                         │sequence      │
                         └──────────────┘
                                │
                                ▼
                         ┌──────────────┐
                         │Convert to    │
                         │DataFrame     │
                         │              │
                         │Returns:      │
                         │pandas        │
                         │DataFrame     │
                         └──────────────┘
                                │
                ┌───────────────┼───────────────┐
                │               │               │
                ▼               ▼               ▼
         ┌───────────┐   ┌───────────┐  ┌───────────┐
         │Filter &   │   │Group &    │  │Export     │
         │Sort       │   │Aggregate  │  │Results    │
         │           │   │           │  │           │
         │df[df[...]]│   │.groupby() │  │.to_csv()  │
         │.sort_     │   │.agg()     │  │.to_excel()│
         │values()   │   │           │  │           │
         └───────────┘   └───────────┘  └───────────┘

┌────────────────────────────────────────────────────────────────────────┐
│                         EXAMPLE CODE                                    │
└────────────────────────────────────────────────────────────────────────┘

from utils.parsers import (
    parse_prodigal_faa,
    parse_diamond_blastp_with_positions,
    retrieve_sequences_from_hits,
    diamond_hits_to_dataframe
)

# Step 1: Parse Prodigal
position_map = parse_prodigal_faa("genome.faa")
# {'contig_1_1': (100, 300, 1), 'contig_1_2': (500, 800, -1), ...}

# Step 2: Parse Diamond with positions
hits = parse_diamond_blastp_with_positions(
    "diamond_results.m8",
    position_map,
    file_basename="my_genome"
)
# [DiamondBlastHit(...), DiamondBlastHit(...), ...]

# Step 3: Retrieve sequences
hits = retrieve_sequences_from_hits(hits, ["genome.fna"])
# Adds genomic_sequence to each hit

# Step 4: Convert to DataFrame
df = diamond_hits_to_dataframe(hits)
# pandas DataFrame with 21 columns

# Step 5: Analyze!
high_confidence = df[(df['evalue'] < 1e-50) & (df['pident'] > 90)]
by_subject = df.groupby('subject_id').agg({
    'protein_id': 'count',
    'pident': 'mean'
})
df.to_csv('results.csv', index=False)

┌────────────────────────────────────────────────────────────────────────┐
│                    WHAT YOU GET IN THE DATAFRAME                        │
└────────────────────────────────────────────────────────────────────────┘

  file_basename  contig_id   protein_id   start  end  strand  subject_id        pident  evalue     genomic_sequence
0 my_genome      contig_1    contig_1_1   100    300  1       IS110_transposase 95.5    1.00e-100  ATCGATCG...
1 my_genome      contig_1    contig_1_2   500    800  -1      IS110_orfB        88.2    1.00e-75   GCTAGCTA...
2 my_genome      contig_2    contig_2_1   50     250  1       IS110_transposase 92.0    1.00e-90   TACGTACG...

All information in one place:
✓ Source file (file_basename)
✓ Location (contig_id, start, end, strand)  
✓ Identity (protein_id)
✓ Match info (subject_id, pident, evalue, bitscore)
✓ Sequence (genomic_sequence)

Ready for pandas operations:
✓ Filter: df[df['evalue'] < 1e-50]
✓ Sort: df.sort_values('pident', ascending=False)
✓ Group: df.groupby('subject_id')['pident'].mean()
✓ Stats: df.describe()
✓ Export: df.to_csv(), df.to_excel()

┌────────────────────────────────────────────────────────────────────────┐
│                         SAME FOR HMMER!                                 │
└────────────────────────────────────────────────────────────────────────┘

from utils.parsers import (
    parse_hmm_tblout_with_positions,
    hmm_hits_to_dataframe
)

# Parse HMMER
hits = parse_hmm_tblout_with_positions(
    "hmmer_results.tbl",
    position_map,
    file_basename="my_genome"
)

# Convert to DataFrame
df = hmm_hits_to_dataframe(hits)

# Analyze
high_score = df[df['score'] > 100]
by_profile = df.groupby('query_name').size()
```
