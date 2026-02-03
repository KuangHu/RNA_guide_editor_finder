"""
Demo script showing pandas DataFrame analysis with parsers

This demonstrates:
1. Converting Diamond BLASTP results to DataFrame
2. Converting HMMER results to DataFrame
3. Filtering and analyzing data with pandas
4. Exporting results to CSV/Excel
5. Common data analysis workflows

Usage:
    conda activate opfi
    python examples/demo_dataframe_analysis.py
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import pandas as pd
    print("✓ pandas is available")
except ImportError:
    print("✗ pandas is not installed")
    print("Install with: pip install pandas")
    sys.exit(1)

from utils.parsers import (
    parse_prodigal_faa,
    parse_diamond_blastp_with_positions,
    parse_hmm_tblout_with_positions,
    retrieve_sequences_from_hits,
    diamond_hits_to_dataframe,
    hmm_hits_to_dataframe,
    save_hits_to_csv,
    DiamondBlastHit,
    HmmHit,
)


def demo_diamond_dataframe():
    """Demo: Convert Diamond BLASTP results to DataFrame"""
    print("=" * 70)
    print("DEMO 1: Diamond BLASTP to DataFrame")
    print("=" * 70)

    # Create sample Diamond hits
    sample_hits = [
        DiamondBlastHit(
            file_basename="genome1",
            contig_id="contig_1",
            protein_id="contig_1_1",
            start=100,
            end=300,
            strand=1,
            query_id="contig_1_1",
            subject_id="IS110_transposase",
            pident=95.5,
            alignment_length=200,
            mismatches=9,
            gap_opens=0,
            query_start=1,
            query_end=200,
            subject_start=15,
            subject_end=214,
            evalue=1e-100,
            bitscore=350.5
        ),
        DiamondBlastHit(
            file_basename="genome1",
            contig_id="contig_1",
            protein_id="contig_1_2",
            start=500,
            end=800,
            strand=-1,
            query_id="contig_1_2",
            subject_id="IS110_orfB",
            pident=88.2,
            alignment_length=150,
            mismatches=18,
            gap_opens=1,
            query_start=1,
            query_end=150,
            subject_start=10,
            subject_end=159,
            evalue=1e-75,
            bitscore=280.3
        ),
        DiamondBlastHit(
            file_basename="genome2",
            contig_id="contig_2",
            protein_id="contig_2_1",
            start=50,
            end=250,
            strand=1,
            query_id="contig_2_1",
            subject_id="IS110_transposase",
            pident=92.0,
            alignment_length=180,
            mismatches=14,
            gap_opens=0,
            query_start=1,
            query_end=180,
            subject_start=5,
            subject_end=184,
            evalue=1e-90,
            bitscore=320.1
        ),
    ]

    print("\n1. Convert to DataFrame:")
    print("-" * 70)
    df = diamond_hits_to_dataframe(sample_hits)
    print(df.head())

    print("\n2. DataFrame info:")
    print("-" * 70)
    print(df.info())

    print("\n3. Basic statistics:")
    print("-" * 70)
    print(df[['pident', 'evalue', 'bitscore', 'length']].describe())

    print("\n4. Filter by E-value:")
    print("-" * 70)
    high_confidence = df[df['evalue'] < 1e-80]
    print(f"Hits with E-value < 1e-80: {len(high_confidence)}")
    print(high_confidence[['protein_id', 'subject_id', 'evalue', 'pident']])

    print("\n5. Filter by identity:")
    print("-" * 70)
    high_identity = df[df['pident'] > 90]
    print(f"Hits with identity > 90%: {len(high_identity)}")
    print(high_identity[['protein_id', 'subject_id', 'pident']])

    print("\n6. Group by subject:")
    print("-" * 70)
    by_subject = df.groupby('subject_id').agg({
        'protein_id': 'count',
        'pident': 'mean',
        'evalue': 'min',
        'bitscore': 'max'
    })
    by_subject.columns = ['count', 'mean_pident', 'min_evalue', 'max_bitscore']
    print(by_subject)

    print("\n7. Group by genome:")
    print("-" * 70)
    by_genome = df.groupby('file_basename').agg({
        'protein_id': 'count',
        'pident': 'mean'
    })
    by_genome.columns = ['num_hits', 'mean_identity']
    print(by_genome)

    print("\n8. Filter by strand:")
    print("-" * 70)
    forward = df[df['is_forward']]
    reverse = df[~df['is_forward']]
    print(f"Forward strand: {len(forward)}")
    print(f"Reverse strand: {len(reverse)}")

    print("\n✓ Diamond DataFrame demo complete!\n")


def demo_hmmer_dataframe():
    """Demo: Convert HMMER results to DataFrame"""
    print("=" * 70)
    print("DEMO 2: HMMER to DataFrame")
    print("=" * 70)

    # Create sample HMMER hits
    sample_hits = [
        HmmHit(
            file_basename="genome1",
            contig_id="contig_1",
            protein_id="contig_1_1",
            start=100,
            end=300,
            strand=1,
            query_name="IS110_transposase",
            target_name="contig_1_1",
            evalue=1.5e-100,
            score=350.5,
            bias=0.1
        ),
        HmmHit(
            file_basename="genome1",
            contig_id="contig_1",
            protein_id="contig_1_2",
            start=500,
            end=800,
            strand=-1,
            query_name="IS110_orfB",
            target_name="contig_1_2",
            evalue=2.3e-75,
            score=280.3,
            bias=0.2
        ),
    ]

    print("\n1. Convert to DataFrame:")
    print("-" * 70)
    df = hmm_hits_to_dataframe(sample_hits)
    print(df.head())

    print("\n2. Filter by score:")
    print("-" * 70)
    high_score = df[df['score'] > 300]
    print(f"Hits with score > 300: {len(high_score)}")
    print(high_score[['protein_id', 'query_name', 'score', 'evalue']])

    print("\n3. Group by HMM profile:")
    print("-" * 70)
    by_profile = df.groupby('query_name').agg({
        'protein_id': 'count',
        'score': 'mean',
        'evalue': 'min'
    })
    by_profile.columns = ['num_hits', 'mean_score', 'min_evalue']
    print(by_profile)

    print("\n✓ HMMER DataFrame demo complete!\n")


def demo_filtering_and_export():
    """Demo: Advanced filtering and export"""
    print("=" * 70)
    print("DEMO 3: Advanced Filtering and Export")
    print("=" * 70)

    # Create sample hits
    sample_hits = [
        DiamondBlastHit(
            file_basename="genome1",
            contig_id="contig_1",
            protein_id=f"contig_1_{i}",
            start=100 * i,
            end=100 * i + 200,
            strand=1 if i % 2 == 0 else -1,
            query_id=f"contig_1_{i}",
            subject_id="IS110_transposase" if i % 3 == 0 else "IS110_orfB",
            pident=90.0 + i,
            alignment_length=200,
            mismatches=10 - i,
            gap_opens=0,
            query_start=1,
            query_end=200,
            subject_start=1,
            subject_end=200,
            evalue=10 ** (-50 - i * 10),
            bitscore=300.0 + i * 10
        )
        for i in range(10)
    ]

    df = diamond_hits_to_dataframe(sample_hits)

    print("\n1. Complex filtering:")
    print("-" * 70)
    filtered = df[
        (df['evalue'] < 1e-70) &
        (df['pident'] > 92) &
        (df['is_forward'])
    ]
    print(f"Hits matching all criteria: {len(filtered)}")
    print(filtered[['protein_id', 'subject_id', 'pident', 'evalue']])

    print("\n2. Sort by multiple columns:")
    print("-" * 70)
    sorted_df = df.sort_values(['subject_id', 'evalue'])
    print(sorted_df[['protein_id', 'subject_id', 'evalue', 'bitscore']].head())

    print("\n3. Add custom columns:")
    print("-" * 70)
    df['log_evalue'] = df['evalue'].apply(lambda x: -1 * pd.np.log10(x) if x > 0 else 300)
    df['identity_category'] = pd.cut(df['pident'], bins=[0, 80, 90, 95, 100],
                                      labels=['low', 'medium', 'high', 'very_high'])
    print(df[['protein_id', 'pident', 'identity_category', 'evalue', 'log_evalue']].head())

    print("\n4. Summary statistics by category:")
    print("-" * 70)
    summary = df.groupby('identity_category').agg({
        'protein_id': 'count',
        'evalue': 'mean',
        'bitscore': 'mean'
    })
    print(summary)

    print("\n5. Export to CSV:")
    print("-" * 70)
    # Note: This would actually save files in a real scenario
    print("  save_hits_to_csv(sample_hits, 'diamond_results.csv')")
    print("  ✓ Would save to diamond_results.csv")

    print("\n✓ Filtering and export demo complete!\n")


def demo_comparison_analysis():
    """Demo: Compare results across genomes"""
    print("=" * 70)
    print("DEMO 4: Cross-Genome Comparison")
    print("=" * 70)

    # Create sample hits from multiple genomes
    hits_genome1 = [
        DiamondBlastHit(
            file_basename="genome1",
            contig_id=f"contig_{i}",
            protein_id=f"genome1_protein_{i}",
            start=100,
            end=300,
            strand=1,
            query_id=f"genome1_protein_{i}",
            subject_id="IS110_transposase",
            pident=90.0 + i,
            alignment_length=200,
            mismatches=10,
            gap_opens=0,
            query_start=1,
            query_end=200,
            subject_start=1,
            subject_end=200,
            evalue=1e-90,
            bitscore=300.0
        )
        for i in range(5)
    ]

    hits_genome2 = [
        DiamondBlastHit(
            file_basename="genome2",
            contig_id=f"contig_{i}",
            protein_id=f"genome2_protein_{i}",
            start=100,
            end=300,
            strand=1,
            query_id=f"genome2_protein_{i}",
            subject_id="IS110_transposase",
            pident=85.0 + i,
            alignment_length=200,
            mismatches=15,
            gap_opens=0,
            query_start=1,
            query_end=200,
            subject_start=1,
            subject_end=200,
            evalue=1e-80,
            bitscore=280.0
        )
        for i in range(3)
    ]

    # Combine all hits
    all_hits = hits_genome1 + hits_genome2
    df = diamond_hits_to_dataframe(all_hits)

    print("\n1. Hits per genome:")
    print("-" * 70)
    print(df['file_basename'].value_counts())

    print("\n2. Compare identity across genomes:")
    print("-" * 70)
    comparison = df.groupby('file_basename').agg({
        'pident': ['mean', 'min', 'max', 'std'],
        'protein_id': 'count'
    })
    print(comparison)

    print("\n3. Find best hit per genome:")
    print("-" * 70)
    best_hits = df.loc[df.groupby('file_basename')['bitscore'].idxmax()]
    print(best_hits[['file_basename', 'protein_id', 'pident', 'bitscore']])

    print("\n✓ Comparison analysis demo complete!\n")


def main():
    """Run all demos"""
    print("\n" + "=" * 70)
    print("PANDAS DATAFRAME ANALYSIS DEMO")
    print("=" * 70)
    print("\nThis demo shows how to convert parser results to DataFrames")
    print("and perform common analysis tasks with pandas.")
    print()

    # Run demos
    demo_diamond_dataframe()
    demo_hmmer_dataframe()
    demo_filtering_and_export()
    demo_comparison_analysis()

    print("=" * 70)
    print("All demos complete!")
    print("=" * 70)
    print("\nKey takeaways:")
    print("1. Use diamond_hits_to_dataframe() or hmm_hits_to_dataframe()")
    print("2. DataFrames support powerful filtering and grouping")
    print("3. Easy export to CSV/Excel with save_hits_to_csv()")
    print("4. Leverage pandas for statistical analysis")
    print()


if __name__ == "__main__":
    main()
