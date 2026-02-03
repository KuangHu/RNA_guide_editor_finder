"""
Demo script showing how to use the parsers module

This demonstrates typical workflows for:
1. Diamond BLASTP results with sequence retrieval
2. HMMER results with sequence retrieval
3. Integration with Prodigal output

Usage:
    python examples/demo_parsers.py
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.parsers import (
    # Prodigal parsers
    parse_prodigal_faa,
    parse_prodigal_faa_full,
    # Diamond BLASTP parsers
    parse_diamond_blastp_with_positions,
    # HMMER parsers
    parse_hmm_tblout_with_positions,
    # Sequence retrieval
    retrieve_sequences_from_hits,
    retrieve_sequence_from_fasta,
)


def demo_diamond_workflow():
    """
    Demo: Complete Diamond BLASTP workflow

    Steps:
    1. Parse Prodigal FAA to get protein positions
    2. Parse Diamond BLASTP results with positions
    3. Retrieve sequences from FASTA files
    """
    print("=" * 70)
    print("DEMO 1: Diamond BLASTP Workflow")
    print("=" * 70)

    # File paths (adjust these to your data)
    prodigal_faa = "data/genome.faa"
    diamond_output = "data/diamond_results.m8"
    genome_fasta = "data/genome.fna"

    print("\nStep 1: Parse Prodigal FAA file to get protein positions...")
    print(f"  File: {prodigal_faa}")

    try:
        position_map = parse_prodigal_faa(prodigal_faa)
        print(f"  Found positions for {len(position_map)} proteins")

        # Show example
        if position_map:
            protein_id = list(position_map.keys())[0]
            start, end, strand = position_map[protein_id]
            print(f"  Example: {protein_id} -> {start}..{end} (strand={strand})")
    except FileNotFoundError:
        print(f"  ⚠️  File not found (this is just a demo)")
        return

    print("\nStep 2: Parse Diamond BLASTP results with genomic positions...")
    print(f"  File: {diamond_output}")

    hits = parse_diamond_blastp_with_positions(
        diamond_output,
        position_map,
        file_basename="my_genome"
    )

    print(f"  Found {len(hits)} hits")

    if hits:
        hit = hits[0]
        print(f"  Example hit:")
        print(f"    File: {hit.file_basename}")
        print(f"    Contig: {hit.contig_id}")
        print(f"    Protein: {hit.protein_id}")
        print(f"    Position: {hit.start}..{hit.end} (strand={hit.strand})")
        print(f"    Subject: {hit.subject_id}")
        print(f"    Identity: {hit.pident}%")
        print(f"    E-value: {hit.evalue}")

    print("\nStep 3: Retrieve sequences from FASTA files...")
    print(f"  FASTA file: {genome_fasta}")

    hits_with_sequences = retrieve_sequences_from_hits(
        hits,
        [genome_fasta],
        verbose=True
    )

    # Show sequences
    print("\n  Sequences retrieved:")
    for hit in hits_with_sequences[:3]:  # Show first 3
        if hasattr(hit, 'genomic_sequence') and hit.genomic_sequence:
            seq = hit.genomic_sequence
            print(f"    {hit.protein_id}: {seq[:50]}... ({len(seq)} bp)")

    print("\n✓ Diamond workflow complete!\n")


def demo_hmmer_workflow():
    """
    Demo: Complete HMMER workflow

    Steps:
    1. Parse Prodigal FAA to get protein positions
    2. Parse HMMER tblout results with positions
    3. Retrieve sequences from FASTA files
    """
    print("=" * 70)
    print("DEMO 2: HMMER Workflow")
    print("=" * 70)

    # File paths (adjust these to your data)
    prodigal_faa = "data/genome.faa"
    hmmer_tblout = "data/hmmer_results.tbl"
    genome_fasta = "data/genome.fna"

    print("\nStep 1: Parse Prodigal FAA file to get protein positions...")
    print(f"  File: {prodigal_faa}")

    try:
        position_map = parse_prodigal_faa(prodigal_faa)
        print(f"  Found positions for {len(position_map)} proteins")
    except FileNotFoundError:
        print(f"  ⚠️  File not found (this is just a demo)")
        return

    print("\nStep 2: Parse HMMER tblout results with genomic positions...")
    print(f"  File: {hmmer_tblout}")

    hits = parse_hmm_tblout_with_positions(
        hmmer_tblout,
        position_map,
        file_basename="my_genome"
    )

    print(f"  Found {len(hits)} hits")

    if hits:
        hit = hits[0]
        print(f"  Example hit:")
        print(f"    File: {hit.file_basename}")
        print(f"    Contig: {hit.contig_id}")
        print(f"    Protein: {hit.protein_id}")
        print(f"    Position: {hit.start}..{hit.end} (strand={hit.strand})")
        print(f"    Query HMM: {hit.query_name}")
        print(f"    E-value: {hit.evalue}")
        print(f"    Score: {hit.score}")

    print("\nStep 3: Retrieve sequences from FASTA files...")
    print(f"  FASTA file: {genome_fasta}")

    hits_with_sequences = retrieve_sequences_from_hits(
        hits,
        [genome_fasta],
        verbose=True
    )

    # Show sequences
    print("\n  Sequences retrieved:")
    for hit in hits_with_sequences[:3]:  # Show first 3
        if hasattr(hit, 'genomic_sequence') and hit.genomic_sequence:
            seq = hit.genomic_sequence
            print(f"    {hit.protein_id}: {seq[:50]}... ({len(seq)} bp)")

    print("\n✓ HMMER workflow complete!\n")


def demo_sequence_retrieval():
    """
    Demo: Direct sequence retrieval from FASTA
    """
    print("=" * 70)
    print("DEMO 3: Direct Sequence Retrieval")
    print("=" * 70)

    fasta_file = "data/genome.fna"

    print(f"\nRetrieving sequence from: {fasta_file}")
    print("  Contig: my_contig")
    print("  Coordinates: 1000..1100 (forward strand)")

    try:
        sequence = retrieve_sequence_from_fasta(
            fasta_file,
            contig_id="my_contig",
            start=1000,
            end=1100,
            strand=1
        )

        if sequence:
            print(f"\n  Retrieved sequence ({len(sequence)} bp):")
            print(f"    {sequence[:50]}...")
        else:
            print("\n  ⚠️  Sequence not found")
    except FileNotFoundError:
        print(f"  ⚠️  File not found (this is just a demo)")

    print("\n✓ Sequence retrieval complete!\n")


def demo_data_structures():
    """
    Demo: Show the data structures returned by parsers
    """
    print("=" * 70)
    print("DEMO 4: Data Structures")
    print("=" * 70)

    print("\nDiamondBlastHit fields:")
    print("  - file_basename: Source file identifier")
    print("  - contig_id: Contig/chromosome identifier")
    print("  - protein_id: Protein identifier")
    print("  - start: Genomic start position (1-based)")
    print("  - end: Genomic end position (1-based)")
    print("  - strand: 1 (forward) or -1 (reverse)")
    print("  - query_id: Query sequence ID")
    print("  - subject_id: Subject sequence ID (hit)")
    print("  - pident: Percentage identity")
    print("  - evalue: E-value")
    print("  - bitscore: Bit score")
    print("  - genomic_sequence: Retrieved sequence (added by retrieve_sequences_from_hits)")

    print("\nHmmHit fields:")
    print("  - file_basename: Source file identifier")
    print("  - contig_id: Contig/chromosome identifier")
    print("  - protein_id: Protein identifier")
    print("  - start: Genomic start position (1-based)")
    print("  - end: Genomic end position (1-based)")
    print("  - strand: 1 (forward) or -1 (reverse)")
    print("  - query_name: HMM profile name")
    print("  - target_name: Target sequence name")
    print("  - evalue: Full sequence E-value")
    print("  - score: Full sequence score")
    print("  - genomic_sequence: Retrieved sequence (added by retrieve_sequences_from_hits)")

    print("\nHelper properties:")
    print("  - hit.length: Genomic region length")
    print("  - hit.is_forward: True if on forward strand")

    print("\n✓ Data structure overview complete!\n")


def main():
    """Run all demos"""
    print("\n" + "=" * 70)
    print("PARSERS MODULE DEMO")
    print("=" * 70)
    print("\nThis demo shows how to use the parsers module for common workflows.")
    print("Note: File paths in this demo are examples and may not exist.")
    print()

    # Run demos
    demo_diamond_workflow()
    demo_hmmer_workflow()
    demo_sequence_retrieval()
    demo_data_structures()

    print("=" * 70)
    print("All demos complete!")
    print("=" * 70)
    print("\nTo use these parsers with your own data:")
    print("1. Replace file paths with your actual data files")
    print("2. Run Prodigal to generate .faa files")
    print("3. Run Diamond BLASTP or HMMER to generate results")
    print("4. Use the parsers to extract and analyze results")
    print()


if __name__ == "__main__":
    main()
