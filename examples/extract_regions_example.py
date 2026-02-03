#!/usr/bin/env python3
"""
Example: Using RegionExtractor to extract transposon regions

This demonstrates the single-protein version of region extraction.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from modules.region_extractor import (
    RegionExtractor,
    create_transposon_dict,
    parse_prodigal_faa,
    extract_contig_from_protein_id
)


def example_1_basic_usage():
    """Example 1: Basic usage with RegionExtractor class"""
    print("=" * 80)
    print("Example 1: Basic Usage")
    print("=" * 80)

    # Initialize extractor with genome file
    genome_path = "/path/to/genome.fna"
    extractor = RegionExtractor(genome_path)

    # Extract regions for a specific protein
    result = extractor.extract_protein_regions(
        contig_id="AAXX01000001.1",
        start=1000,
        end=2000,
        strand=1,  # Forward strand
        upstream_length=350,
        downstream_length=250
    )

    print(f"Contig: {result['contig_id']}")
    print(f"Strand: {result['strand']}")
    print(f"Coding length: {result['coding_coords']['length']} bp")
    print(f"Upstream length: {result['upstream_coords']['length']} bp")
    print(f"Downstream length: {result['downstream_coords']['length']} bp")
    print(f"Coding sequence: {result['coding_sequence'][:50]}...")
    print()


def example_2_convenience_function():
    """Example 2: Using the convenience function"""
    print("=" * 80)
    print("Example 2: Convenience Function")
    print("=" * 80)

    transposon = create_transposon_dict(
        genome_fna_path="/path/to/genome.fna",
        contig_id="AAXX01000001.1",
        protein_id="AAXX01000001.1_123",
        start=1000,
        end=2000,
        strand=-1,  # Reverse strand
        upstream_length=350,
        downstream_length=250,
        additional_metadata={
            'domain_info': {
                'DEDD': {'evalue': 1.2e-50, 'start': 10, 'end': 150},
                'Tnp20': {'evalue': 3.4e-40, 'start': 200, 'end': 350}
            }
        }
    )

    print(f"Protein ID: {transposon['protein_id']}")
    print(f"Coordinates:")
    print(f"  Coding: {transposon['coordinates']['coding']['start']}-"
          f"{transposon['coordinates']['coding']['end']}")
    print(f"  Upstream: {transposon['coordinates']['upstream']['start']}-"
          f"{transposon['coordinates']['upstream']['end']}")
    print(f"  Downstream: {transposon['coordinates']['downstream']['start']}-"
          f"{transposon['coordinates']['downstream']['end']}")
    print()


def example_3_with_prodigal():
    """Example 3: Integration with Prodigal FAA parsing"""
    print("=" * 80)
    print("Example 3: With Prodigal FAA File")
    print("=" * 80)

    # Parse Prodigal FAA file to get coordinates
    faa_path = "/path/to/proteins.faa"
    positions = parse_prodigal_faa(faa_path)

    # Initialize extractor
    genome_path = "/path/to/genome.fna"
    extractor = RegionExtractor(genome_path)

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

        print(f"Extracted regions for {protein_id}")
        print(f"  Coding: {result['coding_coords']['length']} bp")
        print(f"  Upstream: {result['upstream_coords']['length']} bp")
        print(f"  Downstream: {result['downstream_coords']['length']} bp")
    print()


def example_4_batch_processing():
    """Example 4: Batch processing multiple proteins"""
    print("=" * 80)
    print("Example 4: Batch Processing")
    print("=" * 80)

    genome_path = "/path/to/genome.fna"
    extractor = RegionExtractor(genome_path)

    # List of proteins to process
    proteins = [
        {"protein_id": "AAXX01000001.1_123", "start": 1000, "end": 2000, "strand": 1},
        {"protein_id": "AAXX01000001.1_124", "start": 3000, "end": 4000, "strand": -1},
        {"protein_id": "AAXX01000002.1_45", "start": 500, "end": 1500, "strand": 1},
    ]

    results = []
    for protein in proteins:
        contig_id = extract_contig_from_protein_id(protein['protein_id'])

        result = extractor.extract_protein_regions(
            contig_id=contig_id,
            start=protein['start'],
            end=protein['end'],
            strand=protein['strand']
        )

        transposon = {
            'protein_id': protein['protein_id'],
            'contig_id': contig_id,
            'strand': protein['strand'],
            'coordinates': {
                'coding': result['coding_coords'],
                'upstream': result['upstream_coords'],
                'downstream': result['downstream_coords']
            },
            'sequences': {
                'coding': result['coding_sequence'],
                'upstream': result['upstream_sequence'],
                'downstream': result['downstream_sequence']
            }
        }

        results.append(transposon)

    print(f"Processed {len(results)} proteins")
    print()


def example_5_real_workflow():
    """
    Example 5: Real workflow similar to the GTDB1 script

    This shows how to integrate everything together.
    """
    print("=" * 80)
    print("Example 5: Complete Workflow")
    print("=" * 80)
    print("This example shows integration with CSV input and JSON output")
    print()

    # Pseudo-code for a complete workflow:
    print("""
    import csv
    import json
    from collections import defaultdict
    from modules.region_extractor import RegionExtractor, parse_prodigal_faa

    # 1. Read input CSV with protein information
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
        # Load genome and Prodigal annotations
        fna_path = f"/path/to/genomes/{genome_name}.fna"
        faa_path = f"/path/to/proteins/{genome_name}.faa"

        extractor = RegionExtractor(fna_path)
        positions = parse_prodigal_faa(faa_path)

        # Process each protein
        for row in protein_rows:
            protein_id = row['protein_id']
            if protein_id in positions:
                start, end, strand = positions[protein_id]
                contig_id = extract_contig_from_protein_id(protein_id)

                result = extractor.extract_protein_regions(
                    contig_id=contig_id,
                    start=start,
                    end=end,
                    strand=strand,
                    upstream_length=350,
                    downstream_length=250
                )

                # Build complete entry with domain info, etc.
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
                    'domains': {
                        'DEDD': {
                            'evalue': row['DEDD_evalue'],
                            'start': row['DEDD_start'],
                            'end': row['DEDD_end']
                        },
                        'Tnp20': {
                            'evalue': row['Tnp20_evalue'],
                            'start': row['Tnp20_start'],
                            'end': row['Tnp20_end']
                        }
                    }
                }

                all_results.append(transposon)

    # 4. Save results
    with open('output.json', 'w') as f:
        json.dump(all_results, f, indent=2)

    print(f"Saved {len(all_results)} transposons to output.json")
    """)


if __name__ == "__main__":
    print("RegionExtractor Usage Examples")
    print()

    print("Note: These examples use placeholder paths.")
    print("Replace with your actual file paths to run.")
    print()

    # Uncomment to run examples with real data
    # example_1_basic_usage()
    # example_2_convenience_function()
    # example_3_with_prodigal()
    # example_4_batch_processing()
    example_5_real_workflow()
