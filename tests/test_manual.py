#!/usr/bin/env python3
"""
Simple manual test for region_extractor module
"""

import tempfile
from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from modules.region_extractor import (
    RegionExtractor,
    reverse_complement,
    extract_contig_from_protein_id,
    create_transposon_dict
)


def test_reverse_complement():
    """Test reverse complement function"""
    print("Testing reverse_complement...")
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("ATGCN") == "NGCAT"
    assert reverse_complement("atgc") == "gcat"
    print("  ✓ reverse_complement works correctly")


def test_extract_contig():
    """Test contig extraction"""
    print("Testing extract_contig_from_protein_id...")
    assert extract_contig_from_protein_id("AAXX01000001.1_281") == "AAXX01000001.1"
    assert extract_contig_from_protein_id("NC_000001.11_123") == "NC_000001.11"
    print("  ✓ extract_contig_from_protein_id works correctly")


def test_region_extractor():
    """Test RegionExtractor with a temporary genome"""
    print("Testing RegionExtractor...")

    # Create a test genome
    sequence = "ATGCATGCATGCATGC" * 100  # 1600 bp
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
        f.write(">test_contig\n")
        f.write(sequence + "\n")
        temp_path = f.name

    try:
        # Initialize extractor
        extractor = RegionExtractor(temp_path)
        print(f"  Loaded genome: {len(extractor.genome_sequences['test_contig'])} bp")

        # Test forward strand extraction
        result = extractor.extract_protein_regions(
            contig_id="test_contig",
            start=500,
            end=700,
            strand=1,
            upstream_length=100,
            downstream_length=100
        )

        print(f"  Forward strand:")
        print(f"    Coding: {result['coding_coords']['start']}-{result['coding_coords']['end']} ({result['coding_coords']['length']} bp)")
        print(f"    Upstream: {result['upstream_coords']['start']}-{result['upstream_coords']['end']} ({result['upstream_coords']['length']} bp)")
        print(f"    Downstream: {result['downstream_coords']['start']}-{result['downstream_coords']['end']} ({result['downstream_coords']['length']} bp)")

        assert result['coding_coords']['length'] == 201  # 1-based inclusive: 700-500+1
        assert result['upstream_coords']['length'] == 100
        assert result['downstream_coords']['length'] == 100

        # Test reverse strand extraction
        result_rev = extractor.extract_protein_regions(
            contig_id="test_contig",
            start=500,
            end=700,
            strand=-1,
            upstream_length=100,
            downstream_length=100
        )

        print(f"  Reverse strand:")
        print(f"    Coding: {result_rev['coding_coords']['start']}-{result_rev['coding_coords']['end']} ({result_rev['coding_coords']['length']} bp)")
        print(f"    Upstream: {result_rev['upstream_coords']['start']}-{result_rev['upstream_coords']['end']} ({result_rev['upstream_coords']['length']} bp)")
        print(f"    Downstream: {result_rev['downstream_coords']['start']}-{result_rev['downstream_coords']['end']} ({result_rev['downstream_coords']['length']} bp)")

        # Verify reverse complementation
        assert result_rev['coding_sequence'] != result['coding_sequence']
        assert len(result_rev['coding_sequence']) == len(result['coding_sequence'])

        print("  ✓ RegionExtractor works correctly")

    finally:
        Path(temp_path).unlink()


def test_create_transposon_dict_func():
    """Test create_transposon_dict convenience function"""
    print("Testing create_transposon_dict...")

    # Create a test genome
    sequence = "ATGCATGCATGCATGC" * 100  # 1600 bp
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
        f.write(">test_contig\n")
        f.write(sequence + "\n")
        temp_path = f.name

    try:
        result = create_transposon_dict(
            genome_fna_path=temp_path,
            contig_id="test_contig",
            protein_id="test_protein_123",
            start=500,
            end=700,
            strand=1,
            upstream_length=100,
            downstream_length=100,
            additional_metadata={
                'annotation': 'IS110 transposase',
                'domains': {'DEDD': {'evalue': 1e-50}}
            }
        )

        print(f"  Created transposon dict:")
        print(f"    Protein ID: {result['protein_id']}")
        print(f"    Contig ID: {result['contig_id']}")
        print(f"    Strand: {result['strand']}")
        print(f"    Coding length: {result['coordinates']['coding']['length']} bp")
        print(f"    Has metadata: {('annotation' in result) and ('domains' in result)}")

        assert result['protein_id'] == "test_protein_123"
        assert result['contig_id'] == "test_contig"
        assert result['annotation'] == 'IS110 transposase'
        assert 'DEDD' in result['domains']

        print("  ✓ create_transposon_dict works correctly")

    finally:
        Path(temp_path).unlink()


def main():
    print("=" * 80)
    print("Manual Tests for region_extractor Module")
    print("=" * 80)
    print()

    try:
        test_reverse_complement()
        print()
        test_extract_contig()
        print()
        test_region_extractor()
        print()
        test_create_transposon_dict_func()
        print()
        print("=" * 80)
        print("All tests passed! ✓")
        print("=" * 80)
    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        return 1
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
