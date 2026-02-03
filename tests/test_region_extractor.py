#!/usr/bin/env python3
"""
Unit tests for region_extractor module
"""

import pytest
import tempfile
from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from modules.region_extractor import (
    RegionExtractor,
    reverse_complement,
    parse_prodigal_faa,
    extract_contig_from_protein_id,
    create_transposon_dict
)


class TestReverseComplement:
    """Tests for reverse_complement function"""

    def test_simple_sequence(self):
        assert reverse_complement("ATGC") == "GCAT"

    def test_with_n(self):
        assert reverse_complement("ATGCN") == "NGCAT"

    def test_lowercase(self):
        assert reverse_complement("atgc") == "gcat"

    def test_mixed_case(self):
        assert reverse_complement("ATgc") == "gcAT"

    def test_ambiguous_bases(self):
        assert reverse_complement("RYWSMK") == "MKSWRY"


class TestExtractContigFromProteinId:
    """Tests for extract_contig_from_protein_id function"""

    def test_standard_format(self):
        assert extract_contig_from_protein_id("AAXX01000001.1_281") == "AAXX01000001.1"

    def test_ncbi_format(self):
        assert extract_contig_from_protein_id("NC_000001.11_123") == "NC_000001.11"

    def test_multiple_underscores(self):
        assert extract_contig_from_protein_id("contig_name_test_456") == "contig_name_test"


class TestParseProdigalFaa:
    """Tests for parse_prodigal_faa function"""

    def test_parse_prodigal_header(self):
        # Create a temporary FAA file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.faa') as f:
            f.write(">protein_1 # 100 # 500 # 1 # ID=1_1;partial=00\n")
            f.write("MTKQVLAAAA\n")
            f.write(">protein_2 # 1000 # 2000 # -1 # ID=1_2;partial=00\n")
            f.write("MKKQVLAAAA\n")
            temp_path = f.name

        try:
            positions = parse_prodigal_faa(temp_path)

            assert "protein_1" in positions
            assert positions["protein_1"] == (100, 500, 1)

            assert "protein_2" in positions
            assert positions["protein_2"] == (1000, 2000, -1)
        finally:
            Path(temp_path).unlink()


class TestRegionExtractor:
    """Tests for RegionExtractor class"""

    @pytest.fixture
    def test_genome_file(self):
        """Create a test genome FASTA file"""
        # Create a simple test genome
        # Contig with 1000 bp
        sequence = "ATGC" * 250  # 1000 bp total

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            f.write(">test_contig_1\n")
            f.write(sequence + "\n")
            temp_path = f.name

        yield temp_path

        # Cleanup
        Path(temp_path).unlink()

    def test_load_genome(self, test_genome_file):
        """Test genome loading"""
        extractor = RegionExtractor(test_genome_file)
        assert "test_contig_1" in extractor.genome_sequences
        assert len(extractor.genome_sequences["test_contig_1"]) == 1000

    def test_extract_forward_strand(self, test_genome_file):
        """Test extraction on forward strand"""
        extractor = RegionExtractor(test_genome_file)

        result = extractor.extract_protein_regions(
            contig_id="test_contig_1",
            start=100,
            end=200,
            strand=1,
            upstream_length=50,
            downstream_length=50
        )

        assert result['strand'] == 1
        assert result['coding_coords']['start'] == 100
        assert result['coding_coords']['end'] == 200
        # 1-based inclusive: positions 100-200 = 101 bases
        assert result['coding_coords']['length'] == 101

        # Check upstream (positions 50-99 = 50 bases)
        assert result['upstream_coords']['start'] == 50
        assert result['upstream_coords']['end'] == 99
        assert result['upstream_coords']['length'] == 50

        # Check downstream (positions 201-250 = 50 bases)
        assert result['downstream_coords']['start'] == 201
        assert result['downstream_coords']['end'] == 250
        assert result['downstream_coords']['length'] == 50

    def test_extract_reverse_strand(self, test_genome_file):
        """Test extraction on reverse strand"""
        extractor = RegionExtractor(test_genome_file)

        result = extractor.extract_protein_regions(
            contig_id="test_contig_1",
            start=100,
            end=200,
            strand=-1,
            upstream_length=50,
            downstream_length=50
        )

        assert result['strand'] == -1
        assert result['coding_coords']['start'] == 100
        assert result['coding_coords']['end'] == 200

        # For reverse strand:
        # - upstream is genomically downstream (after end)
        # - downstream is genomically upstream (before start)
        assert result['upstream_coords']['start'] == 201
        assert result['upstream_coords']['end'] == 250

        assert result['downstream_coords']['start'] == 50
        assert result['downstream_coords']['end'] == 99

        # Sequences should be reverse complemented
        assert result['coding_sequence'] != extractor.genome_sequences["test_contig_1"][99:200]

    def test_boundary_conditions_start(self, test_genome_file):
        """Test extraction at start of contig"""
        extractor = RegionExtractor(test_genome_file)

        result = extractor.extract_protein_regions(
            contig_id="test_contig_1",
            start=1,
            end=100,
            strand=1,
            upstream_length=500,  # More than available
            downstream_length=50
        )

        # Upstream should be limited to what's available
        assert result['upstream_coords']['length'] == 0
        assert result['coding_coords']['length'] == 100

    def test_boundary_conditions_end(self, test_genome_file):
        """Test extraction at end of contig"""
        extractor = RegionExtractor(test_genome_file)

        result = extractor.extract_protein_regions(
            contig_id="test_contig_1",
            start=900,
            end=1000,
            strand=1,
            upstream_length=50,
            downstream_length=500  # More than available
        )

        # Downstream should be limited to what's available
        assert result['downstream_coords']['length'] == 0
        # 1-based inclusive: positions 900-1000 = 101 bases
        assert result['coding_coords']['length'] == 101

    def test_invalid_contig(self, test_genome_file):
        """Test error handling for invalid contig"""
        extractor = RegionExtractor(test_genome_file)

        with pytest.raises(KeyError):
            extractor.extract_protein_regions(
                contig_id="nonexistent_contig",
                start=100,
                end=200,
                strand=1
            )

    def test_invalid_coordinates(self, test_genome_file):
        """Test error handling for invalid coordinates"""
        extractor = RegionExtractor(test_genome_file)

        with pytest.raises(ValueError):
            extractor.extract_protein_regions(
                contig_id="test_contig_1",
                start=200,
                end=100,  # End before start
                strand=1
            )

    def test_invalid_strand(self, test_genome_file):
        """Test error handling for invalid strand"""
        extractor = RegionExtractor(test_genome_file)

        with pytest.raises(ValueError):
            extractor.extract_protein_regions(
                contig_id="test_contig_1",
                start=100,
                end=200,
                strand=0  # Invalid strand
            )


class TestCreateTransposonDict:
    """Tests for create_transposon_dict convenience function"""

    @pytest.fixture
    def test_genome_file(self):
        """Create a test genome FASTA file"""
        sequence = "ATGC" * 250  # 1000 bp total

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fna') as f:
            f.write(">test_contig_1\n")
            f.write(sequence + "\n")
            temp_path = f.name

        yield temp_path
        Path(temp_path).unlink()

    def test_create_transposon_dict_basic(self, test_genome_file):
        """Test basic transposon dict creation"""
        result = create_transposon_dict(
            genome_fna_path=test_genome_file,
            contig_id="test_contig_1",
            protein_id="test_protein_1",
            start=100,
            end=200,
            strand=1
        )

        assert result['protein_id'] == "test_protein_1"
        assert result['contig_id'] == "test_contig_1"
        assert result['strand'] == 1
        assert 'coordinates' in result
        assert 'sequences' in result

    def test_create_transposon_dict_with_metadata(self, test_genome_file):
        """Test transposon dict with additional metadata"""
        metadata = {
            'domain_info': {
                'DEDD': {'evalue': 1e-50}
            },
            'annotation': 'IS110 transposase'
        }

        result = create_transposon_dict(
            genome_fna_path=test_genome_file,
            contig_id="test_contig_1",
            protein_id="test_protein_1",
            start=100,
            end=200,
            strand=1,
            additional_metadata=metadata
        )

        assert result['domain_info'] == metadata['domain_info']
        assert result['annotation'] == metadata['annotation']


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
