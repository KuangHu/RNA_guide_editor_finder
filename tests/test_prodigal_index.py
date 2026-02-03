#!/usr/bin/env python3
"""
Tests for the Prodigal Index module.
"""

import pytest
import tempfile
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.prodigal_index import (
    ProdigalIndex,
    ProteinPosition,
    build_prodigal_index,
    build_index_from_gff,
    _extract_contig_id,
    _parse_prodigal_header,
)


class TestHelperFunctions:
    """Test helper functions."""

    def test_extract_contig_id_standard(self):
        """Test contig extraction from standard protein ID."""
        assert _extract_contig_id("BA000021.3_1") == "BA000021.3"
        assert _extract_contig_id("NC_000001.11_123") == "NC_000001.11"
        assert _extract_contig_id("contig_1_456") == "contig_1"

    def test_extract_contig_id_no_underscore_number(self):
        """Test contig extraction when no gene number suffix."""
        assert _extract_contig_id("contig_name") == "contig_name"
        assert _extract_contig_id("simple") == "simple"

    def test_parse_prodigal_header_valid(self):
        """Test parsing valid Prodigal header."""
        header = ">BA000021.3_1 # 182 # 2068 # 1 # ID=1_1;partial=00"
        result = _parse_prodigal_header(header)

        assert result is not None
        protein_id, start, end, strand, contig_id = result
        assert protein_id == "BA000021.3_1"
        assert start == 182
        assert end == 2068
        assert strand == 1
        assert contig_id == "BA000021.3"

    def test_parse_prodigal_header_reverse_strand(self):
        """Test parsing reverse strand header."""
        header = ">contig_2 # 1000 # 2000 # -1 # ID=1_2"
        result = _parse_prodigal_header(header)

        assert result is not None
        _, start, end, strand, _ = result
        assert start == 1000
        assert end == 2000
        assert strand == -1

    def test_parse_prodigal_header_invalid(self):
        """Test parsing invalid headers."""
        assert _parse_prodigal_header("not a header") is None
        assert _parse_prodigal_header(">simple_header") is None
        assert _parse_prodigal_header(">id # not # numbers") is None


class TestBuildAndQuery:
    """Test index building and querying."""

    @pytest.fixture
    def sample_faa_file(self):
        """Create a sample FAA file for testing."""
        content = """>BA000021.3_1 # 182 # 2068 # 1 # ID=1_1;partial=00
MTKQVLAAAA
>BA000021.3_2 # 2230 # 3066 # 1 # ID=1_2;partial=00
MKKQVLBBBB
>BA000021.3_3 # 3076 # 3315 # -1 # ID=1_3;partial=00
MLLQVLCCCC
>NC_000001.11_1 # 100 # 500 # 1 # ID=2_1;partial=00
MDDQVLDDDD
>NC_000001.11_2 # 600 # 1200 # -1 # ID=2_2;partial=00
MEEQVLEEEE
"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.faa') as f:
            f.write(content)
            faa_path = f.name

        yield faa_path

        # Cleanup
        os.unlink(faa_path)
        db_path = faa_path + ".sqlite"
        if os.path.exists(db_path):
            os.unlink(db_path)

    def test_build_index(self, sample_faa_file):
        """Test building an index from FAA file."""
        db_path = sample_faa_file + ".sqlite"

        result = build_prodigal_index(
            sample_faa_file,
            db_path,
            show_progress=False
        )

        assert result == db_path
        assert os.path.exists(db_path)
        assert os.path.getsize(db_path) > 0

    def test_query_single(self, sample_faa_file):
        """Test single protein lookup."""
        db_path = sample_faa_file + ".sqlite"
        build_prodigal_index(sample_faa_file, db_path, show_progress=False)

        with ProdigalIndex(db_path) as idx:
            pos = idx.get_position("BA000021.3_1")

            assert pos is not None
            assert pos.protein_id == "BA000021.3_1"
            assert pos.start == 182
            assert pos.end == 2068
            assert pos.strand == 1
            assert pos.contig_id == "BA000021.3"
            assert pos.length == 2068 - 182 + 1

    def test_query_tuple(self, sample_faa_file):
        """Test tuple lookup."""
        db_path = sample_faa_file + ".sqlite"
        build_prodigal_index(sample_faa_file, db_path, show_progress=False)

        with ProdigalIndex(db_path) as idx:
            result = idx.get_position_tuple("BA000021.3_2")

            assert result == (2230, 3066, 1)

    def test_query_not_found(self, sample_faa_file):
        """Test lookup for non-existent protein."""
        db_path = sample_faa_file + ".sqlite"
        build_prodigal_index(sample_faa_file, db_path, show_progress=False)

        with ProdigalIndex(db_path) as idx:
            assert idx.get_position("nonexistent") is None
            assert idx.get_position_tuple("nonexistent") is None

    def test_batch_query(self, sample_faa_file):
        """Test batch lookup."""
        db_path = sample_faa_file + ".sqlite"
        build_prodigal_index(sample_faa_file, db_path, show_progress=False)

        with ProdigalIndex(db_path) as idx:
            ids = ["BA000021.3_1", "BA000021.3_3", "NC_000001.11_1", "nonexistent"]
            results = idx.get_positions(ids)

            assert len(results) == 3  # 3 found, 1 not found
            assert "BA000021.3_1" in results
            assert "BA000021.3_3" in results
            assert "NC_000001.11_1" in results
            assert "nonexistent" not in results

            # Check values
            assert results["BA000021.3_3"].strand == -1
            assert results["NC_000001.11_1"].start == 100

    def test_batch_query_tuples(self, sample_faa_file):
        """Test batch lookup with tuples."""
        db_path = sample_faa_file + ".sqlite"
        build_prodigal_index(sample_faa_file, db_path, show_progress=False)

        with ProdigalIndex(db_path) as idx:
            ids = ["BA000021.3_1", "BA000021.3_2"]
            results = idx.get_positions_tuples(ids)

            assert results["BA000021.3_1"] == (182, 2068, 1)
            assert results["BA000021.3_2"] == (2230, 3066, 1)

    def test_contains(self, sample_faa_file):
        """Test contains check."""
        db_path = sample_faa_file + ".sqlite"
        build_prodigal_index(sample_faa_file, db_path, show_progress=False)

        with ProdigalIndex(db_path) as idx:
            assert idx.contains("BA000021.3_1")
            assert "BA000021.3_1" in idx
            assert not idx.contains("nonexistent")
            assert "nonexistent" not in idx

    def test_len(self, sample_faa_file):
        """Test length."""
        db_path = sample_faa_file + ".sqlite"
        build_prodigal_index(sample_faa_file, db_path, show_progress=False)

        with ProdigalIndex(db_path) as idx:
            assert len(idx) == 5

    def test_get_by_contig(self, sample_faa_file):
        """Test getting all proteins from a contig."""
        db_path = sample_faa_file + ".sqlite"
        build_prodigal_index(sample_faa_file, db_path, show_progress=False)

        with ProdigalIndex(db_path) as idx:
            proteins = idx.get_proteins_by_contig("BA000021.3")

            assert len(proteins) == 3
            # Should be sorted by start position
            assert proteins[0].start < proteins[1].start < proteins[2].start

    def test_stats(self, sample_faa_file):
        """Test stats method."""
        db_path = sample_faa_file + ".sqlite"
        build_prodigal_index(sample_faa_file, db_path, show_progress=False)

        with ProdigalIndex(db_path) as idx:
            stats = idx.stats()

            assert stats['total_proteins'] == 5
            assert stats['total_contigs'] == 2
            assert stats['db_size_mb'] > 0


class TestGffIndex:
    """Test GFF-based index building."""

    @pytest.fixture
    def sample_gff_file(self):
        """Create a sample GFF file for testing."""
        content = """##gff-version 3
BA000021.3\tProdigal_v2.6.3\tCDS\t182\t2068\t372.0\t+\t0\tID=1_1;partial=00
BA000021.3\tProdigal_v2.6.3\tCDS\t2230\t3066\t104.1\t+\t0\tID=1_2;partial=00
BA000021.3\tProdigal_v2.6.3\tCDS\t3076\t3315\t49.1\t-\t0\tID=1_3;partial=00
"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.gff') as f:
            f.write(content)
            gff_path = f.name

        yield gff_path

        os.unlink(gff_path)
        db_path = gff_path + ".sqlite"
        if os.path.exists(db_path):
            os.unlink(db_path)

    def test_build_from_gff(self, sample_gff_file):
        """Test building index from GFF file."""
        db_path = sample_gff_file + ".sqlite"

        result = build_index_from_gff(
            sample_gff_file,
            db_path,
            show_progress=False
        )

        assert result == db_path
        assert os.path.exists(db_path)

        with ProdigalIndex(db_path) as idx:
            assert len(idx) == 3


class TestProteinPosition:
    """Test ProteinPosition dataclass."""

    def test_properties(self):
        """Test ProteinPosition properties."""
        pos = ProteinPosition(
            protein_id="test_1",
            start=100,
            end=500,
            strand=1,
            contig_id="test"
        )

        assert pos.length == 401
        assert pos.as_tuple() == (100, 500, 1)


class TestErrorHandling:
    """Test error handling."""

    def test_index_not_found(self):
        """Test opening non-existent index."""
        with pytest.raises(FileNotFoundError):
            ProdigalIndex("/nonexistent/path.sqlite")

    def test_faa_not_found(self):
        """Test building from non-existent FAA."""
        with pytest.raises(FileNotFoundError):
            build_prodigal_index("/nonexistent/file.faa", "/tmp/test.sqlite")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
