"""
Unit tests for the IS Element Processor module.

Tests cover:
- Sequence extraction from IS element dicts
- Flanking x noncoding alignment with provenance tracking
- Key mapping to FilterEngine format
- File loading and batch processing
- Edge cases (empty NC regions, missing flanking, etc.)

Author: Kuang Hu
Date: 2026-02-09
"""

import json
import os
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules.is_element_processor import (
    ISElementProcessor,
    load_is_elements,
    process_is_element,
    process_is_elements_file,
)
from modules.alignment_filter import create_default_pipeline


# ------------------------------------------------------------------
# Fixtures: synthetic IS element data
# ------------------------------------------------------------------

def _make_element(is_id="TEST_s0", sample="TEST",
                  is_length=200,
                  upstream_seq="", downstream_seq="",
                  noncoding_regions=None):
    """Build a minimal IS element dict for testing."""
    if noncoding_regions is None:
        noncoding_regions = []
    return {
        "is_id": is_id,
        "sample": sample,
        "seqid": "s0",
        "cluster": "c0",
        "group": "g1",
        "is_element": {
            "sequence": "A" * is_length,
            "length": is_length,
            "contig": "contig1",
            "start": 1000,
            "end": 1000 + is_length,
            "strand": "+",
        },
        "flanking_upstream": {
            "sequence": upstream_seq,
            "length": len(upstream_seq),
            "contig_start": 920,
            "contig_end": 1000,
        },
        "flanking_downstream": {
            "sequence": downstream_seq,
            "length": len(downstream_seq),
            "contig_start": 1000 + is_length,
            "contig_end": 1000 + is_length + len(downstream_seq),
        },
        "orf_annotation": {
            "is_id": is_id,
            "is_length": is_length,
            "num_orfs": 1,
            "num_noncoding_regions": len(noncoding_regions),
            "noncoding_regions": noncoding_regions,
        },
    }


# Shared motif: 12bp that will appear in both flanking and noncoding
SHARED_MOTIF = "TTTAAACCCGGG"


@pytest.fixture
def element_with_shared_motif():
    """Element where upstream flanking shares a 12bp motif with noncoding region 0."""
    upstream = "ATGCATGC" + SHARED_MOTIF + "GCTAGCTAGCTAGCTA"  # 36bp
    downstream = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"  # 35bp, no shared motif
    nc_regions = [
        {
            "start": 1,
            "end": 40,
            "length": 40,
            "type": "5_prime_utr",
            "sequence": "CCCCCCCCCCCC" + SHARED_MOTIF + "ACACACACACACACAC",  # 40bp
        },
        {
            "start": 100,
            "end": 130,
            "length": 31,
            "type": "intergenic",
            "sequence": "GTACGTACGTACGTACGTACGTACGTACGTG",  # 31bp, no shared motif
        },
    ]
    return _make_element(
        is_id="SAMPLE_s0", sample="SAMPLE", is_length=200,
        upstream_seq=upstream, downstream_seq=downstream,
        noncoding_regions=nc_regions,
    )


@pytest.fixture
def element_empty_noncoding():
    """Element with no noncoding regions (100% coding)."""
    return _make_element(
        is_id="CODING_s0", sample="CODING",
        upstream_seq="ACGTACGTACGTACGTACGT",
        downstream_seq="GCTAGCTAGCTAGCTAGCTA",
        noncoding_regions=[],
    )


@pytest.fixture
def element_missing_flanking():
    """Element with empty upstream flanking."""
    return _make_element(
        is_id="NOFLANKING_s0", sample="NOFLANKING",
        upstream_seq="",
        downstream_seq="GCTAGCTAGCTAGCTAGCTA",
        noncoding_regions=[{
            "start": 1, "end": 20, "length": 20,
            "type": "5_prime_utr",
            "sequence": "GCTAGCTAGCTAGCTAGCTA",
        }],
    )


@pytest.fixture
def json_file_two_elements(tmp_path):
    """Write a temp JSON file with two IS elements."""
    motif = "TTTAAACCCGGG"
    elements = [
        _make_element(
            is_id="S1_s0", sample="S1", is_length=150,
            upstream_seq="ATGCATGC" + motif + "GCTA",
            downstream_seq="CCCCCCCCCCCCCCCCCCCCCCCC",
            noncoding_regions=[{
                "start": 1, "end": 30, "length": 30,
                "type": "5_prime_utr",
                "sequence": "CCCCCCCC" + motif + "ACACACACAC",
            }],
        ),
        _make_element(
            is_id="S1_s1", sample="S1", is_length=300,
            upstream_seq="AAAAAAAAAAAAAAAAAAAAAAAA",
            downstream_seq="TTTTTTTTTTTTTTTTTTTTTTTT",
            noncoding_regions=[{
                "start": 50, "end": 80, "length": 31,
                "type": "intergenic",
                "sequence": "GTACGTACGTACGTACGTACGTACGTACGTG",
            }],
        ),
    ]
    path = tmp_path / "is_elements.json"
    path.write_text(json.dumps(elements))
    return str(path)


# ------------------------------------------------------------------
# Tests: _extract_sequences
# ------------------------------------------------------------------

class TestExtractSequences:

    def test_extracts_both_flanking_and_nc(self, element_with_shared_motif):
        proc = ISElementProcessor()
        flanking, nc = proc._extract_sequences(element_with_shared_motif)

        assert "upstream" in flanking
        assert "downstream" in flanking
        assert len(flanking["upstream"]) == 36
        assert len(flanking["downstream"]) == 35
        assert len(nc) == 2
        assert nc[0]["type"] == "5_prime_utr"
        assert nc[1]["type"] == "intergenic"

    def test_empty_noncoding(self, element_empty_noncoding):
        proc = ISElementProcessor()
        flanking, nc = proc._extract_sequences(element_empty_noncoding)

        assert len(nc) == 0
        assert flanking["upstream"] != ""

    def test_missing_flanking(self, element_missing_flanking):
        proc = ISElementProcessor()
        flanking, nc = proc._extract_sequences(element_missing_flanking)

        assert flanking["upstream"] == ""
        assert flanking["downstream"] != ""


# ------------------------------------------------------------------
# Tests: process_element core flow
# ------------------------------------------------------------------

class TestProcessElement:

    def test_finds_shared_motif(self, element_with_shared_motif):
        """Upstream flanking shares a 12bp motif with NC region 0."""
        proc = ISElementProcessor(min_length=9, max_mismatches=0,
                                  check_forward=True, check_revcomp=False)
        result = proc.process_element(element_with_shared_motif)

        assert result["is_id"] == "SAMPLE_s0"
        assert result["sample"] == "SAMPLE"
        assert result["num_noncoding_regions"] == 2
        assert result["total_hits"] >= 1

        # At least one hit should be the 12bp shared motif
        motif_hits = [a for a in result["alignments"]
                      if a["ungapped"]["length"] >= 12
                      and a["flanking_source"] == "upstream"
                      and a["noncoding_region_index"] == 0]
        assert len(motif_hits) >= 1

    def test_provenance_fields(self, element_with_shared_motif):
        proc = ISElementProcessor(min_length=9, max_mismatches=0,
                                  check_forward=True, check_revcomp=False)
        result = proc.process_element(element_with_shared_motif)

        for aln in result["alignments"]:
            assert "is_id" in aln
            assert "sample" in aln
            assert "flanking_source" in aln
            assert aln["flanking_source"] in ("upstream", "downstream")
            assert "noncoding_region_index" in aln
            assert "noncoding_region_type" in aln
            assert "noncoding_start" in aln
            assert "noncoding_end" in aln
            assert "noncoding_length" in aln
            assert "ungapped" in aln
            assert "gapped" in aln

    def test_ungapped_alignment_fields(self, element_with_shared_motif):
        proc = ISElementProcessor(min_length=9, max_mismatches=0,
                                  check_forward=True, check_revcomp=False)
        result = proc.process_element(element_with_shared_motif)

        for aln in result["alignments"]:
            ug = aln["ungapped"]
            assert "pos_in_flanking" in ug
            assert "pos_in_noncoding" in ug
            assert "length" in ug
            assert "seq_flanking" in ug
            assert "seq_noncoding" in ug
            assert "mismatches" in ug
            assert "mismatch_positions" in ug
            assert "orientation" in ug

    def test_empty_noncoding_produces_no_hits(self, element_empty_noncoding):
        proc = ISElementProcessor(min_length=9)
        result = proc.process_element(element_empty_noncoding)

        assert result["total_hits"] == 0
        assert result["alignments"] == []

    def test_missing_flanking_skips_direction(self, element_missing_flanking):
        proc = ISElementProcessor(min_length=9)
        result = proc.process_element(element_missing_flanking)

        # Only downstream should produce hits (upstream is empty)
        for aln in result["alignments"]:
            assert aln["flanking_source"] == "downstream"

    def test_both_directions_searched(self, element_with_shared_motif):
        """Both upstream and downstream flanking should be tested."""
        proc = ISElementProcessor(min_length=9, max_mismatches=0)
        result = proc.process_element(element_with_shared_motif)

        sources = {a["flanking_source"] for a in result["alignments"]}
        # upstream should have hits (shared motif); downstream might not
        assert "upstream" in sources

    def test_gapped_extension(self, element_with_shared_motif):
        proc = ISElementProcessor(min_length=9, max_mismatches=0,
                                  check_forward=True, check_revcomp=False,
                                  extend_with_gaps=True)
        result = proc.process_element(element_with_shared_motif)

        # At least one hit should have gapped extension attempted
        # (it may be None if identity is below threshold)
        for aln in result["alignments"]:
            assert "gapped" in aln  # key always present

    def test_result_structure(self, element_with_shared_motif):
        proc = ISElementProcessor(min_length=9)
        result = proc.process_element(element_with_shared_motif)

        assert "is_id" in result
        assert "sample" in result
        assert "is_element_length" in result
        assert "num_noncoding_regions" in result
        assert "total_hits" in result
        assert "alignments" in result
        assert isinstance(result["alignments"], list)
        assert result["total_hits"] == len(result["alignments"])


# ------------------------------------------------------------------
# Tests: filter pipeline integration
# ------------------------------------------------------------------

class TestFilterIntegration:

    def test_filter_pipeline_produces_filtered_key(self, element_with_shared_motif):
        pipeline = create_default_pipeline()
        proc = ISElementProcessor(min_length=9, max_mismatches=0,
                                  check_forward=True, check_revcomp=False,
                                  filter_pipeline=pipeline)
        result = proc.process_element(element_with_shared_motif)

        assert "filtered" in result
        filt = result["filtered"]
        assert "total" in filt
        assert "passed" in filt
        assert "failed" in filt
        assert "by_confidence" in filt
        assert "results" in filt
        assert filt["total"] == result["total_hits"]
        assert filt["passed"] + filt["failed"] == filt["total"]

    def test_filter_result_per_alignment(self, element_with_shared_motif):
        pipeline = create_default_pipeline()
        proc = ISElementProcessor(min_length=9, max_mismatches=0,
                                  check_forward=True, check_revcomp=False,
                                  filter_pipeline=pipeline)
        result = proc.process_element(element_with_shared_motif)

        for entry in result["filtered"]["results"]:
            assert "alignment" in entry
            assert "filter_result" in entry
            fr = entry["filter_result"]
            assert "pass" in fr
            assert "confidence" in fr

    def test_no_filter_when_pipeline_is_none(self, element_with_shared_motif):
        proc = ISElementProcessor(min_length=9, filter_pipeline=None)
        result = proc.process_element(element_with_shared_motif)

        assert "filtered" not in result


# ------------------------------------------------------------------
# Tests: _map_to_filter_format
# ------------------------------------------------------------------

class TestMapToFilterFormat:

    def test_ungapped_mapping(self):
        proc = ISElementProcessor()
        aln = {
            "flanking_source": "upstream",
            "flanking_length": 80,
            "noncoding_start": 417,
            "noncoding_end": 487,
            "noncoding_length": 71,
            "ungapped": {
                "pos_in_flanking": 5,
                "pos_in_noncoding": 12,
                "length": 14,
                "seq_flanking": "ACGTACGTACGTAC",
                "seq_noncoding": "ACGTACGTACGTAC",
                "mismatches": 0,
                "mismatch_positions": [],
                "orientation": "forward",
            },
            "gapped": None,
        }
        mapped = proc._map_to_filter_format(aln)

        assert mapped["aligned_sequence"] == "ACGTACGTACGTAC"
        assert mapped["non_coding_start"] == 417
        assert mapped["non_coding_end"] == 487
        assert mapped["mismatches"] == 0
        assert mapped["gaps"] == 0
        assert mapped["alignment_string"] == ""
        assert mapped["query_length"] == 80
        assert mapped["target_length"] == 71

    def test_gapped_mapping(self):
        proc = ISElementProcessor()
        aln = {
            "flanking_source": "upstream",
            "flanking_length": 80,
            "noncoding_start": 100,
            "noncoding_end": 150,
            "noncoding_length": 51,
            "ungapped": {
                "pos_in_flanking": 5,
                "pos_in_noncoding": 12,
                "length": 14,
                "seq_flanking": "ACGTACGTACGTAC",
                "seq_noncoding": "ACGTACGTACGTAC",
                "mismatches": 0,
                "mismatch_positions": [],
                "orientation": "forward",
            },
            "gapped": {
                "pos_in_flanking": 3,
                "pos_in_noncoding": 10,
                "end_in_flanking": 20,
                "end_in_noncoding": 28,
                "length_in_flanking": 17,
                "length_in_noncoding": 18,
                "seq1_aligned": "ACGT-ACGTACGTAC",
                "seq2_aligned": "ACGTAACGTACGTAC",
                "alignment_length": 15,
                "matches": 14,
                "mismatches": 0,
                "gaps": 1,
                "identity": 0.933,
                "alignment_string": "|||| ||||||||||",
                "orientation": "forward",
            },
        }
        mapped = proc._map_to_filter_format(aln)

        # aligned_sequence should be the flanking portion without gap chars
        assert mapped["aligned_sequence"] == "ACGTACGTACGTAC"
        assert mapped["mismatches"] == 0
        assert mapped["gaps"] == 1
        assert mapped["alignment_string"] == "|||| ||||||||||"


# ------------------------------------------------------------------
# Tests: _build_transposon_data
# ------------------------------------------------------------------

class TestBuildTransposonData:

    def test_basic(self):
        proc = ISElementProcessor()
        element = {"is_element": {"length": 1527}}
        td = proc._build_transposon_data(element)

        assert td["start"] == 1
        assert td["end"] == 1527

    def test_missing_is_element(self):
        proc = ISElementProcessor()
        td = proc._build_transposon_data({})

        assert td["start"] == 1
        assert td["end"] == 0


# ------------------------------------------------------------------
# Tests: file I/O
# ------------------------------------------------------------------

class TestFileIO:

    def test_load_is_elements(self, json_file_two_elements):
        elements = load_is_elements(json_file_two_elements)

        assert len(elements) == 2
        assert elements[0]["is_id"] == "S1_s0"
        assert elements[1]["is_id"] == "S1_s1"

    def test_load_is_elements_invalid_format(self, tmp_path):
        bad_path = tmp_path / "bad.json"
        bad_path.write_text(json.dumps({"not": "a list"}))

        with pytest.raises(ValueError, match="Expected JSON array"):
            load_is_elements(str(bad_path))

    def test_process_file(self, json_file_two_elements):
        proc = ISElementProcessor(min_length=9, check_forward=True,
                                  check_revcomp=False)
        results = proc.process_file(json_file_two_elements)

        assert len(results) == 2
        assert results[0]["is_id"] == "S1_s0"
        assert results[1]["is_id"] == "S1_s1"
        # First element should have hits (shared motif)
        assert results[0]["total_hits"] >= 1

    def test_process_batch(self, json_file_two_elements, tmp_path):
        # Create a second file with one element
        elements = [_make_element(
            is_id="S2_s0", sample="S2", is_length=100,
            upstream_seq="ACGTACGTACGTACGTACGT",
            downstream_seq="GCTAGCTAGCTAGCTAGCTA",
            noncoding_regions=[],
        )]
        path2 = tmp_path / "is_elements2.json"
        path2.write_text(json.dumps(elements))

        proc = ISElementProcessor(min_length=9)
        results = proc.process_batch([json_file_two_elements, str(path2)])

        assert len(results) == 3  # 2 from file1 + 1 from file2


# ------------------------------------------------------------------
# Tests: convenience functions
# ------------------------------------------------------------------

class TestConvenienceFunctions:

    def test_process_is_element(self, element_with_shared_motif):
        result = process_is_element(
            element_with_shared_motif,
            min_length=9, max_mismatches=0,
            check_forward=True, check_revcomp=False,
        )
        assert result["is_id"] == "SAMPLE_s0"
        assert result["total_hits"] >= 1

    def test_process_is_elements_file(self, json_file_two_elements):
        results = process_is_elements_file(
            json_file_two_elements,
            min_length=9, check_forward=True, check_revcomp=False,
        )
        assert len(results) == 2


# ------------------------------------------------------------------
# Tests: edge cases
# ------------------------------------------------------------------

class TestEdgeCases:

    def test_nc_region_shorter_than_min_length(self):
        """NC region shorter than min_length should produce no hits."""
        element = _make_element(
            upstream_seq="ACGTACGTACGTACGTACGT",
            noncoding_regions=[{
                "start": 1, "end": 5, "length": 5,
                "type": "5_prime_utr",
                "sequence": "ACGTA",  # Only 5bp
            }],
        )
        proc = ISElementProcessor(min_length=9)
        result = proc.process_element(element)

        assert result["total_hits"] == 0

    def test_noncoding_type_passthrough(self):
        """Noncoding type field should be passed through as-is."""
        for nc_type in ["5_prime_utr", "3_prime_utr", "intergenic", "5prime_utr"]:
            element = _make_element(
                upstream_seq="ATGCATGC" + SHARED_MOTIF + "GCTA",
                noncoding_regions=[{
                    "start": 1, "end": 30, "length": 30,
                    "type": nc_type,
                    "sequence": "CCCCCCCC" + SHARED_MOTIF + "ACACACACAC",
                }],
            )
            proc = ISElementProcessor(min_length=9, check_forward=True,
                                      check_revcomp=False)
            result = proc.process_element(element)

            if result["total_hits"] > 0:
                assert result["alignments"][0]["noncoding_region_type"] == nc_type

    def test_multiple_nc_regions_indexed_correctly(self):
        """Each noncoding region should get its correct 0-based index."""
        element = _make_element(
            upstream_seq="ATGCATGC" + SHARED_MOTIF + "GCTA",
            noncoding_regions=[
                {
                    "start": 1, "end": 20, "length": 20,
                    "type": "5_prime_utr",
                    "sequence": "GTACGTACGTACGTACGTAC",  # no shared motif
                },
                {
                    "start": 80, "end": 110, "length": 30,
                    "type": "intergenic",
                    "sequence": "CCCCCCCC" + SHARED_MOTIF + "ACACACACAC",
                },
            ],
        )
        proc = ISElementProcessor(min_length=9, check_forward=True,
                                  check_revcomp=False)
        result = proc.process_element(element)

        # Hits in NC region 1 (the one with the shared motif)
        nc1_hits = [a for a in result["alignments"]
                    if a["noncoding_region_index"] == 1]
        assert len(nc1_hits) >= 1
        assert nc1_hits[0]["noncoding_region_type"] == "intergenic"
        assert nc1_hits[0]["noncoding_start"] == 80

    def test_revcomp_orientation_tracked(self):
        """Reverse complement matches should have correct orientation."""
        from utils.parsers import reverse_complement
        motif = "ACGTACGTACGT"
        rc_motif = reverse_complement(motif)

        element = _make_element(
            upstream_seq="GGGG" + motif + "GGGG",
            noncoding_regions=[{
                "start": 1, "end": 24, "length": 24,
                "type": "5_prime_utr",
                "sequence": "CCCCCC" + rc_motif + "CCCCCC",
            }],
        )
        proc = ISElementProcessor(min_length=9, check_forward=False,
                                  check_revcomp=True)
        result = proc.process_element(element)

        rc_hits = [a for a in result["alignments"]
                   if a["ungapped"]["orientation"] == "reverse_complement"]
        assert len(rc_hits) >= 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
