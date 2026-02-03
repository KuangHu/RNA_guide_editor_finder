"""
Unit tests for parsers module

Tests all parser functions including:
- Prodigal parsers
- FASTA parsers
- GFF parsers
- BED parsers
- Diamond BLASTP parsers
- HMMER parsers
- Sequence retrieval functions

Run with: python -m pytest tests/test_parsers.py -v
"""

import os
import sys
import tempfile
import pytest
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.parsers import (
    # Data classes
    ProdigalGene,
    FastaRecord,
    GffFeature,
    BedRecord,
    DiamondBlastHit,
    HmmHit,
    # Prodigal parsers
    parse_prodigal_faa,
    parse_prodigal_faa_full,
    parse_prodigal_gff,
    # FASTA parsers
    parse_fasta,
    parse_fasta_records,
    iter_fasta,
    # GFF parsers
    parse_gff,
    parse_gff_line,
    iter_gff,
    # BED parsers
    parse_bed,
    iter_bed,
    # Diamond BLASTP parsers
    parse_diamond_blastp,
    parse_diamond_blastp_with_positions,
    iter_diamond_blastp_with_positions,
    # HMMER parsers
    parse_hmm_tblout,
    parse_hmm_tblout_with_positions,
    iter_hmm_tblout_with_positions,
    # Sequence retrieval
    retrieve_sequence_from_fasta,
    retrieve_sequences_from_hits,
    reverse_complement,
    # Utilities
    detect_file_format,
)


# ============================================
# Test Fixtures - Sample Files
# ============================================

@pytest.fixture
def sample_fasta_file(tmp_path):
    """Create a sample FASTA file"""
    fasta_content = """>contig_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>contig_2
AAAATTTTCCCCGGGGAAAATTTTCCCCGGGG
"""
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)
    return str(fasta_file)


@pytest.fixture
def sample_prodigal_faa_file(tmp_path):
    """Create a sample Prodigal FAA file"""
    faa_content = """>contig_1_1 # 100 # 300 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=None
MKLVPQRSTAVILGKLMNPQ
>contig_1_2 # 500 # 800 # -1 # ID=1_2;partial=00;start_type=ATG;rbs_motif=GGA
MLKIPVRSTQAVILGKLMN
>contig_2_1 # 50 # 250 # 1 # ID=2_1;partial=00;start_type=ATG;rbs_motif=None
MKLVPQRSTAVILGKLMNPQRST
"""
    faa_file = tmp_path / "test.faa"
    faa_file.write_text(faa_content)
    return str(faa_file)


@pytest.fixture
def sample_gff_file(tmp_path):
    """Create a sample GFF file"""
    gff_content = """##gff-version 3
contig_1	Prodigal	CDS	100	300	.	+	0	ID=contig_1_1;product=hypothetical protein
contig_1	Prodigal	CDS	500	800	.	-	0	ID=contig_1_2;product=transposase
contig_2	Prodigal	CDS	50	250	.	+	0	ID=contig_2_1;product=DNA polymerase
"""
    gff_file = tmp_path / "test.gff"
    gff_file.write_text(gff_content)
    return str(gff_file)


@pytest.fixture
def sample_bed_file(tmp_path):
    """Create a sample BED file"""
    bed_content = """contig_1	99	300	gene1	100	+
contig_1	499	800	gene2	200	-
contig_2	49	250	gene3	150	+
"""
    bed_file = tmp_path / "test.bed"
    bed_file.write_text(bed_content)
    return str(bed_file)


@pytest.fixture
def sample_diamond_blastp_file(tmp_path):
    """Create a sample Diamond BLASTP output file"""
    blastp_content = """contig_1_1	IS110_transposase	95.5	200	9	0	1	200	15	214	1e-100	350.5
contig_1_2	IS110_orfB	88.2	150	18	1	1	150	10	159	1e-75	280.3
contig_2_1	DNA_polymerase	92.0	180	14	0	1	180	5	184	1e-90	320.1
"""
    blastp_file = tmp_path / "test_blastp.m8"
    blastp_file.write_text(blastp_content)
    return str(blastp_file)


@pytest.fixture
def sample_hmm_tblout_file(tmp_path):
    """Create a sample HMMER tblout file"""
    hmm_content = """#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
contig_1_1           -            200 IS110_transposase    PF03400.15   250  1.5e-100  350.5   0.1   1   1  1.5e-100  1.5e-100  350.5   0.1     1   250     1   200     1   200 0.98 -
contig_1_2           -            150 IS110_orfB           PF12345.10   180  2.3e-75   280.3   0.2   1   1  2.3e-75   2.3e-75   280.3   0.2     1   180     1   150     1   150 0.95 -
contig_2_1           -            180 DNA_polymerase       PF00136.22   200  5.1e-90   320.1   0.0   1   1  5.1e-90   5.1e-90   320.1   0.0     1   200     1   180     1   180 0.99 -
"""
    hmm_file = tmp_path / "test_hmmer.tbl"
    hmm_file.write_text(hmm_content)
    return str(hmm_file)


@pytest.fixture
def sample_position_map():
    """Create a sample position map for testing"""
    return {
        'contig_1_1': (100, 300, 1),
        'contig_1_2': (500, 800, -1),
        'contig_2_1': (50, 250, 1),
    }


# ============================================
# Test Prodigal Parsers
# ============================================

class TestProdigalParsers:
    """Test Prodigal file parsers"""

    def test_parse_prodigal_faa(self, sample_prodigal_faa_file):
        """Test parsing Prodigal FAA file"""
        positions = parse_prodigal_faa(sample_prodigal_faa_file)

        assert len(positions) == 3
        assert 'contig_1_1' in positions
        assert positions['contig_1_1'] == (100, 300, 1)
        assert positions['contig_1_2'] == (500, 800, -1)
        assert positions['contig_2_1'] == (50, 250, 1)

    def test_parse_prodigal_faa_full(self, sample_prodigal_faa_file):
        """Test parsing Prodigal FAA file with full information"""
        genes = parse_prodigal_faa_full(sample_prodigal_faa_file)

        assert len(genes) == 3
        assert 'contig_1_1' in genes

        gene = genes['contig_1_1']
        assert isinstance(gene, ProdigalGene)
        assert gene.start == 100
        assert gene.end == 300
        assert gene.strand == 1
        assert gene.is_forward
        assert gene.length == 201
        assert gene.sequence == "MKLVPQRSTAVILGKLMNPQ"

    def test_parse_prodigal_gff(self, sample_gff_file):
        """Test parsing Prodigal GFF file"""
        features = parse_prodigal_gff(sample_gff_file)

        assert len(features) == 3
        assert all(isinstance(f, GffFeature) for f in features)

        feature = features[0]
        assert feature.seqid == "contig_1"
        assert feature.feature_type == "CDS"
        assert feature.start == 100
        assert feature.end == 300
        assert feature.strand == "+"


# ============================================
# Test FASTA Parsers
# ============================================

class TestFastaParsers:
    """Test FASTA file parsers"""

    def test_parse_fasta(self, sample_fasta_file):
        """Test parsing FASTA file"""
        sequences = parse_fasta(sample_fasta_file)

        assert len(sequences) == 2
        assert 'contig_1' in sequences
        assert 'contig_2' in sequences
        assert sequences['contig_1'].startswith('ATCGATCG')
        assert len(sequences['contig_2']) == 32

    def test_parse_fasta_records(self, sample_fasta_file):
        """Test parsing FASTA file into records"""
        records = parse_fasta_records(sample_fasta_file)

        assert len(records) == 2
        assert all(isinstance(r, FastaRecord) for r in records)

        record = records[0]
        assert record.id == "contig_1"
        assert record.header == "contig_1"
        assert len(record.sequence) > 0
        assert record.length == len(record.sequence)

    def test_iter_fasta(self, sample_fasta_file):
        """Test iterating over FASTA records"""
        records = list(iter_fasta(sample_fasta_file))

        assert len(records) == 2
        assert records[0].id == "contig_1"
        assert records[1].id == "contig_2"


# ============================================
# Test GFF Parsers
# ============================================

class TestGffParsers:
    """Test GFF file parsers"""

    def test_parse_gff_line(self):
        """Test parsing a single GFF line"""
        line = "contig_1\tProdigal\tCDS\t100\t300\t.\t+\t0\tID=gene1;product=transposase"
        feature = parse_gff_line(line)

        assert feature is not None
        assert feature.seqid == "contig_1"
        assert feature.source == "Prodigal"
        assert feature.feature_type == "CDS"
        assert feature.start == 100
        assert feature.end == 300
        assert feature.strand == "+"
        assert feature.length == 201
        assert 'ID' in feature.attributes
        assert feature.attributes['ID'] == 'gene1'

    def test_parse_gff(self, sample_gff_file):
        """Test parsing GFF file"""
        features = parse_gff(sample_gff_file)

        assert len(features) == 3
        assert all(isinstance(f, GffFeature) for f in features)

    def test_iter_gff(self, sample_gff_file):
        """Test iterating over GFF features"""
        features = list(iter_gff(sample_gff_file))

        assert len(features) == 3
        assert features[0].seqid == "contig_1"


# ============================================
# Test BED Parsers
# ============================================

class TestBedParsers:
    """Test BED file parsers"""

    def test_parse_bed(self, sample_bed_file):
        """Test parsing BED file"""
        records = parse_bed(sample_bed_file)

        assert len(records) == 3
        assert all(isinstance(r, BedRecord) for r in records)

        record = records[0]
        assert record.chrom == "contig_1"
        assert record.start == 99
        assert record.end == 300
        assert record.name == "gene1"
        assert record.score == 100
        assert record.strand == "+"
        assert record.length == 201

    def test_iter_bed(self, sample_bed_file):
        """Test iterating over BED records"""
        records = list(iter_bed(sample_bed_file))

        assert len(records) == 3
        assert records[0].chrom == "contig_1"


# ============================================
# Test Diamond BLASTP Parsers
# ============================================

class TestDiamondBlastpParsers:
    """Test Diamond BLASTP parsers"""

    def test_parse_diamond_blastp(self, sample_diamond_blastp_file):
        """Test parsing Diamond BLASTP output"""
        hits = parse_diamond_blastp(sample_diamond_blastp_file)

        assert len(hits) == 3
        assert hits[0][0] == "contig_1_1"  # query_id
        assert hits[0][1] == "IS110_transposase"  # subject_id
        assert hits[0][2] == 95.5  # pident
        assert hits[0][10] == 1e-100  # evalue

    def test_parse_diamond_blastp_with_positions(
        self, sample_diamond_blastp_file, sample_position_map
    ):
        """Test parsing Diamond BLASTP with genomic positions"""
        hits = parse_diamond_blastp_with_positions(
            sample_diamond_blastp_file,
            sample_position_map,
            file_basename="test_genome"
        )

        assert len(hits) == 3
        assert all(isinstance(h, DiamondBlastHit) for h in hits)

        hit = hits[0]
        assert hit.file_basename == "test_genome"
        assert hit.contig_id == "contig_1"
        assert hit.protein_id == "contig_1_1"
        assert hit.start == 100
        assert hit.end == 300
        assert hit.strand == 1
        assert hit.is_forward
        assert hit.length == 201
        assert hit.query_id == "contig_1_1"
        assert hit.subject_id == "IS110_transposase"
        assert hit.pident == 95.5
        assert hit.evalue == 1e-100

    def test_iter_diamond_blastp_with_positions(
        self, sample_diamond_blastp_file, sample_position_map
    ):
        """Test iterating over Diamond BLASTP hits"""
        hits = list(iter_diamond_blastp_with_positions(
            sample_diamond_blastp_file,
            sample_position_map
        ))

        assert len(hits) == 3
        assert all(isinstance(h, DiamondBlastHit) for h in hits)


# ============================================
# Test HMMER Parsers
# ============================================

class TestHmmerParsers:
    """Test HMMER parsers"""

    def test_parse_hmm_tblout(self, sample_hmm_tblout_file):
        """Test parsing HMMER tblout output"""
        hits = parse_hmm_tblout(sample_hmm_tblout_file)

        assert len(hits) == 3
        assert hits[0][0] == "contig_1_1"  # target_name
        assert hits[0][1] == "IS110_transposase"  # query_name
        assert hits[0][3] == 1.5e-100  # evalue

    def test_parse_hmm_tblout_with_positions(
        self, sample_hmm_tblout_file, sample_position_map
    ):
        """Test parsing HMMER tblout with genomic positions"""
        hits = parse_hmm_tblout_with_positions(
            sample_hmm_tblout_file,
            sample_position_map,
            file_basename="test_genome"
        )

        assert len(hits) == 3
        assert all(isinstance(h, HmmHit) for h in hits)

        hit = hits[0]
        assert hit.file_basename == "test_genome"
        assert hit.contig_id == "contig_1"
        assert hit.protein_id == "contig_1_1"
        assert hit.start == 100
        assert hit.end == 300
        assert hit.strand == 1
        assert hit.is_forward
        assert hit.length == 201
        assert hit.query_name == "IS110_transposase"
        assert hit.target_name == "contig_1_1"
        assert hit.evalue == 1.5e-100

    def test_iter_hmm_tblout_with_positions(
        self, sample_hmm_tblout_file, sample_position_map
    ):
        """Test iterating over HMMER hits"""
        hits = list(iter_hmm_tblout_with_positions(
            sample_hmm_tblout_file,
            sample_position_map
        ))

        assert len(hits) == 3
        assert all(isinstance(h, HmmHit) for h in hits)


# ============================================
# Test Sequence Retrieval
# ============================================

class TestSequenceRetrieval:
    """Test sequence retrieval functions"""

    def test_reverse_complement(self):
        """Test reverse complement function"""
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("GCTAGCTA") == "TAGCTAGC"
        assert reverse_complement("atcg") == "cgat"  # lowercase

    def test_retrieve_sequence_from_fasta_forward(self, sample_fasta_file):
        """Test retrieving sequence from FASTA on forward strand"""
        seq = retrieve_sequence_from_fasta(
            sample_fasta_file,
            "contig_1",
            start=1,
            end=10,
            strand=1
        )

        assert seq is not None
        assert len(seq) == 10
        assert seq == "ATCGATCGAT"

    def test_retrieve_sequence_from_fasta_reverse(self, sample_fasta_file):
        """Test retrieving sequence from FASTA on reverse strand"""
        seq = retrieve_sequence_from_fasta(
            sample_fasta_file,
            "contig_1",
            start=1,
            end=10,
            strand=-1
        )

        assert seq is not None
        assert len(seq) == 10
        # Should be reverse complement of ATCGATCGAT
        assert seq == "ATCGATCGAT"  # Palindrome

    def test_retrieve_sequence_not_found(self, sample_fasta_file):
        """Test retrieving sequence for non-existent contig"""
        seq = retrieve_sequence_from_fasta(
            sample_fasta_file,
            "nonexistent_contig",
            start=1,
            end=10,
            strand=1
        )

        assert seq is None

    def test_retrieve_sequences_from_diamond_hits(
        self, sample_diamond_blastp_file, sample_position_map, sample_fasta_file
    ):
        """Test retrieving sequences for Diamond BLASTP hits"""
        # Parse Diamond hits
        hits = parse_diamond_blastp_with_positions(
            sample_diamond_blastp_file,
            sample_position_map,
            file_basename="test_genome"
        )

        # Retrieve sequences
        hits_with_seq = retrieve_sequences_from_hits(
            hits,
            [sample_fasta_file],
            verbose=False
        )

        assert len(hits_with_seq) > 0
        # At least some hits should have sequences
        # (depends on whether contig IDs match between files)

    def test_retrieve_sequences_from_hmm_hits(
        self, sample_hmm_tblout_file, sample_position_map, sample_fasta_file
    ):
        """Test retrieving sequences for HMMER hits"""
        # Parse HMMER hits
        hits = parse_hmm_tblout_with_positions(
            sample_hmm_tblout_file,
            sample_position_map,
            file_basename="test_genome"
        )

        # Retrieve sequences
        hits_with_seq = retrieve_sequences_from_hits(
            hits,
            [sample_fasta_file],
            verbose=False
        )

        assert len(hits_with_seq) > 0


# ============================================
# Test Utilities
# ============================================

class TestUtilities:
    """Test utility functions"""

    def test_detect_file_format_fasta(self, sample_fasta_file):
        """Test detecting FASTA format"""
        format_type = detect_file_format(sample_fasta_file)
        assert format_type == 'fasta'

    def test_detect_file_format_gff(self, sample_gff_file):
        """Test detecting GFF format"""
        format_type = detect_file_format(sample_gff_file)
        assert format_type == 'gff'

    def test_detect_file_format_bed(self, sample_bed_file):
        """Test detecting BED format"""
        format_type = detect_file_format(sample_bed_file)
        assert format_type == 'bed'


# ============================================
# Integration Tests
# ============================================

class TestParserIntegration:
    """Integration tests for parser workflows"""

    def test_complete_diamond_workflow(
        self, sample_prodigal_faa_file, sample_diamond_blastp_file, sample_fasta_file
    ):
        """Test complete workflow: Prodigal -> Diamond -> Sequence retrieval"""

        # Step 1: Parse Prodigal positions
        position_map = parse_prodigal_faa(sample_prodigal_faa_file)
        assert len(position_map) == 3

        # Step 2: Parse Diamond BLASTP with positions
        hits = parse_diamond_blastp_with_positions(
            sample_diamond_blastp_file,
            position_map,
            file_basename="test_genome"
        )
        assert len(hits) == 3
        assert all(h.start > 0 for h in hits)

        # Step 3: Retrieve sequences
        hits_with_seq = retrieve_sequences_from_hits(
            hits,
            [sample_fasta_file],
            verbose=False
        )
        assert len(hits_with_seq) == 3

    def test_complete_hmmer_workflow(
        self, sample_prodigal_faa_file, sample_hmm_tblout_file, sample_fasta_file
    ):
        """Test complete workflow: Prodigal -> HMMER -> Sequence retrieval"""

        # Step 1: Parse Prodigal positions
        position_map = parse_prodigal_faa(sample_prodigal_faa_file)

        # Step 2: Parse HMMER with positions
        hits = parse_hmm_tblout_with_positions(
            sample_hmm_tblout_file,
            position_map,
            file_basename="test_genome"
        )
        assert len(hits) == 3

        # Step 3: Retrieve sequences
        hits_with_seq = retrieve_sequences_from_hits(
            hits,
            [sample_fasta_file],
            verbose=False
        )
        assert len(hits_with_seq) == 3


class TestDataFrameConversion:
    """Test DataFrame conversion functions"""

    def test_diamond_hits_to_dataframe(self, sample_diamond_blastp_file, sample_position_map):
        """Test converting Diamond hits to DataFrame"""
        try:
            import pandas as pd
        except ImportError:
            pytest.skip("pandas not installed")

        from utils.parsers import diamond_hits_to_dataframe

        hits = parse_diamond_blastp_with_positions(
            sample_diamond_blastp_file,
            sample_position_map,
            file_basename="test_genome"
        )

        df = diamond_hits_to_dataframe(hits)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        assert 'file_basename' in df.columns
        assert 'contig_id' in df.columns
        assert 'protein_id' in df.columns
        assert 'start' in df.columns
        assert 'end' in df.columns
        assert 'strand' in df.columns
        assert 'pident' in df.columns
        assert 'evalue' in df.columns
        assert 'length' in df.columns
        assert 'is_forward' in df.columns

        # Check data types
        assert df['start'].dtype == 'int64'
        assert df['strand'].dtype == 'int8'
        assert df['pident'].dtype == 'float64'
        assert df['is_forward'].dtype == 'bool'

    def test_hmm_hits_to_dataframe(self, sample_hmm_tblout_file, sample_position_map):
        """Test converting HMMER hits to DataFrame"""
        try:
            import pandas as pd
        except ImportError:
            pytest.skip("pandas not installed")

        from utils.parsers import hmm_hits_to_dataframe

        hits = parse_hmm_tblout_with_positions(
            sample_hmm_tblout_file,
            sample_position_map,
            file_basename="test_genome"
        )

        df = hmm_hits_to_dataframe(hits)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        assert 'file_basename' in df.columns
        assert 'contig_id' in df.columns
        assert 'protein_id' in df.columns
        assert 'query_name' in df.columns
        assert 'evalue' in df.columns
        assert 'score' in df.columns

    def test_empty_hits_to_dataframe(self):
        """Test converting empty hit list to DataFrame"""
        try:
            import pandas as pd
        except ImportError:
            pytest.skip("pandas not installed")

        from utils.parsers import diamond_hits_to_dataframe, hmm_hits_to_dataframe

        # Empty Diamond hits
        df_diamond = diamond_hits_to_dataframe([])
        assert isinstance(df_diamond, pd.DataFrame)
        assert len(df_diamond) == 0
        assert 'file_basename' in df_diamond.columns

        # Empty HMMER hits
        df_hmmer = hmm_hits_to_dataframe([])
        assert isinstance(df_hmmer, pd.DataFrame)
        assert len(df_hmmer) == 0
        assert 'file_basename' in df_hmmer.columns

    def test_dataframe_filtering(self, sample_diamond_blastp_file, sample_position_map):
        """Test filtering DataFrame"""
        try:
            import pandas as pd
        except ImportError:
            pytest.skip("pandas not installed")

        from utils.parsers import diamond_hits_to_dataframe

        hits = parse_diamond_blastp_with_positions(
            sample_diamond_blastp_file,
            sample_position_map,
            file_basename="test_genome"
        )

        df = diamond_hits_to_dataframe(hits)

        # Test filtering
        high_identity = df[df['pident'] > 90]
        assert len(high_identity) >= 0

        forward_strand = df[df['is_forward']]
        assert all(forward_strand['strand'] == 1)

    def test_save_hits_to_csv(self, sample_diamond_blastp_file, sample_position_map, tmp_path):
        """Test saving hits to CSV"""
        try:
            import pandas as pd
        except ImportError:
            pytest.skip("pandas not installed")

        from utils.parsers import save_hits_to_csv

        hits = parse_diamond_blastp_with_positions(
            sample_diamond_blastp_file,
            sample_position_map,
            file_basename="test_genome"
        )

        output_file = tmp_path / "test_output.csv"
        save_hits_to_csv(hits, str(output_file), include_sequence=False)

        assert output_file.exists()

        # Read back and verify
        df = pd.read_csv(output_file)
        assert len(df) == 3
        assert 'file_basename' in df.columns


# ============================================
# Test ProdigalIndex Integration
# ============================================

class TestProdigalIndexIntegration:
    """Test that parsers work with ProdigalIndex as position source"""

    @pytest.fixture
    def sample_faa_file(self, tmp_path):
        """Create a sample FAA file for building ProdigalIndex"""
        content = """>contig_1_1 # 100 # 300 # 1 # ID=1_1;partial=00
MTKQVLAAAA
>contig_1_2 # 500 # 800 # -1 # ID=1_2;partial=00
MKKQVLBBBB
>contig_2_1 # 50 # 250 # 1 # ID=1_3;partial=00
MLLQVLCCCC
"""
        faa_file = tmp_path / "test.faa"
        faa_file.write_text(content)
        return str(faa_file)

    @pytest.fixture
    def prodigal_index(self, sample_faa_file, tmp_path):
        """Build a ProdigalIndex for testing"""
        from utils.prodigal_index import ProdigalIndex, build_prodigal_index
        db_path = str(tmp_path / "test.sqlite")
        build_prodigal_index(sample_faa_file, db_path, show_progress=False)
        return ProdigalIndex(db_path)

    def test_diamond_with_prodigal_index(
        self, sample_diamond_blastp_file, prodigal_index
    ):
        """Test parsing Diamond BLASTP with ProdigalIndex"""
        hits = parse_diamond_blastp_with_positions(
            sample_diamond_blastp_file,
            prodigal_index,
            file_basename="test_genome"
        )

        assert len(hits) == 3
        assert all(isinstance(h, DiamondBlastHit) for h in hits)

        hit = hits[0]
        assert hit.protein_id == "contig_1_1"
        assert hit.start == 100
        assert hit.end == 300
        assert hit.strand == 1

    def test_iter_diamond_with_prodigal_index(
        self, sample_diamond_blastp_file, prodigal_index
    ):
        """Test iterating Diamond BLASTP with ProdigalIndex"""
        hits = list(iter_diamond_blastp_with_positions(
            sample_diamond_blastp_file,
            prodigal_index
        ))

        assert len(hits) == 3
        assert all(isinstance(h, DiamondBlastHit) for h in hits)

    def test_hmmer_with_prodigal_index(
        self, sample_hmm_tblout_file, prodigal_index
    ):
        """Test parsing HMMER tblout with ProdigalIndex"""
        hits = parse_hmm_tblout_with_positions(
            sample_hmm_tblout_file,
            prodigal_index,
            file_basename="test_genome"
        )

        assert len(hits) == 3
        assert all(isinstance(h, HmmHit) for h in hits)

        hit = hits[0]
        assert hit.protein_id == "contig_1_1"
        assert hit.start == 100
        assert hit.end == 300
        assert hit.strand == 1

    def test_iter_hmmer_with_prodigal_index(
        self, sample_hmm_tblout_file, prodigal_index
    ):
        """Test iterating HMMER with ProdigalIndex"""
        hits = list(iter_hmm_tblout_with_positions(
            sample_hmm_tblout_file,
            prodigal_index
        ))

        assert len(hits) == 3
        assert all(isinstance(h, HmmHit) for h in hits)

    def test_dict_and_index_produce_same_results(
        self, sample_diamond_blastp_file, sample_position_map, prodigal_index
    ):
        """Verify dict and ProdigalIndex produce identical results"""
        hits_dict = parse_diamond_blastp_with_positions(
            sample_diamond_blastp_file,
            sample_position_map
        )
        hits_index = parse_diamond_blastp_with_positions(
            sample_diamond_blastp_file,
            prodigal_index
        )

        assert len(hits_dict) == len(hits_index)
        for h_dict, h_idx in zip(hits_dict, hits_index):
            assert h_dict.protein_id == h_idx.protein_id
            assert h_dict.start == h_idx.start
            assert h_dict.end == h_idx.end
            assert h_dict.strand == h_idx.strand


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
