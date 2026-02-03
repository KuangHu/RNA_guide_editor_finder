"""
Parsers Module for RNA Guide Editor Finder

A collection of parsers for common bioinformatics file formats:
- Prodigal output (.faa, .gff)
- FASTA files
- GFF/GFF3 files
- BED files

Usage:
    from utils.parsers import parse_prodigal_faa, parse_fasta, parse_gff
"""

from typing import Dict, List, Tuple, Optional, Iterator, Union

# Import ProdigalIndex for type hints (lazy import to avoid circular dependency)
# The actual class is imported within functions that need it
try:
    from .prodigal_index import ProdigalIndex
    PRODIGAL_INDEX_AVAILABLE = True
except ImportError:
    PRODIGAL_INDEX_AVAILABLE = False
    ProdigalIndex = None  # type: ignore

# Type alias for position source: either a dict or ProdigalIndex
PositionSource = Union[Dict[str, Tuple[int, int, int]], 'ProdigalIndex']
from dataclasses import dataclass, asdict
from pathlib import Path

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False


# ============================================
# Data Classes
# ============================================

@dataclass
class ProdigalGene:
    """Represents a gene predicted by Prodigal"""
    protein_id: str
    start: int
    end: int
    strand: int  # 1 for forward, -1 for reverse
    sequence: Optional[str] = None

    @property
    def length(self) -> int:
        return abs(self.end - self.start) + 1

    @property
    def is_forward(self) -> bool:
        return self.strand == 1


@dataclass
class FastaRecord:
    """Represents a FASTA sequence record"""
    header: str
    sequence: str

    @property
    def id(self) -> str:
        """Extract ID (first word) from header"""
        return self.header.split()[0]

    @property
    def description(self) -> str:
        """Extract description (everything after ID) from header"""
        parts = self.header.split(maxsplit=1)
        return parts[1] if len(parts) > 1 else ""

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass
class GffFeature:
    """Represents a GFF/GFF3 feature"""
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    score: Optional[float]
    strand: str  # '+', '-', or '.'
    phase: Optional[int]
    attributes: Dict[str, str]

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class BedRecord:
    """Represents a BED format record"""
    chrom: str
    start: int  # 0-based
    end: int    # exclusive
    name: Optional[str] = None
    score: Optional[int] = None
    strand: Optional[str] = None

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class DiamondBlastHit:
    """Represents a Diamond BLASTP hit with genomic coordinates"""
    file_basename: str      # Source file identifier
    contig_id: str          # Contig/chromosome identifier
    protein_id: str         # Protein identifier
    start: int              # Genomic start position
    end: int                # Genomic end position
    strand: int             # 1 for forward, -1 for reverse
    query_id: str           # Query sequence ID
    subject_id: str         # Subject sequence ID
    pident: float           # Percentage identity
    alignment_length: int   # Alignment length
    mismatches: int         # Number of mismatches
    gap_opens: int          # Number of gap openings
    query_start: int        # Query alignment start
    query_end: int          # Query alignment end
    subject_start: int      # Subject alignment start
    subject_end: int        # Subject alignment end
    evalue: float           # E-value
    bitscore: float         # Bit score

    @property
    def length(self) -> int:
        return abs(self.end - self.start) + 1

    @property
    def is_forward(self) -> bool:
        return self.strand == 1


@dataclass
class HmmHit:
    """Represents an HMM search hit with genomic coordinates"""
    file_basename: str      # Source file identifier
    contig_id: str          # Contig/chromosome identifier
    protein_id: str         # Protein identifier (target name)
    start: int              # Genomic start position
    end: int                # Genomic end position
    strand: int             # 1 for forward, -1 for reverse
    query_name: str         # HMM profile name
    target_name: str        # Target sequence name
    evalue: float           # Full sequence E-value
    score: float            # Full sequence score
    bias: float             # Bias correction
    domain_evalue: Optional[float] = None  # Best domain E-value
    domain_score: Optional[float] = None   # Best domain score
    domain_number: Optional[int] = None    # Number of domains

    @property
    def length(self) -> int:
        return abs(self.end - self.start) + 1

    @property
    def is_forward(self) -> bool:
        return self.strand == 1


# ============================================
# Prodigal Parsers
# ============================================

def parse_prodigal_faa(faa_path: str) -> Dict[str, Tuple[int, int, int]]:
    """
    Parse Prodigal .faa (protein FASTA) output file.

    Prodigal header format:
    >contig_1 # start # end # strand # ID=1;partial=00;...

    Args:
        faa_path: Path to Prodigal .faa file

    Returns:
        Dict mapping protein_id -> (start, end, strand)
        strand is 1 for forward, -1 for reverse
    """
    position_map = {}

    with open(faa_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract the protein ID (everything after '>' until first space)
                protein_id = line[1:].split()[0]

                # Extract start, end, strand from the header
                # Format: >ID # start # end # strand # ...
                parts = line.split('#')
                if len(parts) >= 4:
                    try:
                        start = int(parts[1].strip())
                        end = int(parts[2].strip())
                        strand = int(parts[3].strip())
                        position_map[protein_id] = (start, end, strand)
                    except ValueError:
                        print(f"Warning: Could not parse positions for {protein_id}")

    return position_map


def parse_prodigal_faa_full(faa_path: str) -> Dict[str, ProdigalGene]:
    """
    Parse Prodigal .faa file and return full gene information including sequences.

    Args:
        faa_path: Path to Prodigal .faa file

    Returns:
        Dict mapping protein_id -> ProdigalGene object
    """
    genes = {}
    current_gene = None
    current_seq = []

    with open(faa_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Save previous gene
                if current_gene is not None:
                    current_gene.sequence = ''.join(current_seq)
                    genes[current_gene.protein_id] = current_gene

                # Parse new gene header
                protein_id = line[1:].split()[0]
                parts = line.split('#')

                if len(parts) >= 4:
                    try:
                        start = int(parts[1].strip())
                        end = int(parts[2].strip())
                        strand = int(parts[3].strip())
                        current_gene = ProdigalGene(
                            protein_id=protein_id,
                            start=start,
                            end=end,
                            strand=strand
                        )
                        current_seq = []
                    except ValueError:
                        print(f"Warning: Could not parse positions for {protein_id}")
                        current_gene = None
                else:
                    current_gene = None
            elif current_gene is not None:
                current_seq.append(line.strip())

        # Save last gene
        if current_gene is not None:
            current_gene.sequence = ''.join(current_seq)
            genes[current_gene.protein_id] = current_gene

    return genes


def parse_prodigal_gff(gff_path: str) -> List[GffFeature]:
    """
    Parse Prodigal GFF output file.

    Args:
        gff_path: Path to Prodigal .gff file

    Returns:
        List of GffFeature objects
    """
    features = []

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            feature = parse_gff_line(line)
            if feature is not None:
                features.append(feature)

    return features


# ============================================
# FASTA Parsers
# ============================================

def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """
    Parse a FASTA file into a dictionary.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Dict mapping sequence_id -> sequence
    """
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def parse_fasta_records(fasta_path: str) -> List[FastaRecord]:
    """
    Parse a FASTA file into FastaRecord objects.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        List of FastaRecord objects
    """
    records = []
    current_header = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_header is not None:
                    records.append(FastaRecord(
                        header=current_header,
                        sequence=''.join(current_seq)
                    ))
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_header is not None:
            records.append(FastaRecord(
                header=current_header,
                sequence=''.join(current_seq)
            ))

    return records


def iter_fasta(fasta_path: str) -> Iterator[FastaRecord]:
    """
    Iterate over FASTA records without loading entire file into memory.

    Args:
        fasta_path: Path to FASTA file

    Yields:
        FastaRecord objects
    """
    current_header = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_header is not None:
                    yield FastaRecord(
                        header=current_header,
                        sequence=''.join(current_seq)
                    )
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_header is not None:
            yield FastaRecord(
                header=current_header,
                sequence=''.join(current_seq)
            )


# ============================================
# GFF/GFF3 Parsers
# ============================================

def parse_gff_line(line: str) -> Optional[GffFeature]:
    """
    Parse a single GFF line.

    Args:
        line: A single line from a GFF file

    Returns:
        GffFeature object or None if line is invalid
    """
    if line.startswith('#') or not line.strip():
        return None

    fields = line.strip().split('\t')
    if len(fields) < 8:
        return None

    try:
        # Parse score
        score = None if fields[5] == '.' else float(fields[5])

        # Parse phase
        phase = None if fields[7] == '.' else int(fields[7])

        # Parse attributes (column 9)
        attributes = {}
        if len(fields) >= 9 and fields[8] != '.':
            for attr in fields[8].split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attributes[key.strip()] = value.strip()
                elif attr.strip():
                    # Handle GFF2 style attributes
                    parts = attr.strip().split(' ', 1)
                    if len(parts) == 2:
                        attributes[parts[0]] = parts[1].strip('"')

        return GffFeature(
            seqid=fields[0],
            source=fields[1],
            feature_type=fields[2],
            start=int(fields[3]),
            end=int(fields[4]),
            score=score,
            strand=fields[6],
            phase=phase,
            attributes=attributes
        )
    except (ValueError, IndexError) as e:
        print(f"Warning: Could not parse GFF line: {line.strip()}")
        return None


def parse_gff(gff_path: str) -> List[GffFeature]:
    """
    Parse a GFF/GFF3 file.

    Args:
        gff_path: Path to GFF file

    Returns:
        List of GffFeature objects
    """
    features = []

    with open(gff_path, 'r') as f:
        for line in f:
            feature = parse_gff_line(line)
            if feature is not None:
                features.append(feature)

    return features


def iter_gff(gff_path: str) -> Iterator[GffFeature]:
    """
    Iterate over GFF features without loading entire file into memory.

    Args:
        gff_path: Path to GFF file

    Yields:
        GffFeature objects
    """
    with open(gff_path, 'r') as f:
        for line in f:
            feature = parse_gff_line(line)
            if feature is not None:
                yield feature


# ============================================
# BED Parsers
# ============================================

def parse_bed(bed_path: str) -> List[BedRecord]:
    """
    Parse a BED file.

    Args:
        bed_path: Path to BED file

    Returns:
        List of BedRecord objects
    """
    records = []

    with open(bed_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track') or not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue

            try:
                record = BedRecord(
                    chrom=fields[0],
                    start=int(fields[1]),
                    end=int(fields[2]),
                    name=fields[3] if len(fields) > 3 else None,
                    score=int(fields[4]) if len(fields) > 4 and fields[4] != '.' else None,
                    strand=fields[5] if len(fields) > 5 else None
                )
                records.append(record)
            except (ValueError, IndexError):
                print(f"Warning: Could not parse BED line: {line.strip()}")

    return records


def iter_bed(bed_path: str) -> Iterator[BedRecord]:
    """
    Iterate over BED records without loading entire file into memory.

    Args:
        bed_path: Path to BED file

    Yields:
        BedRecord objects
    """
    with open(bed_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track') or not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue

            try:
                yield BedRecord(
                    chrom=fields[0],
                    start=int(fields[1]),
                    end=int(fields[2]),
                    name=fields[3] if len(fields) > 3 else None,
                    score=int(fields[4]) if len(fields) > 4 and fields[4] != '.' else None,
                    strand=fields[5] if len(fields) > 5 else None
                )
            except (ValueError, IndexError):
                print(f"Warning: Could not parse BED line: {line.strip()}")


# ============================================
# Diamond BLASTP Parsers
# ============================================

def parse_diamond_blastp(
    blastp_file: str,
    file_basename: Optional[str] = None
) -> List[Tuple[str, str, float, int, int, int, int, int, int, int, float, float]]:
    """
    Parse Diamond BLASTP output in tabular format (-outfmt 6).

    Standard columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

    Args:
        blastp_file: Path to Diamond BLASTP output file
        file_basename: Optional basename to identify source file

    Returns:
        List of tuples with (query_id, subject_id, pident, length, mismatch,
                            gapopen, qstart, qend, sstart, send, evalue, bitscore)
    """
    if file_basename is None:
        file_basename = Path(blastp_file).stem

    hits = []

    with open(blastp_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            try:
                hit = (
                    fields[0],              # query_id
                    fields[1],              # subject_id
                    float(fields[2]),       # pident
                    int(fields[3]),         # alignment_length
                    int(fields[4]),         # mismatches
                    int(fields[5]),         # gap_opens
                    int(fields[6]),         # query_start
                    int(fields[7]),         # query_end
                    int(fields[8]),         # subject_start
                    int(fields[9]),         # subject_end
                    float(fields[10]),      # evalue
                    float(fields[11])       # bitscore
                )
                hits.append(hit)
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse BLASTP line: {line.strip()}")

    return hits


def _lookup_position(
    protein_id: str,
    position_source: PositionSource
) -> Optional[Tuple[int, int, int]]:
    """
    Look up protein position from either a dict or ProdigalIndex.

    Args:
        protein_id: Protein identifier to look up
        position_source: Either a dict mapping protein_id -> (start, end, strand)
                        or a ProdigalIndex object

    Returns:
        Tuple of (start, end, strand) or None if not found
    """
    # Check if it's a dict
    if isinstance(position_source, dict):
        return position_source.get(protein_id)

    # Otherwise assume it's a ProdigalIndex
    if PRODIGAL_INDEX_AVAILABLE and hasattr(position_source, 'get_position_tuple'):
        return position_source.get_position_tuple(protein_id)

    # Fallback: try dict-like access
    try:
        return position_source[protein_id]
    except (KeyError, TypeError):
        return None


def _check_position_exists(
    protein_id: str,
    position_source: PositionSource
) -> bool:
    """
    Check if a protein ID exists in the position source.

    Args:
        protein_id: Protein identifier to check
        position_source: Either a dict or ProdigalIndex

    Returns:
        True if protein_id exists in position_source
    """
    if isinstance(position_source, dict):
        return protein_id in position_source

    # ProdigalIndex supports 'in' operator
    if PRODIGAL_INDEX_AVAILABLE and hasattr(position_source, 'contains'):
        return position_source.contains(protein_id)

    # Fallback
    try:
        return protein_id in position_source
    except TypeError:
        return False


def parse_diamond_blastp_with_positions(
    blastp_file: str,
    position_source: PositionSource,
    file_basename: Optional[str] = None,
    extract_contig_from_protein_id: bool = True
) -> List[DiamondBlastHit]:
    """
    Parse Diamond BLASTP output and add genomic positions from a position source.

    Args:
        blastp_file: Path to Diamond BLASTP output file
        position_source: Either a Dict mapping protein_id -> (start, end, strand)
                        or a ProdigalIndex object for fast SQLite-based lookups.
                        Dict can be obtained from parse_prodigal_faa().
                        ProdigalIndex is recommended for large databases (>1M proteins).
        file_basename: Optional basename to identify source file
        extract_contig_from_protein_id: If True, extract contig_id from protein_id
                                       Assumes format like "contig_name_geneID"

    Returns:
        List of DiamondBlastHit objects with genomic coordinates
    """
    if file_basename is None:
        file_basename = Path(blastp_file).stem

    hits = []

    with open(blastp_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            try:
                query_id = fields[0]
                subject_id = fields[1]

                # Get genomic position from position source
                position = _lookup_position(query_id, position_source)
                if position is None:
                    print(f"Warning: {query_id} not found in position source, skipping")
                    continue

                start, end, strand = position

                # Extract contig ID
                if extract_contig_from_protein_id:
                    # Try to extract contig from protein ID
                    # Common formats: "contig_1", "contig_name_gene_1"
                    parts = query_id.rsplit('_', 1)
                    contig_id = parts[0] if len(parts) > 1 else query_id
                else:
                    contig_id = query_id

                hit = DiamondBlastHit(
                    file_basename=file_basename,
                    contig_id=contig_id,
                    protein_id=query_id,
                    start=start,
                    end=end,
                    strand=strand,
                    query_id=query_id,
                    subject_id=subject_id,
                    pident=float(fields[2]),
                    alignment_length=int(fields[3]),
                    mismatches=int(fields[4]),
                    gap_opens=int(fields[5]),
                    query_start=int(fields[6]),
                    query_end=int(fields[7]),
                    subject_start=int(fields[8]),
                    subject_end=int(fields[9]),
                    evalue=float(fields[10]),
                    bitscore=float(fields[11])
                )
                hits.append(hit)

            except (ValueError, IndexError, KeyError) as e:
                print(f"Warning: Could not parse BLASTP line: {line.strip()}")

    return hits


def iter_diamond_blastp_with_positions(
    blastp_file: str,
    position_source: PositionSource,
    file_basename: Optional[str] = None,
    extract_contig_from_protein_id: bool = True
) -> Iterator[DiamondBlastHit]:
    """
    Iterate over Diamond BLASTP hits without loading entire file into memory.

    Args:
        blastp_file: Path to Diamond BLASTP output file
        position_source: Either a Dict mapping protein_id -> (start, end, strand)
                        or a ProdigalIndex object for fast SQLite-based lookups.
        file_basename: Optional basename to identify source file
        extract_contig_from_protein_id: If True, extract contig_id from protein_id

    Yields:
        DiamondBlastHit objects
    """
    if file_basename is None:
        file_basename = Path(blastp_file).stem

    with open(blastp_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            try:
                query_id = fields[0]
                subject_id = fields[1]

                position = _lookup_position(query_id, position_source)
                if position is None:
                    continue

                start, end, strand = position

                if extract_contig_from_protein_id:
                    parts = query_id.rsplit('_', 1)
                    contig_id = parts[0] if len(parts) > 1 else query_id
                else:
                    contig_id = query_id

                yield DiamondBlastHit(
                    file_basename=file_basename,
                    contig_id=contig_id,
                    protein_id=query_id,
                    start=start,
                    end=end,
                    strand=strand,
                    query_id=query_id,
                    subject_id=subject_id,
                    pident=float(fields[2]),
                    alignment_length=int(fields[3]),
                    mismatches=int(fields[4]),
                    gap_opens=int(fields[5]),
                    query_start=int(fields[6]),
                    query_end=int(fields[7]),
                    subject_start=int(fields[8]),
                    subject_end=int(fields[9]),
                    evalue=float(fields[10]),
                    bitscore=float(fields[11])
                )

            except (ValueError, IndexError, KeyError):
                continue


# ============================================
# HMMER Parsers
# ============================================

def parse_hmm_tblout(hmm_file: str) -> List[Tuple[str, str, str, float, float, float]]:
    """
    Parse HMMER tblout format output (--tblout).

    tblout format columns (space-delimited):
    0. target name
    1. target accession
    2. target length
    3. query name
    4. query accession
    5. query length
    6. full sequence E-value
    7. full sequence score
    8. full sequence bias
    9-18. domain info (best 1 domain, etc.)

    Args:
        hmm_file: Path to HMMER tblout file

    Returns:
        List of tuples with (target_name, query_name, accession,
                            full_evalue, full_score, full_bias)
    """
    hits = []

    with open(hmm_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            # tblout uses spaces, but description can have spaces
            # So we need to be careful with parsing
            fields = line.split()
            if len(fields) < 9:
                continue

            try:
                hit = (
                    fields[0],              # target_name
                    fields[3],              # query_name
                    fields[1],              # target_accession
                    float(fields[6]),       # full_evalue
                    float(fields[7]),       # full_score
                    float(fields[8])        # full_bias
                )
                hits.append(hit)
            except (ValueError, IndexError):
                print(f"Warning: Could not parse HMM line: {line.strip()}")

    return hits


def parse_hmm_tblout_with_positions(
    hmm_file: str,
    position_source: PositionSource,
    file_basename: Optional[str] = None,
    extract_contig_from_protein_id: bool = True
) -> List[HmmHit]:
    """
    Parse HMMER tblout output and add genomic positions from a position source.

    Args:
        hmm_file: Path to HMMER tblout file
        position_source: Either a Dict mapping protein_id -> (start, end, strand)
                        or a ProdigalIndex object for fast SQLite-based lookups.
                        Dict can be obtained from parse_prodigal_faa().
                        ProdigalIndex is recommended for large databases (>1M proteins).
        file_basename: Optional basename to identify source file
        extract_contig_from_protein_id: If True, extract contig_id from protein_id

    Returns:
        List of HmmHit objects with genomic coordinates
    """
    if file_basename is None:
        file_basename = Path(hmm_file).stem

    hits = []

    with open(hmm_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.split()
            if len(fields) < 9:
                continue

            try:
                target_name = fields[0]
                query_name = fields[3]

                # Get genomic position from position source
                position = _lookup_position(target_name, position_source)
                if position is None:
                    print(f"Warning: {target_name} not found in position source, skipping")
                    continue

                start, end, strand = position

                # Extract contig ID
                if extract_contig_from_protein_id:
                    parts = target_name.rsplit('_', 1)
                    contig_id = parts[0] if len(parts) > 1 else target_name
                else:
                    contig_id = target_name

                # Try to parse domain information if available
                domain_evalue = None
                domain_score = None
                domain_number = None
                if len(fields) >= 19:
                    try:
                        domain_number = int(fields[9])  # number of domains
                        domain_evalue = float(fields[12])  # best domain c-Evalue
                        domain_score = float(fields[14])  # best domain score
                    except (ValueError, IndexError):
                        pass

                hit = HmmHit(
                    file_basename=file_basename,
                    contig_id=contig_id,
                    protein_id=target_name,
                    start=start,
                    end=end,
                    strand=strand,
                    query_name=query_name,
                    target_name=target_name,
                    evalue=float(fields[6]),
                    score=float(fields[7]),
                    bias=float(fields[8]),
                    domain_evalue=domain_evalue,
                    domain_score=domain_score,
                    domain_number=domain_number
                )
                hits.append(hit)

            except (ValueError, IndexError, KeyError) as e:
                print(f"Warning: Could not parse HMM line: {line.strip()}")

    return hits


def iter_hmm_tblout_with_positions(
    hmm_file: str,
    position_source: PositionSource,
    file_basename: Optional[str] = None,
    extract_contig_from_protein_id: bool = True
) -> Iterator[HmmHit]:
    """
    Iterate over HMMER tblout hits without loading entire file into memory.

    Args:
        hmm_file: Path to HMMER tblout file
        position_source: Either a Dict mapping protein_id -> (start, end, strand)
                        or a ProdigalIndex object for fast SQLite-based lookups.
        file_basename: Optional basename to identify source file
        extract_contig_from_protein_id: If True, extract contig_id from protein_id

    Yields:
        HmmHit objects
    """
    if file_basename is None:
        file_basename = Path(hmm_file).stem

    with open(hmm_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.split()
            if len(fields) < 9:
                continue

            try:
                target_name = fields[0]
                query_name = fields[3]

                position = _lookup_position(target_name, position_source)
                if position is None:
                    continue

                start, end, strand = position

                if extract_contig_from_protein_id:
                    parts = target_name.rsplit('_', 1)
                    contig_id = parts[0] if len(parts) > 1 else target_name
                else:
                    contig_id = target_name

                domain_evalue = None
                domain_score = None
                domain_number = None
                if len(fields) >= 19:
                    try:
                        domain_number = int(fields[9])
                        domain_evalue = float(fields[12])
                        domain_score = float(fields[14])
                    except (ValueError, IndexError):
                        pass

                yield HmmHit(
                    file_basename=file_basename,
                    contig_id=contig_id,
                    protein_id=target_name,
                    start=start,
                    end=end,
                    strand=strand,
                    query_name=query_name,
                    target_name=target_name,
                    evalue=float(fields[6]),
                    score=float(fields[7]),
                    bias=float(fields[8]),
                    domain_evalue=domain_evalue,
                    domain_score=domain_score,
                    domain_number=domain_number
                )

            except (ValueError, IndexError, KeyError):
                continue


# ============================================
# Sequence Retrieval Functions
# ============================================

def retrieve_sequence_from_fasta(
    fasta_path: str,
    contig_id: str,
    start: int,
    end: int,
    strand: int = 1
) -> Optional[str]:
    """
    Retrieve a sequence from a FASTA file based on coordinates.

    Args:
        fasta_path: Path to FASTA file (.fna, .fasta, .fa)
        contig_id: Contig/chromosome identifier
        start: Start position (1-based, inclusive)
        end: End position (1-based, inclusive)
        strand: 1 for forward, -1 for reverse

    Returns:
        DNA sequence string, or None if not found
    """
    # Parse FASTA file
    sequences = parse_fasta(fasta_path)

    # Try to find the contig
    sequence = None
    for seq_id, seq in sequences.items():
        if seq_id == contig_id or seq_id.startswith(contig_id):
            sequence = seq
            break

    if sequence is None:
        return None

    # Extract region (convert to 0-based)
    region_start = max(0, start - 1)
    region_end = min(len(sequence), end)

    extracted = sequence[region_start:region_end]

    # Reverse complement if on reverse strand
    if strand == -1:
        extracted = reverse_complement(extracted)

    return extracted


def retrieve_sequences_from_hits(
    hits: List,
    fasta_files: List[str],
    verbose: bool = True
) -> List:
    """
    Retrieve sequences for a list of hits from FASTA files.

    Works with DiamondBlastHit and HmmHit objects.

    Args:
        hits: List of hit objects (DiamondBlastHit or HmmHit)
        fasta_files: List of FASTA file paths to search
        verbose: Print progress messages

    Returns:
        List of hits with sequences added (creates new attribute 'genomic_sequence')
    """
    if not hits:
        return hits

    # Build a mapping of contig_id to hits
    contig_to_hits = {}
    for hit in hits:
        if hit.contig_id not in contig_to_hits:
            contig_to_hits[hit.contig_id] = []
        contig_to_hits[hit.contig_id].append(hit)

    if verbose:
        print(f"Searching for {len(hits)} hits across {len(fasta_files)} FASTA files...")

    sequences_found = 0

    # Search through FASTA files
    for fasta_file in fasta_files:
        if verbose:
            print(f"  Processing {Path(fasta_file).name}...")

        try:
            sequences = parse_fasta(fasta_file)

            for seq_id, sequence in sequences.items():
                # Check if this sequence matches any contig we're looking for
                for contig_id, contig_hits in contig_to_hits.items():
                    if seq_id == contig_id or seq_id.startswith(contig_id):
                        # Extract sequences for all hits on this contig
                        for hit in contig_hits:
                            if hasattr(hit, 'genomic_sequence') and hit.genomic_sequence:
                                continue  # Already have sequence

                            # Extract region
                            region_start = max(0, hit.start - 1)
                            region_end = min(len(sequence), hit.end)
                            extracted = sequence[region_start:region_end]

                            # Reverse complement if needed
                            if hit.strand == -1:
                                extracted = reverse_complement(extracted)

                            # Add to hit object
                            hit.genomic_sequence = extracted
                            sequences_found += 1

        except Exception as e:
            if verbose:
                print(f"  Warning: Could not process {fasta_file}: {e}")

    if verbose:
        print(f"Retrieved sequences for {sequences_found}/{len(hits)} hits")

    return hits


def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    Args:
        sequence: DNA sequence string

    Returns:
        Reverse complement sequence
    """
    # Standard bases
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n',
                  # IUPAC ambiguous bases
                  'R': 'Y', 'Y': 'R',  # R=A/G, Y=C/T
                  'W': 'W', 'S': 'S',  # W=A/T (self-complement), S=G/C (self-complement)
                  'M': 'K', 'K': 'M',  # M=A/C, K=G/T
                  'r': 'y', 'y': 'r',
                  'w': 'w', 's': 's',
                  'm': 'k', 'k': 'm',
                  'B': 'V', 'V': 'B',  # B=C/G/T, V=A/C/G
                  'D': 'H', 'H': 'D',  # D=A/G/T, H=A/C/T
                  'b': 'v', 'v': 'b',
                  'd': 'h', 'h': 'd'}

    return ''.join(complement.get(base, base) for base in reversed(sequence))


# ============================================
# DataFrame Conversion Functions
# ============================================

def diamond_hits_to_dataframe(hits: List[DiamondBlastHit]) -> 'pd.DataFrame':
    """
    Convert list of DiamondBlastHit objects to pandas DataFrame.

    Args:
        hits: List of DiamondBlastHit objects

    Returns:
        pandas DataFrame with all hit information

    Raises:
        ImportError: If pandas is not installed
    """
    if not PANDAS_AVAILABLE:
        raise ImportError("pandas is required for DataFrame conversion. Install with: pip install pandas")

    if not hits:
        # Return empty DataFrame with expected columns
        return pd.DataFrame(columns=[
            'file_basename', 'contig_id', 'protein_id', 'start', 'end', 'strand',
            'query_id', 'subject_id', 'pident', 'alignment_length', 'mismatches',
            'gap_opens', 'query_start', 'query_end', 'subject_start', 'subject_end',
            'evalue', 'bitscore', 'length', 'is_forward', 'genomic_sequence'
        ])

    # Convert to list of dicts
    data = []
    for hit in hits:
        hit_dict = asdict(hit)
        # Add computed properties
        hit_dict['length'] = hit.length
        hit_dict['is_forward'] = hit.is_forward
        # Add genomic_sequence if it exists
        if hasattr(hit, 'genomic_sequence'):
            hit_dict['genomic_sequence'] = hit.genomic_sequence
        else:
            hit_dict['genomic_sequence'] = None
        data.append(hit_dict)

    df = pd.DataFrame(data)

    # Set appropriate data types
    df['start'] = df['start'].astype('int64')
    df['end'] = df['end'].astype('int64')
    df['strand'] = df['strand'].astype('int8')
    df['pident'] = df['pident'].astype('float64')
    df['evalue'] = df['evalue'].astype('float64')
    df['bitscore'] = df['bitscore'].astype('float64')

    return df


def hmm_hits_to_dataframe(hits: List[HmmHit]) -> 'pd.DataFrame':
    """
    Convert list of HmmHit objects to pandas DataFrame.

    Args:
        hits: List of HmmHit objects

    Returns:
        pandas DataFrame with all hit information

    Raises:
        ImportError: If pandas is not installed
    """
    if not PANDAS_AVAILABLE:
        raise ImportError("pandas is required for DataFrame conversion. Install with: pip install pandas")

    if not hits:
        # Return empty DataFrame with expected columns
        return pd.DataFrame(columns=[
            'file_basename', 'contig_id', 'protein_id', 'start', 'end', 'strand',
            'query_name', 'target_name', 'evalue', 'score', 'bias',
            'domain_evalue', 'domain_score', 'domain_number',
            'length', 'is_forward', 'genomic_sequence'
        ])

    # Convert to list of dicts
    data = []
    for hit in hits:
        hit_dict = asdict(hit)
        # Add computed properties
        hit_dict['length'] = hit.length
        hit_dict['is_forward'] = hit.is_forward
        # Add genomic_sequence if it exists
        if hasattr(hit, 'genomic_sequence'):
            hit_dict['genomic_sequence'] = hit.genomic_sequence
        else:
            hit_dict['genomic_sequence'] = None
        data.append(hit_dict)

    df = pd.DataFrame(data)

    # Set appropriate data types
    df['start'] = df['start'].astype('int64')
    df['end'] = df['end'].astype('int64')
    df['strand'] = df['strand'].astype('int8')
    df['evalue'] = df['evalue'].astype('float64')
    df['score'] = df['score'].astype('float64')
    df['bias'] = df['bias'].astype('float64')

    return df


def save_hits_to_csv(hits: List, output_file: str, include_sequence: bool = True) -> None:
    """
    Save hits to CSV file.

    Works with both DiamondBlastHit and HmmHit objects.

    Args:
        hits: List of hit objects
        output_file: Path to output CSV file
        include_sequence: Whether to include genomic_sequence column
    """
    if not PANDAS_AVAILABLE:
        raise ImportError("pandas is required for CSV export. Install with: pip install pandas")

    if not hits:
        print("Warning: No hits to save")
        return

    # Determine hit type and convert to DataFrame
    if isinstance(hits[0], DiamondBlastHit):
        df = diamond_hits_to_dataframe(hits)
    elif isinstance(hits[0], HmmHit):
        df = hmm_hits_to_dataframe(hits)
    else:
        raise ValueError(f"Unsupported hit type: {type(hits[0])}")

    # Optionally exclude sequence column (can be very large)
    if not include_sequence and 'genomic_sequence' in df.columns:
        df = df.drop(columns=['genomic_sequence'])

    df.to_csv(output_file, index=False)
    print(f"Saved {len(df)} hits to {output_file}")


def save_hits_to_excel(hits: List, output_file: str, include_sequence: bool = False) -> None:
    """
    Save hits to Excel file.

    Works with both DiamondBlastHit and HmmHit objects.

    Args:
        hits: List of hit objects
        output_file: Path to output Excel file (.xlsx)
        include_sequence: Whether to include genomic_sequence column (default: False for Excel)

    Note:
        Requires openpyxl: pip install openpyxl
    """
    if not PANDAS_AVAILABLE:
        raise ImportError("pandas is required for Excel export. Install with: pip install pandas openpyxl")

    if not hits:
        print("Warning: No hits to save")
        return

    # Determine hit type and convert to DataFrame
    if isinstance(hits[0], DiamondBlastHit):
        df = diamond_hits_to_dataframe(hits)
    elif isinstance(hits[0], HmmHit):
        df = hmm_hits_to_dataframe(hits)
    else:
        raise ValueError(f"Unsupported hit type: {type(hits[0])}")

    # Optionally exclude sequence column
    if not include_sequence and 'genomic_sequence' in df.columns:
        df = df.drop(columns=['genomic_sequence'])

    df.to_excel(output_file, index=False, engine='openpyxl')
    print(f"Saved {len(df)} hits to {output_file}")


# ============================================
# Utility Functions
# ============================================

def detect_file_format(file_path: str) -> str:
    """
    Detect file format based on extension and content.

    Args:
        file_path: Path to file

    Returns:
        Format string: 'fasta', 'gff', 'bed', or 'unknown'
    """
    path = Path(file_path)
    ext = path.suffix.lower()

    # Check by extension first
    if ext in ['.fa', '.fasta', '.fna', '.faa', '.ffn']:
        return 'fasta'
    elif ext in ['.gff', '.gff3', '.gtf']:
        return 'gff'
    elif ext in ['.bed']:
        return 'bed'

    # Check by content
    with open(file_path, 'r') as f:
        first_line = f.readline()
        if first_line.startswith('>'):
            return 'fasta'
        elif first_line.startswith('##gff') or '\t' in first_line:
            return 'gff'

    return 'unknown'


# ============================================
# Example Usage
# ============================================

if __name__ == "__main__":
    print("Parsers module - example usage")
    print("=" * 50)

    # Example: Parse Prodigal output
    example_prodigal = """
    from utils.parsers import parse_prodigal_faa

    positions = parse_prodigal_faa("proteins.faa")
    for protein_id, (start, end, strand) in positions.items():
        print(f"{protein_id}: {start}-{end} ({'+' if strand == 1 else '-'})")
    """
    print("Prodigal parsing:")
    print(example_prodigal)

    # Example: Parse FASTA
    example_fasta = """
    from utils.parsers import parse_fasta, iter_fasta

    # Load all sequences
    sequences = parse_fasta("sequences.fasta")

    # Or iterate for large files
    for record in iter_fasta("large_file.fasta"):
        print(f"{record.id}: {record.length} bp")
    """
    print("FASTA parsing:")
    print(example_fasta)
