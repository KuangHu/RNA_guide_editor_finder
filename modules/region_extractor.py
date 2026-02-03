#!/usr/bin/env python3
"""
Region Extractor Module

Extracts coding regions and flanking non-coding regions from genomic sequences.
Handles both forward and reverse strand orientations.

Typical usage:
    extractor = RegionExtractor(genome_fna_path)
    result = extractor.extract_protein_regions(
        contig_id="AAXX01000001.1",
        start=1000,
        end=2000,
        strand=1,
        upstream_length=350,
        downstream_length=250
    )
"""

from pathlib import Path
from typing import Dict, Tuple, Optional
from Bio import SeqIO
import sys

# Import utility functions from parsers module
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.parsers import parse_prodigal_faa, reverse_complement


class RegionExtractor:
    """
    Extracts coding and non-coding regions from genomic sequences.

    Attributes:
        genome_sequences (Dict[str, str]): Loaded genomic sequences {contig_id: sequence}
    """

    def __init__(self, genome_fna_path: Optional[str] = None):
        """
        Initialize the RegionExtractor.

        Args:
            genome_fna_path: Path to genomic .fna file. If provided, sequences are loaded.
        """
        self.genome_sequences = {}
        if genome_fna_path:
            self.load_genome(genome_fna_path)

    def load_genome(self, fna_path: str) -> None:
        """
        Load genomic sequences from .fna file.

        Args:
            fna_path: Path to the genomic FASTA file

        Raises:
            FileNotFoundError: If the file doesn't exist
        """
        if not Path(fna_path).exists():
            raise FileNotFoundError(f"Genome file not found: {fna_path}")

        self.genome_sequences = {}
        for record in SeqIO.parse(fna_path, "fasta"):
            self.genome_sequences[record.id] = str(record.seq)

    def extract_protein_regions(
        self,
        contig_id: str,
        start: int,
        end: int,
        strand: int,
        upstream_length: int = 350,
        downstream_length: int = 250
    ) -> Dict:
        """
        Extract coding region and flanking sequences for a protein.

        Coordinates are 1-based (Prodigal/GFF format).
        All sequences are oriented according to the gene strand.

        Args:
            contig_id: Contig/chromosome ID
            start: Gene/protein coding region start position (1-based, inclusive)
                   For transposons: this is the transposase gene start
            end: Gene/protein coding region end position (1-based, inclusive)
                 For transposons: this is the transposase gene end
            strand: Gene strand (1 for forward, -1 for reverse)
            upstream_length: Length of upstream region to extract (default: 350)
            downstream_length: Length of downstream region to extract (default: 250)

        Returns:
            Dictionary containing:
                - coding_sequence: Gene sequence (oriented by strand)
                - upstream_sequence: Upstream flanking sequence
                - downstream_sequence: Downstream flanking sequence
                - coding_coords: Dict with start, end, length
                - upstream_coords: Dict with start, end, length
                - downstream_coords: Dict with start, end, length
                - contig_id: Contig identifier
                - strand: Gene strand

        Raises:
            KeyError: If contig_id not found in loaded sequences
            ValueError: If coordinates are invalid
        """
        if contig_id not in self.genome_sequences:
            raise KeyError(f"Contig '{contig_id}' not found in loaded genome sequences")

        genome_seq = self.genome_sequences[contig_id]
        genome_len = len(genome_seq)

        # Validate coordinates
        if start < 1 or end > genome_len or start > end:
            raise ValueError(
                f"Invalid coordinates: start={start}, end={end}, genome_length={genome_len}"
            )

        # Convert to 0-based indexing for Python
        start_0 = start - 1
        end_0 = end

        if strand == 1:  # Forward strand
            regions = self._extract_forward_strand(
                genome_seq, start_0, end_0, upstream_length, downstream_length, genome_len
            )
        elif strand == -1:  # Reverse strand
            regions = self._extract_reverse_strand(
                genome_seq, start_0, end_0, upstream_length, downstream_length, genome_len
            )
        else:
            raise ValueError(f"Invalid strand value: {strand}. Must be 1 or -1")

        # Add metadata
        regions['contig_id'] = contig_id
        regions['strand'] = strand

        return regions

    def _extract_forward_strand(
        self,
        genome_seq: str,
        start_0: int,
        end_0: int,
        upstream_len: int,
        downstream_len: int,
        genome_len: int
    ) -> Dict:
        """Extract regions for forward strand gene."""
        # Coding region
        coding_seq = genome_seq[start_0:end_0]

        # Upstream (before start)
        upstream_start = max(0, start_0 - upstream_len)
        upstream_seq = genome_seq[upstream_start:start_0]

        # Downstream (after end)
        downstream_end = min(genome_len, end_0 + downstream_len)
        downstream_seq = genome_seq[end_0:downstream_end]

        return {
            'coding_sequence': coding_seq,
            'upstream_sequence': upstream_seq,
            'downstream_sequence': downstream_seq,
            'coding_coords': {
                'start': start_0 + 1,  # Convert back to 1-based
                'end': end_0,
                'length': end_0 - start_0
            },
            'upstream_coords': {
                'start': upstream_start + 1,  # 1-based
                'end': start_0,  # This is start - 1 in 1-based coords
                'length': len(upstream_seq)
            },
            'downstream_coords': {
                'start': end_0 + 1,  # 1-based
                'end': downstream_end,
                'length': len(downstream_seq)
            }
        }

    def _extract_reverse_strand(
        self,
        genome_seq: str,
        start_0: int,
        end_0: int,
        upstream_len: int,
        downstream_len: int,
        genome_len: int
    ) -> Dict:
        """
        Extract regions for reverse strand gene.

        For reverse strand:
        - "upstream" is genomically downstream (after end)
        - "downstream" is genomically upstream (before start)
        - All sequences are reverse complemented
        """
        # Coding region (reverse complement)
        coding_seq = genome_seq[start_0:end_0]
        coding_seq_rc = reverse_complement(coding_seq)

        # Upstream (in gene orientation, genomically downstream)
        upstream_end = min(genome_len, end_0 + upstream_len)
        upstream_seq = genome_seq[end_0:upstream_end]
        upstream_seq_rc = reverse_complement(upstream_seq)

        # Downstream (in gene orientation, genomically upstream)
        downstream_start = max(0, start_0 - downstream_len)
        downstream_seq = genome_seq[downstream_start:start_0]
        downstream_seq_rc = reverse_complement(downstream_seq)

        return {
            'coding_sequence': coding_seq_rc,
            'upstream_sequence': upstream_seq_rc,
            'downstream_sequence': downstream_seq_rc,
            'coding_coords': {
                'start': start_0 + 1,  # 1-based
                'end': end_0,
                'length': end_0 - start_0
            },
            'upstream_coords': {
                'start': end_0 + 1,  # 1-based
                'end': upstream_end,
                'length': len(upstream_seq)
            },
            'downstream_coords': {
                'start': downstream_start + 1,  # 1-based
                'end': start_0,  # This is start - 1 in 1-based
                'length': len(downstream_seq)
            }
        }


# Note: reverse_complement and parse_prodigal_faa are imported from utils.parsers
# to avoid code duplication


def extract_contig_from_protein_id(protein_id: str) -> str:
    """
    Extract contig ID from protein ID.

    Args:
        protein_id: Protein identifier

    Returns:
        Contig identifier

    Example:
        'AAXX01000001.1_281' -> 'AAXX01000001.1'
        'NC_000001.11_123' -> 'NC_000001.11'
    """
    return '_'.join(protein_id.split('_')[:-1])


def create_transposon_dict(
    genome_fna_path: str,
    contig_id: str,
    protein_id: str,
    start: int,
    end: int,
    strand: int,
    upstream_length: int = 350,
    downstream_length: int = 250,
    additional_metadata: Optional[Dict] = None
) -> Dict:
    """
    Create a complete transposon dictionary entry.

    This is a convenience function that creates a RegionExtractor,
    extracts regions, and formats everything into a dictionary.

    Args:
        genome_fna_path: Path to genomic .fna file
        contig_id: Contig identifier
        protein_id: Protein identifier (e.g., transposase)
        start: Coding region start position (1-based) - the transposase gene start
        end: Coding region end position (1-based) - the transposase gene end
        strand: Gene strand (1 or -1)
        upstream_length: Length of upstream region (default: 350)
        downstream_length: Length of downstream region (default: 250)
        additional_metadata: Optional dict of additional fields to include

    Returns:
        Dictionary with protein information and extracted sequences
    """
    extractor = RegionExtractor(genome_fna_path)

    regions = extractor.extract_protein_regions(
        contig_id=contig_id,
        start=start,
        end=end,
        strand=strand,
        upstream_length=upstream_length,
        downstream_length=downstream_length
    )

    # Build complete dictionary
    transposon_dict = {
        'protein_id': protein_id,
        'contig_id': contig_id,
        'strand': strand,
        'coordinates': {
            'coding': regions['coding_coords'],
            'upstream': regions['upstream_coords'],
            'downstream': regions['downstream_coords']
        },
        'sequences': {
            'coding': regions['coding_sequence'],
            'upstream': regions['upstream_sequence'],
            'downstream': regions['downstream_sequence']
        }
    }

    # Add any additional metadata
    if additional_metadata:
        transposon_dict.update(additional_metadata)

    return transposon_dict
