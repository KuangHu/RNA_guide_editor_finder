"""
Sequence Search Module for RNA Guide Editor Finder

This module handles:
1. Query preparation (DNA sequence -> FASTA)
2. MMseqs2 search against GTDB database
3. Hit sequence retrieval from .fna files
4. Deduplication of results
"""

import subprocess
import tempfile
import os
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass


@dataclass
class SearchHit:
    """Represents a single search hit from MMseqs2."""
    query_id: str
    target_id: str
    sequence_identity: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    evalue: float
    bit_score: float
    sequence: Optional[str] = None


class SequenceSearcher:
    """
    Handles sequence search using MMseqs2 against genomic databases.
    """

    def __init__(
        self,
        database_path: str,
        fna_folder: str,
        min_seq_identity: float = 0.8,
        sensitivity: float = 7.5,
        threads: int = 4
    ):
        """
        Initialize the sequence searcher.

        Args:
            database_path: Path to the MMseqs2 database (e.g., GTDB)
            fna_folder: Path to folder containing .fna files with full sequences
            min_seq_identity: Minimum sequence identity threshold (default: 0.8)
            sensitivity: MMseqs2 sensitivity parameter (default: 7.5)
            threads: Number of threads for MMseqs2 (default: 4)
        """
        self.database_path = Path(database_path)
        self.fna_folder = Path(fna_folder)
        self.min_seq_identity = min_seq_identity
        self.sensitivity = sensitivity
        self.threads = threads

        # Validate paths
        if not self.database_path.exists():
            raise FileNotFoundError(f"Database not found: {database_path}")
        if not self.fna_folder.exists():
            raise FileNotFoundError(f"FNA folder not found: {fna_folder}")

    def create_query_fasta(self, sequence: str, query_id: str = "query") -> str:
        """
        Create a temporary FASTA file from a DNA sequence.

        Args:
            sequence: DNA sequence string (50-200bp)
            query_id: Identifier for the query sequence

        Returns:
            Path to the temporary FASTA file
        """
        # Validate sequence
        sequence = sequence.upper().strip()
        valid_bases = set("ATCGN")
        if not all(base in valid_bases for base in sequence):
            raise ValueError("Invalid DNA sequence: contains non-ATCGN characters")

        if not (50 <= len(sequence) <= 200):
            print(f"Warning: Sequence length ({len(sequence)}bp) is outside recommended range (50-200bp)")

        # Create temp file
        temp_fasta = tempfile.NamedTemporaryFile(
            mode='w', suffix='.fasta', delete=False
        )
        temp_fasta.write(f">{query_id}\n{sequence}\n")
        temp_fasta.close()

        return temp_fasta.name

    def run_mmseqs_search(
        self,
        query_fasta: str,
        output_file: Optional[str] = None,
        tmp_dir: Optional[str] = None
    ) -> str:
        """
        Run MMseqs2 easy-search against the database.

        Args:
            query_fasta: Path to query FASTA file
            output_file: Path for output file (default: temp file)
            tmp_dir: Temporary directory for MMseqs2 (default: auto-created)

        Returns:
            Path to the results file (m8 format)
        """
        # Set up output and temp paths
        if output_file is None:
            output_file = tempfile.NamedTemporaryFile(
                suffix='.m8', delete=False
            ).name

        if tmp_dir is None:
            tmp_dir = tempfile.mkdtemp(prefix='mmseqs_tmp_')

        # Build MMseqs2 command
        cmd = [
            'mmseqs', 'easy-search',
            query_fasta,
            str(self.database_path),
            output_file,
            tmp_dir,
            '-s', str(self.sensitivity),
            '--min-seq-id', str(self.min_seq_identity),
            '--threads', str(self.threads)
        ]

        print(f"Running MMseqs2 search: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            print("MMseqs2 search completed successfully")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"MMseqs2 search failed: {e.stderr}")

        return output_file

    def parse_mmseqs_results(self, results_file: str) -> List[SearchHit]:
        """
        Parse MMseqs2 m8 format results.

        Args:
            results_file: Path to the m8 results file

        Returns:
            List of SearchHit objects
        """
        hits = []

        with open(results_file, 'r') as f:
            for line in f:
                if line.strip():
                    fields = line.strip().split('\t')
                    if len(fields) >= 12:
                        hit = SearchHit(
                            query_id=fields[0],
                            target_id=fields[1],
                            sequence_identity=float(fields[2]),
                            alignment_length=int(fields[3]),
                            mismatches=int(fields[4]),
                            gap_opens=int(fields[5]),
                            query_start=int(fields[6]),
                            query_end=int(fields[7]),
                            target_start=int(fields[8]),
                            target_end=int(fields[9]),
                            evalue=float(fields[10]),
                            bit_score=float(fields[11])
                        )
                        hits.append(hit)

        print(f"Parsed {len(hits)} hits from results")
        return hits

    def retrieve_sequences_from_fna(
        self,
        hits: List[SearchHit],
        context_bp: int = 0
    ) -> List[SearchHit]:
        """
        Retrieve full sequences for hits from .fna files.

        Args:
            hits: List of SearchHit objects
            context_bp: Additional context bases to include on each side

        Returns:
            List of SearchHit objects with sequences populated
        """
        # Build index of target IDs to hits
        target_to_hits: Dict[str, List[SearchHit]] = {}
        for hit in hits:
            # Extract genome/contig identifier from target_id
            # Format may vary, common patterns: "genome_id_contig" or "accession"
            target_key = hit.target_id.split('_')[0] if '_' in hit.target_id else hit.target_id
            if target_key not in target_to_hits:
                target_to_hits[target_key] = []
            target_to_hits[target_key].append(hit)

        # Search through .fna files
        fna_files = list(self.fna_folder.glob("**/*.fna")) + \
                    list(self.fna_folder.glob("**/*.fasta")) + \
                    list(self.fna_folder.glob("**/*.fa"))

        print(f"Searching through {len(fna_files)} sequence files...")

        sequences_found = 0
        for fna_file in fna_files:
            current_header = None
            current_seq = []

            with open(fna_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        # Process previous sequence
                        if current_header and current_seq:
                            self._assign_sequence_to_hits(
                                current_header, ''.join(current_seq),
                                hits, context_bp
                            )

                        current_header = line[1:].strip().split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line.strip())

                # Process last sequence in file
                if current_header and current_seq:
                    count = self._assign_sequence_to_hits(
                        current_header, ''.join(current_seq),
                        hits, context_bp
                    )
                    sequences_found += count

        # Count how many hits have sequences
        with_seq = sum(1 for h in hits if h.sequence is not None)
        print(f"Retrieved sequences for {with_seq}/{len(hits)} hits")

        return hits

    def _assign_sequence_to_hits(
        self,
        header: str,
        full_sequence: str,
        hits: List[SearchHit],
        context_bp: int
    ) -> int:
        """
        Assign sequence to matching hits.

        Returns:
            Number of hits that received sequences
        """
        count = 0
        for hit in hits:
            if hit.target_id == header or header.startswith(hit.target_id):
                # Extract the target region with optional context
                start = max(0, hit.target_start - 1 - context_bp)
                end = min(len(full_sequence), hit.target_end + context_bp)
                hit.sequence = full_sequence[start:end]
                count += 1
        return count

    def deduplicate_hits(
        self,
        hits: List[SearchHit],
        identity_threshold: float = 0.95
    ) -> List[SearchHit]:
        """
        Deduplicate hits based on sequence similarity.

        Uses a simple greedy clustering approach:
        1. Sort by bit score (descending)
        2. Keep hit if not too similar to any kept hit

        Args:
            hits: List of SearchHit objects with sequences
            identity_threshold: Maximum identity to consider as duplicate

        Returns:
            Deduplicated list of SearchHit objects
        """
        if not hits:
            return []

        # Filter to only hits with sequences
        hits_with_seq = [h for h in hits if h.sequence]

        if not hits_with_seq:
            print("Warning: No hits have sequences, returning all hits")
            return hits

        # Sort by bit score
        sorted_hits = sorted(hits_with_seq, key=lambda h: h.bit_score, reverse=True)

        deduplicated = []

        for hit in sorted_hits:
            is_duplicate = False

            for kept_hit in deduplicated:
                identity = self._calculate_sequence_identity(
                    hit.sequence, kept_hit.sequence
                )
                if identity >= identity_threshold:
                    is_duplicate = True
                    break

            if not is_duplicate:
                deduplicated.append(hit)

        print(f"Deduplicated {len(hits_with_seq)} -> {len(deduplicated)} hits")
        return deduplicated

    def _calculate_sequence_identity(self, seq1: str, seq2: str) -> float:
        """
        Calculate sequence identity between two sequences.

        Simple implementation using global alignment concept.
        For more accurate results, consider using a proper alignment library.
        """
        if not seq1 or not seq2:
            return 0.0

        # Simple identity calculation for similar-length sequences
        min_len = min(len(seq1), len(seq2))
        max_len = max(len(seq1), len(seq2))

        if min_len == 0:
            return 0.0

        # Count matches at each position
        matches = sum(
            1 for i in range(min_len)
            if seq1[i].upper() == seq2[i].upper()
        )

        return matches / max_len

    def search(
        self,
        sequence: str,
        query_id: str = "query",
        context_bp: int = 0,
        dedup_threshold: float = 0.95,
        cleanup: bool = True
    ) -> List[SearchHit]:
        """
        Complete search pipeline: query -> search -> retrieve -> deduplicate.

        Args:
            sequence: DNA query sequence (50-200bp)
            query_id: Identifier for the query
            context_bp: Additional context to include around hits
            dedup_threshold: Identity threshold for deduplication
            cleanup: Whether to remove temporary files

        Returns:
            Deduplicated list of SearchHit objects with sequences
        """
        temp_files = []

        try:
            # Step 1: Create query FASTA
            print("Step 1: Creating query FASTA...")
            query_fasta = self.create_query_fasta(sequence, query_id)
            temp_files.append(query_fasta)

            # Step 2: Run MMseqs2 search
            print("Step 2: Running MMseqs2 search...")
            results_file = self.run_mmseqs_search(query_fasta)
            temp_files.append(results_file)

            # Step 3: Parse results
            print("Step 3: Parsing results...")
            hits = self.parse_mmseqs_results(results_file)

            if not hits:
                print("No hits found")
                return []

            # Step 4: Retrieve sequences
            print("Step 4: Retrieving sequences from FNA files...")
            hits = self.retrieve_sequences_from_fna(hits, context_bp)

            # Step 5: Deduplicate
            print("Step 5: Deduplicating results...")
            deduplicated = self.deduplicate_hits(hits, dedup_threshold)

            return deduplicated

        finally:
            # Cleanup temporary files
            if cleanup:
                for temp_file in temp_files:
                    if os.path.exists(temp_file):
                        os.remove(temp_file)


def write_results_fasta(
    hits: List[SearchHit],
    output_file: str,
    include_metadata: bool = True
) -> None:
    """
    Write search results to a FASTA file.

    Args:
        hits: List of SearchHit objects
        output_file: Path to output FASTA file
        include_metadata: Whether to include hit metadata in headers
    """
    with open(output_file, 'w') as f:
        for i, hit in enumerate(hits):
            if hit.sequence:
                if include_metadata:
                    header = (
                        f">{hit.target_id} "
                        f"identity={hit.sequence_identity:.2f} "
                        f"evalue={hit.evalue:.2e} "
                        f"score={hit.bit_score:.1f}"
                    )
                else:
                    header = f">{hit.target_id}"

                f.write(f"{header}\n{hit.sequence}\n")

    print(f"Wrote {len([h for h in hits if h.sequence])} sequences to {output_file}")


# Default paths for GTDB database
DEFAULT_DATABASE_PATH = "/groups/rubin/projects/kuang/db/gtdb_mmseq_db"
DEFAULT_FNA_FOLDER = "/groups/rubin/projects/kuang/db/GTDB/"


# Example usage
if __name__ == "__main__":
    # Example query sequence (50-200bp)
    test_sequence = "tatccctccagtgcagagaaaatcggccagttttctctgcctgcagtccgcatgccgtatcgggccttgggttctaacctgttgcgtagatttatgcagcggactgcctttctcccaaagtgataaaccggacagtatcatggaccggttttcccggtaatccgtatttgcaaggttggtttcact"

    # Initialize searcher with default paths
    searcher = SequenceSearcher(
        database_path=DEFAULT_DATABASE_PATH,
        fna_folder=DEFAULT_FNA_FOLDER,
        min_seq_identity=0.8,
        sensitivity=7.5
    )

    # Run search
    results = searcher.search(
        sequence=test_sequence,
        dedup_threshold=0.95
    )

    # Write output
    write_results_fasta(results, "output.fasta")
