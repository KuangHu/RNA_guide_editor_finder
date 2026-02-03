#!/usr/bin/env python3
"""
Prodigal Index Module

Fast SQLite-based index for looking up protein positions from large Prodigal FAA files.
Designed for files with hundreds of millions of proteins (e.g., GTDB combined FAA).

Usage:
    # Build index (one-time, ~30-60 min for 300M proteins)
    build_prodigal_index("combine.faa", "combine.faa.sqlite")

    # Fast lookups
    idx = ProdigalIndex("combine.faa.sqlite")
    pos = idx.get_position("BA000021.3_1")  # Returns (start, end, strand)

    # Batch lookups for DIAMOND/HMMER hits
    positions = idx.get_positions(["prot1", "prot2", "prot3"])

Author: Kuang Hu
Date: 2026-01-30
"""

import sqlite3
import os
import sys
import time
from typing import Dict, List, Tuple, Optional, Iterator
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ProteinPosition:
    """Position information for a protein."""
    protein_id: str
    start: int
    end: int
    strand: int
    contig_id: Optional[str] = None

    @property
    def length(self) -> int:
        """Genomic length of the protein-coding region."""
        return abs(self.end - self.start) + 1

    def as_tuple(self) -> Tuple[int, int, int]:
        """Return (start, end, strand) tuple."""
        return (self.start, self.end, self.strand)


class ProdigalIndex:
    """
    Fast SQLite-based index for Prodigal protein positions.

    Provides O(1) lookups by protein_id with minimal memory usage.
    Supports batch queries for processing DIAMOND/HMMER results.

    Attributes:
        db_path: Path to the SQLite database file

    Example:
        >>> idx = ProdigalIndex("proteins.sqlite")
        >>> pos = idx.get_position("contig_1")
        >>> print(f"Position: {pos.start}-{pos.end}, strand: {pos.strand}")

        >>> # Batch lookup
        >>> hits = ["prot1", "prot2", "prot3"]
        >>> positions = idx.get_positions(hits)
        >>> for pid, pos in positions.items():
        ...     print(f"{pid}: {pos.start}-{pos.end}")
    """

    def __init__(self, db_path: str):
        """
        Open an existing Prodigal index database.

        Args:
            db_path: Path to the SQLite database file

        Raises:
            FileNotFoundError: If database file doesn't exist
        """
        # Initialize attributes first (for safe __del__)
        self._conn = None
        self._cursor = None
        self._count = None

        if not os.path.exists(db_path):
            raise FileNotFoundError(f"Index database not found: {db_path}")

        self.db_path = db_path
        self._connect()

    def _connect(self):
        """Establish database connection with optimized settings."""
        self._conn = sqlite3.connect(self.db_path)
        self._conn.row_factory = sqlite3.Row
        self._cursor = self._conn.cursor()

        # Optimize for read-only access
        self._cursor.execute("PRAGMA query_only = ON")
        self._cursor.execute("PRAGMA cache_size = -64000")  # 64MB cache
        self._cursor.execute("PRAGMA mmap_size = 1073741824")  # 1GB mmap

    def close(self):
        """Close the database connection."""
        if self._conn:
            self._conn.close()
            self._conn = None
            self._cursor = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __len__(self) -> int:
        """Return total number of proteins in the index."""
        if self._count is None:
            self._cursor.execute("SELECT COUNT(*) FROM proteins")
            self._count = self._cursor.fetchone()[0]
        return self._count

    def get_position(self, protein_id: str) -> Optional[ProteinPosition]:
        """
        Get position for a single protein.

        Args:
            protein_id: Protein identifier (e.g., "BA000021.3_1")

        Returns:
            ProteinPosition object or None if not found
        """
        self._cursor.execute(
            "SELECT protein_id, start, end, strand, contig_id FROM proteins WHERE protein_id = ?",
            (protein_id,)
        )
        row = self._cursor.fetchone()

        if row:
            return ProteinPosition(
                protein_id=row['protein_id'],
                start=row['start'],
                end=row['end'],
                strand=row['strand'],
                contig_id=row['contig_id']
            )
        return None

    def get_position_tuple(self, protein_id: str) -> Optional[Tuple[int, int, int]]:
        """
        Get position as tuple (start, end, strand) for a single protein.

        Slightly faster than get_position() when you just need the coordinates.

        Args:
            protein_id: Protein identifier

        Returns:
            Tuple of (start, end, strand) or None if not found
        """
        self._cursor.execute(
            "SELECT start, end, strand FROM proteins WHERE protein_id = ?",
            (protein_id,)
        )
        row = self._cursor.fetchone()

        if row:
            return (row['start'], row['end'], row['strand'])
        return None

    def get_positions(self, protein_ids: List[str]) -> Dict[str, ProteinPosition]:
        """
        Batch lookup for multiple proteins.

        More efficient than calling get_position() in a loop.

        Args:
            protein_ids: List of protein identifiers

        Returns:
            Dict mapping protein_id -> ProteinPosition (only found proteins included)
        """
        if not protein_ids:
            return {}

        # Use chunked queries for very large batches
        chunk_size = 900  # SQLite has a limit on variables
        results = {}

        for i in range(0, len(protein_ids), chunk_size):
            chunk = protein_ids[i:i + chunk_size]
            placeholders = ','.join(['?' for _ in chunk])

            self._cursor.execute(
                f"SELECT protein_id, start, end, strand, contig_id FROM proteins WHERE protein_id IN ({placeholders})",
                chunk
            )

            for row in self._cursor.fetchall():
                results[row['protein_id']] = ProteinPosition(
                    protein_id=row['protein_id'],
                    start=row['start'],
                    end=row['end'],
                    strand=row['strand'],
                    contig_id=row['contig_id']
                )

        return results

    def get_positions_tuples(self, protein_ids: List[str]) -> Dict[str, Tuple[int, int, int]]:
        """
        Batch lookup returning tuples instead of ProteinPosition objects.

        Args:
            protein_ids: List of protein identifiers

        Returns:
            Dict mapping protein_id -> (start, end, strand)
        """
        if not protein_ids:
            return {}

        chunk_size = 900
        results = {}

        for i in range(0, len(protein_ids), chunk_size):
            chunk = protein_ids[i:i + chunk_size]
            placeholders = ','.join(['?' for _ in chunk])

            self._cursor.execute(
                f"SELECT protein_id, start, end, strand FROM proteins WHERE protein_id IN ({placeholders})",
                chunk
            )

            for row in self._cursor.fetchall():
                results[row['protein_id']] = (row['start'], row['end'], row['strand'])

        return results

    def get_proteins_by_contig(self, contig_id: str) -> List[ProteinPosition]:
        """
        Get all proteins from a specific contig.

        Args:
            contig_id: Contig identifier

        Returns:
            List of ProteinPosition objects sorted by start position
        """
        self._cursor.execute(
            "SELECT protein_id, start, end, strand, contig_id FROM proteins WHERE contig_id = ? ORDER BY start",
            (contig_id,)
        )

        return [
            ProteinPosition(
                protein_id=row['protein_id'],
                start=row['start'],
                end=row['end'],
                strand=row['strand'],
                contig_id=row['contig_id']
            )
            for row in self._cursor.fetchall()
        ]

    def contains(self, protein_id: str) -> bool:
        """Check if a protein exists in the index."""
        self._cursor.execute(
            "SELECT 1 FROM proteins WHERE protein_id = ? LIMIT 1",
            (protein_id,)
        )
        return self._cursor.fetchone() is not None

    def __contains__(self, protein_id: str) -> bool:
        """Support 'in' operator: if protein_id in idx: ..."""
        return self.contains(protein_id)

    def stats(self) -> Dict:
        """Get index statistics."""
        self._cursor.execute("SELECT COUNT(*) FROM proteins")
        total = self._cursor.fetchone()[0]

        self._cursor.execute("SELECT COUNT(DISTINCT contig_id) FROM proteins")
        contigs = self._cursor.fetchone()[0]

        return {
            'total_proteins': total,
            'total_contigs': contigs,
            'db_size_mb': os.path.getsize(self.db_path) / (1024 * 1024)
        }


def _extract_contig_id(protein_id: str) -> str:
    """
    Extract contig ID from protein ID.

    Prodigal format: contig_genenum (e.g., "BA000021.3_1" -> "BA000021.3")
    """
    parts = protein_id.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return protein_id


def _parse_prodigal_header(line: str) -> Optional[Tuple[str, int, int, int, str]]:
    """
    Parse a Prodigal FAA header line.

    Format: >protein_id # start # end # strand # attributes

    Returns:
        Tuple of (protein_id, start, end, strand, contig_id) or None if invalid
    """
    if not line.startswith('>'):
        return None

    # Extract protein ID
    protein_id = line[1:].split()[0]

    # Extract positions from # delimited fields
    parts = line.split('#')
    if len(parts) >= 4:
        try:
            start = int(parts[1].strip())
            end = int(parts[2].strip())
            strand = int(parts[3].strip())
            contig_id = _extract_contig_id(protein_id)
            return (protein_id, start, end, strand, contig_id)
        except (ValueError, IndexError):
            pass

    return None


def build_prodigal_index(
    faa_path: str,
    db_path: str,
    chunk_size: int = 100000,
    show_progress: bool = True
) -> str:
    """
    Build a SQLite index from a Prodigal FAA file.

    Creates an optimized SQLite database for fast protein position lookups.
    Designed to handle files with hundreds of millions of proteins.

    Args:
        faa_path: Path to Prodigal .faa file
        db_path: Path for output SQLite database
        chunk_size: Number of records to insert per transaction (default: 100000)
        show_progress: Whether to print progress updates

    Returns:
        Path to the created database

    Example:
        >>> build_prodigal_index("combine.faa", "combine.faa.sqlite")
        Building index: combine.faa -> combine.faa.sqlite
        Processed 100,000,000 proteins...
        Processed 200,000,000 proteins...
        Processed 300,000,000 proteins...
        Done! 322,484,155 proteins indexed in 45.2 minutes
        Creating indexes...
        Index complete: combine.faa.sqlite (12.5 GB)
    """
    if not os.path.exists(faa_path):
        raise FileNotFoundError(f"FAA file not found: {faa_path}")

    # Remove existing database
    if os.path.exists(db_path):
        os.remove(db_path)

    if show_progress:
        print(f"Building index: {faa_path} -> {db_path}")

    start_time = time.time()

    # Create database with optimized settings for bulk insert
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Optimize for bulk insert
    cursor.execute("PRAGMA synchronous = OFF")
    cursor.execute("PRAGMA journal_mode = MEMORY")
    cursor.execute("PRAGMA cache_size = -256000")  # 256MB cache
    cursor.execute("PRAGMA temp_store = MEMORY")

    # Create table (no indexes yet - add after bulk insert for speed)
    cursor.execute("""
        CREATE TABLE proteins (
            protein_id TEXT PRIMARY KEY,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            strand INTEGER NOT NULL,
            contig_id TEXT
        )
    """)

    # Parse FAA and insert in chunks
    count = 0
    chunk = []

    with open(faa_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                parsed = _parse_prodigal_header(line)
                if parsed:
                    chunk.append(parsed)

                    if len(chunk) >= chunk_size:
                        cursor.executemany(
                            "INSERT OR IGNORE INTO proteins (protein_id, start, end, strand, contig_id) VALUES (?, ?, ?, ?, ?)",
                            chunk
                        )
                        conn.commit()
                        count += len(chunk)
                        chunk = []

                        if show_progress and count % 10000000 == 0:
                            elapsed = time.time() - start_time
                            rate = count / elapsed
                            print(f"  Processed {count:,} proteins... ({rate:,.0f}/sec)")

    # Insert remaining
    if chunk:
        cursor.executemany(
            "INSERT OR IGNORE INTO proteins (protein_id, start, end, strand, contig_id) VALUES (?, ?, ?, ?, ?)",
            chunk
        )
        conn.commit()
        count += len(chunk)

    elapsed = time.time() - start_time
    if show_progress:
        print(f"  Inserted {count:,} proteins in {elapsed/60:.1f} minutes")
        print("  Creating indexes...")

    # Create indexes after bulk insert (much faster)
    index_start = time.time()
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_contig ON proteins(contig_id)")
    conn.commit()

    if show_progress:
        index_time = time.time() - index_start
        print(f"  Index created in {index_time:.1f} seconds")

    # Optimize database
    cursor.execute("ANALYZE")
    conn.commit()
    conn.close()

    total_time = time.time() - start_time
    db_size = os.path.getsize(db_path) / (1024 * 1024 * 1024)

    if show_progress:
        print(f"Done! {count:,} proteins indexed in {total_time/60:.1f} minutes")
        print(f"Database size: {db_size:.2f} GB")

    return db_path


def build_index_from_gff(
    gff_path: str,
    db_path: str,
    chunk_size: int = 100000,
    show_progress: bool = True
) -> str:
    """
    Build a SQLite index from a Prodigal GFF file.

    Alternative to FAA parsing - GFF files contain the same position info.

    Args:
        gff_path: Path to Prodigal .gff file
        db_path: Path for output SQLite database
        chunk_size: Number of records per transaction
        show_progress: Whether to print progress

    Returns:
        Path to the created database
    """
    if not os.path.exists(gff_path):
        raise FileNotFoundError(f"GFF file not found: {gff_path}")

    if os.path.exists(db_path):
        os.remove(db_path)

    if show_progress:
        print(f"Building index from GFF: {gff_path} -> {db_path}")

    start_time = time.time()

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Optimize for bulk insert
    cursor.execute("PRAGMA synchronous = OFF")
    cursor.execute("PRAGMA journal_mode = MEMORY")
    cursor.execute("PRAGMA cache_size = -256000")
    cursor.execute("PRAGMA temp_store = MEMORY")

    cursor.execute("""
        CREATE TABLE proteins (
            protein_id TEXT PRIMARY KEY,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            strand INTEGER NOT NULL,
            contig_id TEXT
        )
    """)

    count = 0
    chunk = []

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'CDS':
                try:
                    contig_id = parts[0]
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = 1 if parts[6] == '+' else -1

                    # Extract ID from attributes
                    attrs = parts[8]
                    protein_id = None
                    for attr in attrs.split(';'):
                        if attr.startswith('ID='):
                            # Prodigal GFF ID format: ID=1_1 -> need to combine with contig
                            gene_id = attr[3:]
                            protein_id = f"{contig_id}_{gene_id.split('_')[-1]}"
                            break

                    if protein_id:
                        chunk.append((protein_id, start, end, strand, contig_id))

                        if len(chunk) >= chunk_size:
                            cursor.executemany(
                                "INSERT OR IGNORE INTO proteins VALUES (?, ?, ?, ?, ?)",
                                chunk
                            )
                            conn.commit()
                            count += len(chunk)
                            chunk = []

                            if show_progress and count % 10000000 == 0:
                                print(f"  Processed {count:,} proteins...")

                except (ValueError, IndexError):
                    continue

    if chunk:
        cursor.executemany(
            "INSERT OR IGNORE INTO proteins VALUES (?, ?, ?, ?, ?)",
            chunk
        )
        conn.commit()
        count += len(chunk)

    if show_progress:
        print(f"  Creating indexes...")

    cursor.execute("CREATE INDEX IF NOT EXISTS idx_contig ON proteins(contig_id)")
    cursor.execute("ANALYZE")
    conn.commit()
    conn.close()

    total_time = time.time() - start_time
    db_size = os.path.getsize(db_path) / (1024 * 1024 * 1024)

    if show_progress:
        print(f"Done! {count:,} proteins indexed in {total_time/60:.1f} minutes")
        print(f"Database size: {db_size:.2f} GB")

    return db_path


# Command-line interface
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Build SQLite index for Prodigal FAA/GFF files"
    )
    parser.add_argument(
        "input",
        help="Input Prodigal FAA or GFF file"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output SQLite database path (default: input.sqlite)"
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=100000,
        help="Records per transaction (default: 100000)"
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress progress output"
    )

    args = parser.parse_args()

    # Determine output path
    if args.output:
        db_path = args.output
    else:
        db_path = args.input + ".sqlite"

    # Determine file type and build index
    if args.input.endswith('.gff'):
        build_index_from_gff(
            args.input,
            db_path,
            chunk_size=args.chunk_size,
            show_progress=not args.quiet
        )
    else:
        build_prodigal_index(
            args.input,
            db_path,
            chunk_size=args.chunk_size,
            show_progress=not args.quiet
        )
