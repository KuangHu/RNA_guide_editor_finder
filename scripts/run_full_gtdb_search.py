#!/usr/bin/env python3
"""
Full GTDB Search with Maximum Optimization
- Uses indexed MMseqs database
- Pre-loads database into RAM
- Uses all available CPU cores
- Optimized sensitivity and search parameters
"""

import os
import sys
import time
import subprocess
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from modules.sequence_search import SequenceSearcher, write_results_fasta

# IS110 transposon query sequence
QUERY_SEQUENCE = "tatccctccagtgcagagaaaatcggccagttttctctgcctgcagtccgcatgccgtatcgggccttgggttctaacctgttgcgtagatttatgcagcggactgcctttctcccaaagtgataaaccggacagtatcatggaccggttttcccggtaatccgtatttgcaaggttggtttcact"

# Organized GTDB paths
GTDB_ORGANIZED_DIR = "/groups/rubin/projects/kuang/db/GTDB_organized"
MMSEQS_DB = f"{GTDB_ORGANIZED_DIR}/mmseqs_db/gtdb_mmseq_db"
FNA_FOLDER = f"{GTDB_ORGANIZED_DIR}/genomes"  # symlink to actual GTDB folder
OUTPUT_FILE = "gtdb_full_search_results.fasta"
LOG_FILE = "gtdb_full_search.log"

# Maximum CPU optimization
MAX_THREADS = 48  # All available cores
SENSITIVITY = 5.0  # Balanced sensitivity for thorough search
MIN_IDENTITY = 0.8

def log(message):
    """Log message to both console and file"""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    log_msg = f"[{timestamp}] {message}"
    print(log_msg)
    with open(LOG_FILE, 'a') as f:
        f.write(log_msg + "\n")

def check_index_exists():
    """Check if MMseqs database index exists"""
    idx_file = f"{MMSEQS_DB}.idx"
    if os.path.exists(idx_file):
        log(f"✓ Index file found: {idx_file}")
        return True
    else:
        log(f"✗ Index file not found: {idx_file}")
        log("  Run 'mmseqs createindex' first!")
        return False

def preload_database():
    """Pre-load database into system RAM (Linux page cache)"""
    log("=" * 70)
    log("Step 1: Pre-loading database into RAM (touchdb)")
    log("This will load ~322 GB into RAM for maximum speed...")
    log("=" * 70)

    start_time = time.time()

    cmd = ["mmseqs", "touchdb", MMSEQS_DB, "--threads", str(MAX_THREADS)]
    log(f"Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        elapsed = time.time() - start_time
        log(f"✓ Database pre-loaded into RAM in {elapsed:.1f} seconds")
        log(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        log(f"✗ Failed to pre-load database: {e}")
        log(f"stderr: {e.stderr}")
        log("Continuing anyway (search will be slower)...")
        return False

def run_search():
    """Run the full GTDB search with optimizations"""
    log("=" * 70)
    log("Step 2: Running full GTDB search")
    log(f"Query: {QUERY_SEQUENCE[:50]}... ({len(QUERY_SEQUENCE)}bp)")
    log(f"Database: {MMSEQS_DB}")
    log(f"Threads: {MAX_THREADS}")
    log(f"Sensitivity: {SENSITIVITY}")
    log(f"Min Identity: {MIN_IDENTITY}")
    log("=" * 70)

    start_time = time.time()

    # Initialize searcher
    searcher = SequenceSearcher(
        database_path=MMSEQS_DB,
        fna_folder=FNA_FOLDER,
        min_seq_identity=MIN_IDENTITY,
        sensitivity=SENSITIVITY,
        threads=MAX_THREADS
    )

    # Run search
    log("Searching GTDB database...")
    results = searcher.search(
        sequence=QUERY_SEQUENCE,
        query_id="IS110_transposon",
        context_bp=50,  # Add 50bp flanking context
        dedup_threshold=0.95,
        cleanup=True
    )

    elapsed = time.time() - start_time

    log(f"✓ Search completed in {elapsed:.1f} seconds")
    log(f"Found {len(results)} unique hits after deduplication")

    return results

def save_results(results):
    """Save results to FASTA file"""
    log("=" * 70)
    log("Step 3: Saving results")
    log("=" * 70)

    if not results:
        log("No results to save")
        return

    write_results_fasta(results, OUTPUT_FILE, include_metadata=True)

    # Get file size
    size_mb = os.path.getsize(OUTPUT_FILE) / (1024 * 1024)
    log(f"✓ Wrote {len(results)} sequences to {OUTPUT_FILE} ({size_mb:.2f} MB)")

    # Show top hits
    log("\nTop 5 hits:")
    for i, hit in enumerate(results[:5], 1):
        log(f"  {i}. {hit.target_id}")
        log(f"     Identity: {hit.sequence_identity:.1%}, E-value: {hit.evalue:.2e}, Score: {hit.bit_score:.1f}")

def main():
    """Main workflow"""
    log("=" * 70)
    log("FULL GTDB SEARCH - MAXIMUM OPTIMIZATION MODE")
    log("=" * 70)
    log(f"CPU cores: {MAX_THREADS}")
    log(f"Available RAM: Check with 'free -h'")
    log("")

    # Check index
    if not check_index_exists():
        log("ERROR: Index not ready. Wait for 'mmseqs createindex' to complete.")
        sys.exit(1)

    # Pre-load database into RAM
    preload_database()

    # Run search
    results = run_search()

    # Save results
    save_results(results)

    log("=" * 70)
    log("COMPLETE!")
    log("=" * 70)
    log(f"Results: {OUTPUT_FILE}")
    log(f"Log: {LOG_FILE}")

if __name__ == "__main__":
    main()
