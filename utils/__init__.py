# Utility functions for RNA Guide Editor Finder

from .parsers import (
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
    # DataFrame conversion
    diamond_hits_to_dataframe,
    hmm_hits_to_dataframe,
    save_hits_to_csv,
    save_hits_to_excel,
    # Utilities
    detect_file_format,
)

from .prodigal_index import (
    # Data class
    ProteinPosition,
    # Index class for fast lookups
    ProdigalIndex,
    # Index builders
    build_prodigal_index,
    build_index_from_gff,
)
