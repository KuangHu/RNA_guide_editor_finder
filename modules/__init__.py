# RNA Guide Editor Finder Modules

from .sequence_search import (
    SequenceSearcher,
    SearchHit,
    write_results_fasta,
    DEFAULT_DATABASE_PATH,
    DEFAULT_FNA_FOLDER,
)
from .alignment_filter import (
    FilterEngine,
    FilterPipeline,
    FilterLevel,
    FilterConfig,
    Filters,
    create_default_pipeline,
    create_strict_pipeline,
    create_relaxed_pipeline,
    create_legacy_pipeline,
)
from .region_extractor import (
    RegionExtractor,
    extract_contig_from_protein_id,
    create_transposon_dict,
)
from .short_alignment_finder import (
    ShortAlignmentFinder,
    find_short_alignments,
    find_alignments_between_sequences,
    extend_alignment_with_gaps,
)
from .is_element_processor import (
    ISElementProcessor,
    process_is_element,
    process_is_elements_file,
    load_is_elements,
)
# Re-export utility functions from utils.parsers for convenience
from utils.parsers import (
    reverse_complement,
    parse_prodigal_faa,
)
