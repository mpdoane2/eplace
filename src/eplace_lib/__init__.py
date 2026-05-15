"""
eplace_lib: A library for analyzing eDNA sequences

ePLACE (environmental Phylogenetic Localisation and Clade Estimation)
provides tools for analyzing environmental DNA sequences.
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("eplace")
except PackageNotFoundError:
    __version__ = "unknown"

__author__ = "Rob Edwards"

# Import main classes for convenient access
from .sequences import SequenceAnalyzer
from .ncbi_download import (
    NCBIDownloader,
    MMseqsDownloader,
    setup_ncbi_database,
    setup_mmseqs_database,
    setup_mmseqs_taxonomy,
    check_available_memory_gb
)
from .blast_analysis import (
    BlastHit,
    FastaReader,
    BlastRunner,
    run_blast_search,
    MMseqs2Runner,
    run_mmseqs_search,
    validate_mmseqs_memory_limit,
    normalize_sequence_id
)
from .taxonomy import (
    TaxonomyExtractor,
    SequenceExtractor,
    process_blast_results_for_taxonomy,
    generate_classification_summary
)
from .alignment import (
    SequenceTrimmer,
    MAFFTAligner,
    IQTreeBuilder,
    process_query_alignment_and_tree,
    process_query_alignment_and_tree_parallel,
    check_alignment_consistency,
    group_hits_by_group_rank,
    create_grouped_fasta_with_queries,
    trim_grouped_sequences,
    process_grouped_alignment_and_tree,
    process_grouped_alignment_and_tree_parallel,
    concatenate_all_groups_and_build_tree
)

__all__ = [
    "SequenceAnalyzer",
    "NCBIDownloader",
    "MMseqsDownloader",
    "setup_ncbi_database",
    "setup_mmseqs_database",
    "setup_mmseqs_taxonomy",
    "check_available_memory_gb",
    "BlastHit",
    "FastaReader",
    "BlastRunner",
    "run_blast_search",
    "MMseqs2Runner",
    "run_mmseqs_search",
    "validate_mmseqs_memory_limit",
    "normalize_sequence_id",
    "TaxonomyExtractor",
    "SequenceExtractor",
    "process_blast_results_for_taxonomy",
    "generate_classification_summary",
    "SequenceTrimmer",
    "MAFFTAligner",
    "IQTreeBuilder",
    "process_query_alignment_and_tree",
    "process_query_alignment_and_tree_parallel",
    "check_alignment_consistency",
    "group_hits_by_group_rank",
    "create_grouped_fasta_with_queries",
    "trim_grouped_sequences",
    "process_grouped_alignment_and_tree",
    "process_grouped_alignment_and_tree_parallel",
    "concatenate_all_groups_and_build_tree"
]
