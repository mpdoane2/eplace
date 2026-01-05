"""
eplace_lib: A library for analyzing eDNA sequences

ePLACE (environmental Phylogenetic Localisation and Clade Estimation)
provides tools for analyzing environmental DNA sequences.
"""

__version__ = "0.1.0"
__author__ = "Rob Edwards"

# Import main classes for convenient access
from .sequences import SequenceAnalyzer
from .ncbi_download import NCBIDownloader, setup_ncbi_database
from .blast_analysis import (
    BlastHit,
    FastaReader,
    BlastRunner,
    run_blast_search
)
from .taxonomy import (
    TaxonomyExtractor,
    SequenceExtractor,
    process_blast_results_for_taxonomy
)
from .alignment import (
    SequenceTrimmer,
    MAFFTAligner,
    IQTreeBuilder,
    process_query_alignment_and_tree,
    check_alignment_consistency,
    group_hits_by_group_rank,
    create_grouped_fasta_with_queries,
    trim_grouped_sequences,
    process_grouped_alignment_and_tree
)

__all__ = [
    "SequenceAnalyzer",
    "NCBIDownloader",
    "setup_ncbi_database",
    "BlastHit",
    "FastaReader",
    "BlastRunner",
    "run_blast_search",
    "TaxonomyExtractor",
    "SequenceExtractor",
    "process_blast_results_for_taxonomy",
    "SequenceTrimmer",
    "MAFFTAligner",
    "IQTreeBuilder",
    "process_query_alignment_and_tree",
    "check_alignment_consistency",
    "group_hits_by_group_rank",
    "create_grouped_fasta_with_queries",
    "trim_grouped_sequences",
    "process_grouped_alignment_and_tree"
]
