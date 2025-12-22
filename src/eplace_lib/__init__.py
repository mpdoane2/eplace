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
    TaxonomicInfo,
    TaxonomyExtractor,
    SequenceExtractor,
    process_blast_results_for_taxonomy
)

__all__ = [
    "SequenceAnalyzer",
    "NCBIDownloader",
    "setup_ncbi_database",
    "BlastHit",
    "FastaReader",
    "BlastRunner",
    "run_blast_search",
    "TaxonomicInfo",
    "TaxonomyExtractor",
    "SequenceExtractor",
    "process_blast_results_for_taxonomy"
]
