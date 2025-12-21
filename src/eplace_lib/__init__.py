"""
eplace_lib: A library for analyzing eDNA sequences

ePLACE (environmental Phylogenetic Localisation and Clade Estimation)
provides tools for analyzing environmental DNA sequences.
"""

__version__ = "0.1.0"
__author__ = "Rob Edwards"

# Import main classes for convenient access
from .sequences import SequenceAnalyzer

__all__ = ["SequenceAnalyzer"]
