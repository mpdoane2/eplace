"""
Sequence analysis module for eDNA data.

This module provides basic functionality for analyzing environmental DNA sequences.
"""

from typing import Dict, List, Optional


class SequenceAnalyzer:
    """
    A class for analyzing eDNA sequences.
    
    This class provides methods for basic sequence analysis operations
    commonly used in environmental DNA studies.
    """
    
    def __init__(self):
        """Initialize the SequenceAnalyzer."""
        pass
    
    def calculate_gc_content(self, sequence: str) -> float:
        """
        Calculate the GC content of a DNA sequence.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            GC content as a percentage (0-100)
        """
        if not sequence:
            return 0.0
        
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100
    
    def reverse_complement(self, sequence: str) -> str:
        """
        Calculate the reverse complement of a DNA sequence.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Reverse complement of the input sequence
        """
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                     'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
        
        rev_comp = ''.join(complement.get(base, base) for base in reversed(sequence))
        return rev_comp
    
    def count_bases(self, sequence: str) -> Dict[str, int]:
        """
        Count the occurrence of each base in a DNA sequence.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Dictionary with base counts
        """
        sequence = sequence.upper()
        return {
            'A': sequence.count('A'),
            'T': sequence.count('T'),
            'G': sequence.count('G'),
            'C': sequence.count('C'),
            'N': sequence.count('N')
        }
