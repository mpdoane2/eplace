"""
Taxonomy extraction and sequence retrieval module.

This module provides functionality for extracting taxonomic information from BLAST results,
selecting representative sequences per taxonomic rank, and extracting sequences from databases.
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
from dataclasses import dataclass

from .blast_analysis import BlastHit

# Configure module logger
logger = logging.getLogger(__name__)


@dataclass
class TaxonomicInfo:
    """
    Represents taxonomic information for a sequence.
    
    Attributes:
        sequence_id: Sequence identifier
        taxid: NCBI taxonomy ID
        kingdom: Kingdom name
        phylum: Phylum name
        class_name: Class name
        order: Order name
        family: Family name
        genus: Genus name
        species: Species name
    """
    sequence_id: str
    taxid: Optional[str] = None
    kingdom: Optional[str] = None
    phylum: Optional[str] = None
    class_name: Optional[str] = None
    order: Optional[str] = None
    family: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None


class TaxonomyExtractor:
    """
    Class for extracting taxonomic information from sequence IDs.
    """
    
    VALID_RANKS = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    
    def __init__(self):
        """Initialize the TaxonomyExtractor."""
        pass
    
    @staticmethod
    def parse_sequence_id(sequence_id: str) -> TaxonomicInfo:
        """
        Parse taxonomic information from a sequence ID.
        
        For NCBI sequences, the ID may contain taxonomic information.
        This is a simplified parser that extracts basic information.
        
        Args:
            sequence_id: Sequence identifier
            
        Returns:
            TaxonomicInfo object
        """
        # For now, create a basic TaxonomicInfo object
        # In a real implementation, this would query NCBI taxonomy database
        # or parse information from the sequence header
        
        tax_info = TaxonomicInfo(sequence_id=sequence_id)
        
        # Extract taxonomy ID if present (format: gi|123456|...)
        parts = sequence_id.split('|')
        for i, part in enumerate(parts):
            if part == 'gi' and i + 1 < len(parts):
                tax_info.taxid = parts[i + 1]
                break
        
        return tax_info
    
    def extract_taxonomy_from_hits(
        self,
        hits: List[BlastHit]
    ) -> Dict[str, TaxonomicInfo]:
        """
        Extract taxonomic information from BLAST hits.
        
        Args:
            hits: List of BlastHit objects
            
        Returns:
            Dictionary mapping subject IDs to TaxonomicInfo objects
        """
        taxonomy_info = {}
        
        for hit in hits:
            if hit.subject_id not in taxonomy_info:
                tax_info = self.parse_sequence_id(hit.subject_id)
                taxonomy_info[hit.subject_id] = tax_info
        
        return taxonomy_info
    
    def group_hits_by_query(
        self,
        hits: List[BlastHit]
    ) -> Dict[str, List[BlastHit]]:
        """
        Group BLAST hits by query sequence.
        
        Args:
            hits: List of BlastHit objects
            
        Returns:
            Dictionary mapping query IDs to lists of hits
        """
        grouped = defaultdict(list)
        for hit in hits:
            grouped[hit.query_id].append(hit)
        return dict(grouped)
    
    def select_representatives_by_rank(
        self,
        hits: List[BlastHit],
        rank: str,
        max_per_rank: int = 1
    ) -> List[BlastHit]:
        """
        Select representative sequences per taxonomic rank.
        
        Args:
            hits: List of BlastHit objects for a single query
            rank: Taxonomic rank ('phylum', 'class', 'order', 'family', 'genus', 'species')
            max_per_rank: Maximum number of representatives per rank (default: 1)
            
        Returns:
            List of representative BlastHit objects
            
        Raises:
            ValueError: If rank is not valid
        """
        if rank not in self.VALID_RANKS:
            raise ValueError(
                f"Invalid rank '{rank}'. Must be one of: {', '.join(self.VALID_RANKS)}"
            )
        
        # For now, since we don't have full taxonomy information,
        # we'll use a simplified approach based on sequence ID patterns
        # In a real implementation, this would use NCBI taxonomy database
        
        # Group hits by taxonomic rank (using subject_id as proxy)
        rank_groups = defaultdict(list)
        
        for hit in hits:
            # Extract taxonomic information
            tax_info = self.parse_sequence_id(hit.subject_id)
            
            # Use subject_id prefix as a proxy for taxonomic grouping
            # In a real implementation, this would use actual taxonomy
            rank_key = self._get_rank_key(hit.subject_id, rank)
            rank_groups[rank_key].append(hit)
        
        # Select best representative from each rank
        representatives = []
        for rank_key, rank_hits in rank_groups.items():
            # Sort by bit score (best first)
            rank_hits.sort(key=lambda h: h.bit_score, reverse=True)
            
            # Take top N representatives
            representatives.extend(rank_hits[:max_per_rank])
        
        logger.info(
            f"Selected {len(representatives)} representative sequences "
            f"from {len(hits)} hits at rank '{rank}'"
        )
        
        return representatives
    
    @staticmethod
    def _get_rank_key(sequence_id: str, rank: str) -> str:
        """
        Get a rank key for grouping sequences.
        
        This is a simplified implementation that uses sequence ID patterns.
        In a real implementation, this would query NCBI taxonomy database.
        
        Args:
            sequence_id: Sequence identifier
            rank: Taxonomic rank
            
        Returns:
            Rank key for grouping
        """
        # Simple heuristic: use the first part of the ID
        # In reality, we'd need to query NCBI taxonomy database
        parts = sequence_id.split('|')
        if len(parts) > 1:
            return f"{rank}_{parts[0]}_{parts[1][:8]}"
        else:
            return f"{rank}_{sequence_id[:16]}"


class SequenceExtractor:
    """
    Class for extracting sequences from BLAST databases.
    """
    
    def __init__(self, blastdb_path: Optional[Path] = None):
        """
        Initialize the SequenceExtractor.
        
        Args:
            blastdb_path: Path to BLAST database directory. If None, uses BLASTDB env var.
        """
        self.blastdb_path = blastdb_path
        if self.blastdb_path is None:
            blastdb_env = os.environ.get('BLASTDB')
            if blastdb_env:
                self.blastdb_path = Path(blastdb_env)
            else:
                self.blastdb_path = Path.home() / "blastdb"
    
    def check_blastdbcmd_available(self) -> bool:
        """
        Check if blastdbcmd is available in the system.
        
        Returns:
            True if blastdbcmd is available, False otherwise
        """
        try:
            result = subprocess.run(
                ['blastdbcmd', '-version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError):
            return False
    
    def extract_sequences(
        self,
        sequence_ids: List[str],
        output_fasta: Path,
        database: str = "core_nt"
    ) -> bool:
        """
        Extract sequences from BLAST database using blastdbcmd.
        
        Args:
            sequence_ids: List of sequence IDs to extract
            output_fasta: Path to output FASTA file
            database: Name of BLAST database (default: "core_nt")
            
        Returns:
            True if extraction was successful, False otherwise
            
        Raises:
            RuntimeError: If blastdbcmd is not available
        """
        if not self.check_blastdbcmd_available():
            raise RuntimeError("blastdbcmd is not available. Please install BLAST+ tools.")
        
        if not sequence_ids:
            logger.warning("No sequence IDs provided for extraction")
            return False
        
        # Build database path
        db_path = self.blastdb_path / database
        
        # Create a temporary file with sequence IDs
        id_file = output_fasta.parent / f"{output_fasta.stem}_ids.txt"
        
        try:
            with open(id_file, 'w') as f:
                for seq_id in sequence_ids:
                    f.write(f"{seq_id}\n")
            
            # Run blastdbcmd
            cmd = [
                'blastdbcmd',
                '-db', str(db_path),
                '-entry_batch', str(id_file),
                '-out', str(output_fasta)
            ]
            
            logger.info(f"Extracting {len(sequence_ids)} sequences from database")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600  # 10 minute timeout
            )
            
            if result.returncode != 0:
                logger.error(f"blastdbcmd failed with error: {result.stderr}")
                return False
            
            logger.info(f"Sequences extracted successfully to {output_fasta}")
            return True
            
        except subprocess.TimeoutExpired:
            logger.error("Sequence extraction timed out")
            return False
        except Exception as e:
            logger.error(f"Error extracting sequences: {e}")
            return False
        finally:
            # Clean up temporary ID file
            if id_file.exists():
                id_file.unlink()
    
    def extract_representatives_for_query(
        self,
        query_id: str,
        representative_hits: List[BlastHit],
        output_dir: Path,
        database: str = "core_nt"
    ) -> Optional[Path]:
        """
        Extract representative sequences for a single query to a FASTA file.
        
        Args:
            query_id: Query sequence identifier
            representative_hits: List of representative BlastHit objects
            output_dir: Output directory for FASTA files
            database: Name of BLAST database
            
        Returns:
            Path to output FASTA file if successful, None otherwise
        """
        if not representative_hits:
            logger.warning(f"No representative hits for query {query_id}")
            return None
        
        # Create output directory if it doesn't exist
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate output filename
        safe_query_id = query_id.replace('|', '_').replace('/', '_')
        output_fasta = output_dir / f"{safe_query_id}_representatives.fasta"
        
        # Extract sequence IDs
        sequence_ids = [hit.subject_id for hit in representative_hits]
        
        # Extract sequences
        success = self.extract_sequences(
            sequence_ids=sequence_ids,
            output_fasta=output_fasta,
            database=database
        )
        
        if success:
            return output_fasta
        else:
            return None


def process_blast_results_for_taxonomy(
    blast_hits: List[BlastHit],
    output_dir: Path,
    rank: str = "species",
    database: str = "core_nt",
    blastdb_path: Optional[Path] = None
) -> Dict[str, Optional[Path]]:
    """
    Process BLAST hits to extract representative sequences per taxonomic rank.
    
    Args:
        blast_hits: List of BlastHit objects
        output_dir: Output directory for FASTA files
        rank: Taxonomic rank for representative selection
        database: Name of BLAST database
        blastdb_path: Path to BLAST database directory
        
    Returns:
        Dictionary mapping query IDs to output FASTA file paths
    """
    tax_extractor = TaxonomyExtractor()
    seq_extractor = SequenceExtractor(blastdb_path)
    
    # Group hits by query
    grouped_hits = tax_extractor.group_hits_by_query(blast_hits)
    
    # Process each query
    results = {}
    
    for query_id, query_hits in grouped_hits.items():
        logger.info(f"Processing query {query_id} with {len(query_hits)} hits")
        
        # Select representatives
        representatives = tax_extractor.select_representatives_by_rank(
            hits=query_hits,
            rank=rank
        )
        
        # Create query-specific output directory
        query_output_dir = output_dir / query_id.replace('|', '_').replace('/', '_')
        
        # Extract sequences
        output_fasta = seq_extractor.extract_representatives_for_query(
            query_id=query_id,
            representative_hits=representatives,
            output_dir=query_output_dir,
            database=database
        )
        
        results[query_id] = output_fasta
    
    return results
