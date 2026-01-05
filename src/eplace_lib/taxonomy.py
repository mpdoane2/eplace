"""
Taxonomy extraction and sequence retrieval module.

This module provides functionality for extracting taxonomic information from BLAST results,
selecting representative sequences per taxonomic rank, and extracting sequences from databases.
"""

import os
import subprocess
import logging
import sys
from pathlib import Path
from typing import Optional, List, Dict
from collections import defaultdict

from .blast_analysis import BlastHit

import pytaxonkit

# Configure module logger
logger = logging.getLogger(__name__)


class TaxonomyExtractor:
    """
    Class for extracting taxonomic information from sequence IDs.
    """

    def parse_taxids(self, tax_ids: list[str]) -> dict[str, dict[str, tuple[str, str]]]:
        """
        Parse taxonomic information from the taxonomy IDs from the BLAST hits

        Args:
            tax_ids: the taxonomy IDs reported by BLAST
            
        Returns:
           dictionary containing the rank and a tuple of the taxonomy ID and the name
        """
        # make sure that duplicate taxids are removed before we look them up
        tax_ids = list(set(tax_ids))
        VALID_RANKS = ['phylum', 'class', 'order', 'family', 'genus', 'species']
        taxonomy_dict = {}

        # we need to get the whole lineage, and then convert it to a dict
        try:
            df = pytaxonkit.lineage(tax_ids)
        except Exception:
            logger.exception(f"Error retrieving taxonomic lineages")
            sys.exit(1)

        df['names'] = df['FullLineage'].str.split(';')
        df['taxids'] = df['FullLineageTaxIDs'].str.split(';')
        df['ranks'] = df['FullLineageRanks'].str.split(';')
        long_df = df.explode(['names', 'taxids', 'ranks'])
        filtered = long_df[long_df['ranks'].isin(VALID_RANKS)]

        for tid, rank, taxid, name in (
                filtered[['TaxID', 'ranks', 'taxids', 'names']]
                        .drop_duplicates()
                        .itertuples(index=False, name=None)
        ):
            tid = str(tid)
            taxid = str(taxid)

            taxonomy_dict.setdefault(tid, {})[rank] = (taxid, name)
        return taxonomy_dict
    
    def group_hits_by_query(
        self,
        hits: list[BlastHit]
    ) -> dict[str, list[BlastHit]]:
        """
        Group BLAST hits by query sequence.
        
        Args:
            hits: list of BlastHit objects
            
        Returns:
            dictionary mapping query IDs to lists of hits
        """
        grouped = defaultdict(list)
        for hit in hits:
            grouped[hit.query_id].append(hit)
        return dict(grouped)
    
    def select_representatives_by_rank(
        self,
        hits: list[BlastHit],
        rank: str,
        max_per_rank: int = 1,
        preferred_representatives: Optional[Dict[str, str]] = None
    ) -> list[BlastHit]:
        """
        Select representative sequences per taxonomic rank.
        
        Args:
            hits: list of BlastHit objects for a single query
            rank: Taxonomic rank for representative selection
            max_per_rank: Maximum number of representatives per rank (default: 1)
            preferred_representatives: Optional dictionary mapping rank_tid to preferred subject_id
                                      to ensure consistent representatives across queries
            
        Returns:
            list of representative BlastHit objects
        """
        
        if preferred_representatives is None:
            preferred_representatives = {}
        
        # Group hits by taxonomic rank (using subject_id as proxy)
        rank_groups = defaultdict(list)
        
        reported_hits = set()
        for hit in hits:
            if not hit.subject_taxonomy:
                logger.info(f"Skipping the taxonomy for {hit.query_id} as no subject taxonomy")
                continue
            if not hit.subject_taxonomy[rank]:
                logger.warning(
                    f"Hit {hit.subject_id} for query {hit.query_id} has no taxonomic information at rank {rank}")
                continue

            if isinstance(hit.subject_taxonomy[rank], tuple):
                # Log the first time we see each rank name
                if hit.subject_taxonomy[rank][1] not in reported_hits:
                    logger.info(f"Found a hit for {hit.query_id} at rank {rank}: {hit.subject_taxonomy[rank][1]} ({hit.subject_taxonomy[rank][0]})")
                    reported_hits.add(hit.subject_taxonomy[rank][1])
                # Add all hits with taxonomic information to rank_groups
                rank_groups[hit.subject_taxonomy[rank][1]].append(hit)
            else:
                logger.warning(f"Not really sure what {hit.subject_taxonomy[rank]} of type {type(hit.subject_taxonomy[rank])} is supposed to be")
        
        # Select best representative from each rank
        representatives = []
        for rank_key, rank_hits in rank_groups.items():
            # Check if we have a preferred representative for this rank
            preferred_subject_id = preferred_representatives.get(rank_key)
            
            if preferred_subject_id:
                # Look for the preferred representative in the current hits
                preferred_hit = next(
                    (hit for hit in rank_hits if hit.subject_id == preferred_subject_id),
                    None
                )
                
                if preferred_hit:
                    # Use the preferred representative
                    logger.info(f"Reusing previously selected representative {preferred_subject_id} for rank {rank_key}")
                    representatives.append(preferred_hit)
                    continue
            
            # No preferred representative or it's not in current hits
            # Sort by bit score (best first) and select new representative
            rank_hits.sort(key=lambda h: h.bit_score, reverse=True)
            
            # Take top N representatives
            representatives.extend(rank_hits[:max_per_rank])
        
        logger.info(
            f"Selected {len(representatives)} representative sequences from {len(hits)} hits at rank '{rank}'"
        )
        
        return representatives
    

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
        sequence_ids: list[str],
        output_fasta: Path,
        database: str = "core_nt"
    ) -> bool:
        """
        Extract sequences from BLAST database using blastdbcmd.
        
        Args:
            sequence_ids: list of sequence IDs to extract
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
            logger.error(f"Error extracting sequences (taxonomy): {e}")
            return False
        finally:
            # Clean up temporary ID file
            if id_file.exists():
                id_file.unlink()
    
    def extract_representatives_for_query(
        self,
        query_id: str,
        representative_hits: list[BlastHit],
        output_dir: Path,
        database: str = "core_nt"
    ) -> Optional[Path]:
        """
        Extract representative sequences for a single query to a FASTA file.
        
        Args:
            query_id: Query sequence identifier
            representative_hits: list of representative BlastHit objects
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


def rewrite_blast_hits(
    blast_hits: List[BlastHit],
    output_file: Path,
    header: bool = True) -> bool:
    """
    Rewrite the blast hits when we have annotated them

    Args:
        blast_hits: list of BlastHit objects
        output_file: the file to write to
        header: whether to include a header line in the file
    
    Returns:
        True on success
    """

    fields = [
        "query_id", "subject_id", "percent_identity", "alignment_length",
        "query_length", "subject_length", "query_start", "query_end",
        "subject_start", "subject_end", "evalue", "bit_score",
        "query_coverage", "subject_taxid", "subject_taxids",
        "subject_taxonomy"
    ]

    with open(output_file, 'w') as out:
        if header:
            print("\t".join(fields), file=out)

        for hit in blast_hits:
            print(
                "\t".join(
                    "" if getattr(hit, f) is None else str(getattr(hit, f))
                    for f in fields
                ),
                file=out
            )
    
    return True


def process_blast_results_for_taxonomy(
    blast_hits: List[BlastHit],
    output_dir: Path,
    rank: str = "genus",
    group_rank: str = "class",
    database: str = "core_nt",
    blastdb_path: Optional[Path] = None
) -> Dict[str, Optional[Path]]:
    """
    Process BLAST hits to extract representative sequences per taxonomic rank.
    
    Args:
        blast_hits: list of BlastHit objects
        output_dir: Output directory for FASTA files
        rank: Taxonomic rank for representative selection
        database: Name of BLAST database
        blastdb_path: Path to BLAST database directory
        
    Returns:
        dictionary mapping query IDs to output FASTA file paths
    """
    
    VALID_RANKS = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    if rank not in VALID_RANKS:
        raise ValueError(f"Rank: {rank} is not a valid rank. It must be one of: {VALID_RANKS}")
    if group_rank not in VALID_RANKS:
        raise ValueError(f"Grouping Rank: {group_rank} is not a valid rank. It must be one of: {VALID_RANKS}")
    
    tax_extractor = TaxonomyExtractor()
    seq_extractor = SequenceExtractor(blastdb_path)
    
    # get all the taxonomies
    subject_taxids = {hit.subject_taxid for hit in blast_hits}
    tax_dict = tax_extractor.parse_taxids(list(subject_taxids))

    # add all the ranks to all the hits
    for h in blast_hits:
        h.subject_taxonomy = tax_dict.get(h.subject_taxid)

    # Group hits by query
    grouped_hits = tax_extractor.group_hits_by_query(blast_hits)
    
    # Track selected representatives across queries to ensure consistency
    # Maps rank_tid -> subject_id of the selected representative
    preferred_representatives = {}
    
    # Process each query
    results = {}
    
    for query_id, query_hits in grouped_hits.items():
        logger.info(f"Processing query {query_id} with {len(query_hits)} hits")
        
        # Select representatives, preferring previously selected ones
        representatives = tax_extractor.select_representatives_by_rank(
            hits=query_hits,
            rank=rank,
            preferred_representatives=preferred_representatives
        )
        if len(representatives) == 0:
            logger.warning(f"Error: No representative sequences for {query_id} at rank {rank}")
            continue
        
        # Update the preferred representatives with newly selected ones
        for rep in representatives:
            if isinstance(rep.subject_taxonomy[rank], tuple) and rep.subject_taxonomy[rank][1] not in preferred_representatives:
                preferred_representatives[rep.subject_taxonomy[rank][1]] = rep.subject_id
                logger.info(f"Recording {rep.subject_id} as representative for rank {rep.subject_taxonomy[rank][1]}")
        
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
