"""
BLAST analysis module for sequence comparison.

This module provides functionality for running BLAST searches and filtering results
based on sequence identity and coverage criteria.
"""

import os
import subprocess
import logging
from typing import Optional
from pathlib import Path
from dataclasses import dataclass

# Configure module logger
logger = logging.getLogger(__name__)


@dataclass
class BlastHit:
    """
    Represents a single BLAST hit result.
    
    Attributes:
        query_id: Query sequence identifier
        subject_id: Subject (database) sequence identifier
        percent_identity: Percentage of identical matches
        alignment_length: Length of alignment
        query_length: Length of query sequence
        subject_length: Length of subject sequence
        query_start: Start position in query
        query_end: End position in query
        subject_start: Start position in subject
        subject_end: End position in subject
        evalue: Expectation value
        bit_score: Bit score
        query_coverage: Percentage of query covered by alignment
        subject_taxid: the subject's taxonomy id
        subject_taxids: blast should give a ; separated list of the hierarchy
        subject_rank_tid: the subjects taxonomy ID at our rank
        subject_rank_name: the subjects taxonomy name at our rank
        subject_phylum_tid: the taxonomy ID of the subject's phylum that we will later use for grouping
        subject_phylum_name: the name of the subject's phylum that we will later use for grouping
    """
    query_id: str
    subject_id: str
    percent_identity: float
    alignment_length: int
    query_length: int
    subject_length: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    bit_score: float
    query_coverage: float
    subject_taxid: str
    subject_taxids: str
    subject_rank_tid: Optional[str] = None
    subject_rank_name: Optional[str] = None
    subject_phylum_tid: Optional[str] = None
    subject_phylum_name: Optional[str] = None
    
    def get_accession(self) -> str:
        """
        Extract the accession number from the subject_id.
        
        BLAST IDs can be in various formats:
        - gi|2273658778|gb|MZ387488.1| -> MZ387488.1
        - ref|NZ_CP123456.1| -> NZ_CP123456.1
        - MZ387488.1 -> MZ387488.1 (already in accession format)
        
        Returns:
            The accession number extracted from subject_id
        """
        # If the ID contains pipes, it's in the gi|...|db|accession| format
        if '|' in self.subject_id:
            parts = self.subject_id.split('|')
            # Look for the accession, which is typically after a database identifier
            # Common formats: gi|123|gb|ACC.1|, ref|ACC.1|, gb|ACC.1|
            for i, part in enumerate(parts):
                if part in ['gb', 'ref', 'emb', 'dbj', 'pdb', 'prf', 'sp', 'tr', 'gnl']:
                    # Next part should be the accession
                    if i + 1 < len(parts):
                        return parts[i + 1]
            # If no database identifier found, return the last non-empty part
            non_empty = [p for p in parts if p]
            if non_empty:
                return non_empty[-1]
        
        # If no pipes, assume it's already an accession
        return self.subject_id

class FastaReader:
    """
    Class for reading FASTA files.
    """
    
    @staticmethod
    def read_fasta(fasta_path: Path) -> dict[str, str]:
        """
        Read sequences from a FASTA file.
        
        Args:
            fasta_path: Path to the FASTA file
            
        Returns:
            dictionary mapping sequence IDs to sequences
            
        Raises:
            FileNotFoundError: If FASTA file doesn't exist
            ValueError: If FASTA file is malformed
        """
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
        
        sequences = {}
        current_id = None
        current_seq = []
        
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    
                    # Start new sequence
                    # Remove '>' and take first word as ID
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Save last sequence
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
        
        if not sequences:
            raise ValueError(f"No sequences found in FASTA file: {fasta_path}")
        
        return sequences
    
    @staticmethod
    def get_sequence_lengths(fasta_path: Path) -> dict[str, int]:
        """
        Get the length of each sequence in a FASTA file.
        
        Args:
            fasta_path: Path to the FASTA file
            
        Returns:
            dictionary mapping sequence IDs to their lengths
        """
        sequences = FastaReader.read_fasta(fasta_path)
        return {seq_id: len(seq) for seq_id, seq in sequences.items()}


class BlastRunner:
    """
    Class for running BLAST searches and parsing results.
    """
    
    def __init__(self, blastdb_path: Optional[Path] = None):
        """
        Initialize the BlastRunner.
        
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
    
    def check_blastn_available(self) -> bool:
        """
        Check if blastn is available in the system.
        
        Returns:
            True if blastn is available, False otherwise
        """
        try:
            result = subprocess.run(
                ['blastn', '-version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError):
            return False
    
    def run_blastn(
        self,
        query_fasta: Path,
        output_file: Path,
        database: str = "core_nt",
        num_threads: int = 1,
        max_target_seqs: int = 100,
        evalue: float = 1e-5,
        outfmt: str = "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore staxid staxids"
    ) -> bool:
        """
        Run blastn search.
        
        Args:
            query_fasta: Path to query FASTA file
            output_file: Path to output file
            database: Name of BLAST database (default: "core_nt")
            num_threads: Number of threads to use
            max_target_seqs: Maximum number of target sequences to report
            evalue: E-value threshold
            outfmt: Output format string
            
        Returns:
            True if BLAST ran successfully, False otherwise
            
        Raises:
            FileNotFoundError: If query file doesn't exist
            RuntimeError: If blastn is not available
        """
        if not query_fasta.exists():
            raise FileNotFoundError(f"Query FASTA file not found: {query_fasta}")
        
        if not self.check_blastn_available():
            raise RuntimeError("blastn is not available. Please install BLAST+ tools.")
        
        # Build database path
        db_path = self.blastdb_path / database
        
        # Build blastn command
        cmd = [
            'blastn',
            '-query', str(query_fasta),
            '-db', str(db_path),
            '-out', str(output_file),
            '-outfmt', outfmt,
            '-num_threads', str(num_threads),
            '-max_target_seqs', str(max_target_seqs),
            '-evalue', str(evalue)
        ]
        
        logger.info(f"Running BLAST search: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )
            
            if result.returncode != 0:
                logger.error(f"BLAST failed with error: {result.stderr}")
                return False
            
            logger.info(f"BLAST search completed successfully. Output: {output_file}")
            return True
            
        except subprocess.TimeoutExpired:
            logger.error("BLAST search timed out")
            return False
        except Exception as e:
            logger.error(f"Error running BLAST: {e}")
            return False
    
    def parse_blast_results(
        self,
        blast_output: Path,
        query_lengths: Optional[dict[str, int]] = None
    ) -> list[BlastHit]:
        """
        Parse BLAST tabular output.
        
        Args:
            blast_output: Path to BLAST output file (tabular format)
            query_lengths: dictionary of query sequence lengths. If None, uses qlen from results.
            
        Returns:
            list of BlastHit objects
            
        Raises:
            FileNotFoundError: If BLAST output file doesn't exist
            ValueError: If BLAST output is malformed
        """
        if not blast_output.exists():
            raise FileNotFoundError(f"BLAST output file not found: {blast_output}")
        
        hits = []
        
        with open(blast_output, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 12:
                    raise ValueError(
                        f"Invalid BLAST output format at line {line_num}: "
                        f"expected at least 12 fields, got {len(fields)}"
                    )
                
                try:
                    query_id = fields[0]
                    subject_id = fields[1]
                    percent_identity = float(fields[2])
                    alignment_length = int(fields[3])
                    query_length = int(fields[4])
                    subject_length = int(fields[5])
                    query_start = int(fields[6])
                    query_end = int(fields[7])
                    subject_start = int(fields[8])
                    subject_end = int(fields[9])
                    evalue = float(fields[10])
                    bit_score = float(fields[11])
                    staxid = fields[12]
                    staxids = fields[13]
                    
                    # Calculate query coverage
                    query_coverage = (abs(query_end - query_start) + 1) / query_length * 100
                    
                    hit = BlastHit(
                        query_id=query_id,
                        subject_id=subject_id,
                        percent_identity=percent_identity,
                        alignment_length=alignment_length,
                        query_length=query_length,
                        subject_length=subject_length,
                        query_start=query_start,
                        query_end=query_end,
                        subject_start=subject_start,
                        subject_end=subject_end,
                        evalue=evalue,
                        bit_score=bit_score,
                        query_coverage=query_coverage,
                        subject_taxid=staxid,
                        subject_taxids=staxids
                    )
                    hits.append(hit)
                    
                except (ValueError, IndexError) as e:
                    raise ValueError(
                        f"Error parsing BLAST output at line {line_num}: {e}"
                    )
        
        return hits
    
    def filter_blast_hits(
        self,
        hits: list[BlastHit],
        min_identity: float = 90.0,
        min_coverage: float = 80.0,
        min_alignment_length: Optional[int] = None
    ) -> list[BlastHit]:
        """
        Filter BLAST hits based on identity and coverage thresholds.
        
        Args:
            hits: list of BlastHit objects
            min_identity: Minimum percent identity (default: 90.0)
            min_coverage: Minimum query coverage percentage (default: 80.0)
            min_alignment_length: Minimum alignment length (optional)
            
        Returns:
            Filtered list of BlastHit objects
        """
        filtered_hits = []
        
        for hit in hits:
            if hit.percent_identity < min_identity:
                continue
            if hit.query_coverage < min_coverage:
                continue
            if min_alignment_length and hit.alignment_length < min_alignment_length:
                continue
            
            filtered_hits.append(hit)
        
        logger.info(
            f"Filtered {len(hits)} hits to {len(filtered_hits)} hits "
            f"(min_identity={min_identity}%, min_coverage={min_coverage}%)"
        )
        
        return filtered_hits


def run_blast_search(
    query_fasta: Path,
    output_file: Path,
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
    database: str = "core_nt",
    blastdb_path: Optional[Path] = None,
    num_threads: int = 1,
    skip_existing: bool = True
) -> tuple[bool, list[BlastHit]]:
    """
    Convenience function to run BLAST search and return filtered hits.
    
    Args:
        query_fasta: Path to query FASTA file
        output_file: Path to output file
        min_identity: Minimum percent identity (default: 90.0)
        min_coverage: Minimum query coverage percentage (default: 80.0)
        database: Name of BLAST database (default: "core_nt")
        blastdb_path: Path to BLAST database directory
        num_threads: Number of threads to use
        
    Returns:
        Tuple of (success: bool, filtered_hits: list[BlastHit])
    """

    runner = BlastRunner(blastdb_path)

    if os.path.exists(output_file) and skip_existing:
        logger.info(f"The blast output file {output_file} already exists. Skipping and using these results")
    else:    
        # Run BLAST
        success = runner.run_blastn(
            query_fasta=query_fasta,
            output_file=output_file,
            database=database,
            num_threads=num_threads
        )
        
        if not success:
            return False, []
        
    # Parse results
    try:
        hits = runner.parse_blast_results(output_file)
        
        # Filter results
        filtered_hits = runner.filter_blast_hits(
            hits,
            min_identity=min_identity,
            min_coverage=min_coverage
        )
        
        return True, filtered_hits
        
    except Exception as e:
        logger.error(f"Error processing BLAST results: {e}")
        return False, []
