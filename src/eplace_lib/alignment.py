"""
Sequence alignment and phylogenetic tree building module.

This module provides functionality for trimming sequences based on BLAST alignments,
aligning sequences using MAFFT, and building phylogenetic trees using IQTree.
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import Optional, List, Dict, Tuple
from dataclasses import dataclass

from .blast_analysis import BlastHit, FastaReader

# Configure module logger
logger = logging.getLogger(__name__)


class SequenceTrimmer:
    """
    Class for trimming sequences based on BLAST alignment coordinates.
    """
    
    @staticmethod
    def trim_sequence_by_coordinates(
        sequence: str,
        start: int,
        end: int
    ) -> str:
        """
        Trim a sequence to extract the region between start and end coordinates.
        
        BLAST coordinates are 1-indexed, so we need to adjust for Python's 0-indexing.
        
        Args:
            sequence: The full sequence string
            start: Start position (1-indexed, inclusive)
            end: End position (1-indexed, inclusive)
            
        Returns:
            Trimmed sequence string
        """
        # Convert 1-indexed BLAST coordinates to 0-indexed Python
        # Also handle reverse complement alignments (start > end)
        if start > end:
            # Reverse strand alignment
            python_start = end - 1
            python_end = start
        else:
            # Forward strand alignment
            python_start = start - 1
            python_end = end
        
        # Ensure coordinates are within bounds
        python_start = max(0, python_start)
        python_end = min(len(sequence), python_end)
        
        return sequence[python_start:python_end]
    
    @staticmethod
    def trim_sequences_from_blast_hits(
        fasta_path: Path,
        blast_hits: List[BlastHit],
        output_fasta: Path,
        query_id: str
    ) -> bool:
        """
        Trim sequences in a FASTA file based on BLAST hit coordinates.
        
        This reads the representative sequences, trims them to the aligned regions,
        and writes them to a new FASTA file along with the query sequence.
        
        Args:
            fasta_path: Path to input FASTA file with full-length sequences
            blast_hits: List of BlastHit objects for this query
            output_fasta: Path to output FASTA file with trimmed sequences
            query_id: The query sequence ID to include in output
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Read all sequences from the input FASTA
            sequences = FastaReader.read_fasta(fasta_path)
            
            # Create a mapping of subject accession to blast hits for quick lookup
            # The FASTA file will have accessions (e.g., MZ387488.1) not full IDs
            hit_map = {hit.get_accession(): hit for hit in blast_hits}
            
            # Open output file
            with open(output_fasta, 'w') as out:
                # First, write the query sequence if it exists
                if query_id in sequences:
                    query_seq = sequences[query_id]
                    out.write(f">{query_id}\n")
                    # Write sequence in lines of 60 characters
                    for i in range(0, len(query_seq), 60):
                        out.write(query_seq[i:i+60] + "\n")
                    logger.info(f"Added query sequence {query_id} ({len(query_seq)} bp)")
                
                # Now process subject sequences
                for seq_id, sequence in sequences.items():
                    if seq_id == query_id:
                        continue  # Skip query, already written
                    
                    # Find the corresponding BLAST hit by accession
                    # The seq_id from FASTA should match the accession from the hit
                    hit = hit_map.get(seq_id)
                    
                    if hit is None:
                        # This sequence doesn't have a BLAST hit, skip it
                        logger.warning(f"No BLAST hit found for sequence {seq_id}, skipping")
                        continue
                    
                    # Trim the sequence based on subject coordinates
                    trimmed_seq = SequenceTrimmer.trim_sequence_by_coordinates(
                        sequence,
                        hit.subject_start,
                        hit.subject_end
                    )
                    
                    # Write trimmed sequence with taxonomic information in header
                    header = seq_id
                    if hit.subject_rank_name:
                        header = f"{seq_id} {hit.subject_rank_name}"
                    
                    out.write(f">{header}\n")
                    # Write sequence in lines of 60 characters
                    for i in range(0, len(trimmed_seq), 60):
                        out.write(trimmed_seq[i:i+60] + "\n")
                    
                    logger.info(
                        f"Trimmed {seq_id} from {len(sequence)} bp to {len(trimmed_seq)} bp "
                        f"(coords: {hit.subject_start}-{hit.subject_end})"
                    )
            
            return True
            
        except Exception as e:
            logger.error(f"Error trimming sequences: {e}")
            return False


class MAFFTAligner:
    """
    Class for running MAFFT sequence alignments.
    """
    
    @staticmethod
    def check_mafft_available() -> bool:
        """
        Check if MAFFT is available in the system.
        
        Returns:
            True if MAFFT is available, False otherwise
        """
        try:
            result = subprocess.run(
                ['mafft', '--version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                logging.info(f"Using mafft to build the alignment ({result.stdout})")
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError):
            return False
    
    @staticmethod
    def align_sequences(
        input_fasta: Path,
        output_fasta: Path,
        auto_orient: bool = True,
        num_threads: int = 1
    ) -> bool:
        """
        Align sequences using MAFFT.
        
        Args:
            input_fasta: Path to input FASTA file with sequences to align
            output_fasta: Path to output aligned FASTA file
            auto_orient: Use MAFFT's auto-orient feature (default: True)
            num_threads: Number of threads to use
            
        Returns:
            True if alignment was successful, False otherwise
        """
        if not MAFFTAligner.check_mafft_available():
            logger.error("MAFFT is not available. Please install MAFFT.")
            return False
        
        if not input_fasta.exists():
            logger.error(f"Input FASTA file not found: {input_fasta}")
            return False
        
        # Build MAFFT command
        cmd = ['mafft']
        
        # Add auto-orient option
        if auto_orient:
            cmd.append('--adjustdirection')
        
        # Add threading
        cmd.extend(['--thread', str(num_threads)])
        
        # Add input file
        cmd.append(str(input_fasta))
        
        logger.info(f"Running MAFFT alignment: {' '.join(cmd)}")
        
        try:
            with open(output_fasta, 'w') as out:
                result = subprocess.run(
                    cmd,
                    stdout=out,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=3600  # 1 hour timeout
                )
            
            if result.returncode != 0:
                logger.error(f"MAFFT failed with error: {result.stderr}")
                return False
            
            logger.info(f"MAFFT alignment completed successfully. Output: {output_fasta}")
            return True
            
        except subprocess.TimeoutExpired:
            logger.error("MAFFT alignment timed out")
            return False
        except Exception as e:
            logger.error(f"Error running MAFFT: {e}")
            return False


class IQTreeBuilder:
    """
    Class for building phylogenetic trees using IQTree.
    """
    
    @staticmethod
    def check_iqtree_available() -> Tuple[bool, Optional[str]]:
        """
        Check if IQTree is available in the system.
        
        Returns:
            Tuple of (available: bool, command: str or None)
        """
        # Try 'iqtree', 'iqtree2', and 'iqtree3' commands
        for cmd in ['iqtree3', 'iqtree2', 'iqtree']:
            try:
                result = subprocess.run(
                    [cmd, '--version'],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0:
                    logging.info(f"Using {cmd} to build the trees ({result.stdout})")
                    return True, cmd
            except (subprocess.SubprocessError, FileNotFoundError):
                continue
        
        return False, None
    
    @staticmethod
    def build_tree(
        alignment_fasta: Path,
        output_prefix: Path,
        num_threads: int = 1,
        model: str = "MFP"
    ) -> bool:
        """
        Build a phylogenetic tree using IQTree.
        
        Args:
            alignment_fasta: Path to aligned FASTA file
            output_prefix: Prefix for output files
            num_threads: Number of threads to use
            model: Substitution model (default: "MFP" for automatic ModelFinder Plus selection)
            
        Returns:
            True if tree building was successful, False otherwise
        """
        available, iqtree_cmd = IQTreeBuilder.check_iqtree_available()
        if not available:
            logger.error("IQTree is not available. Please install IQTree or IQTree2.")
            return False
        
        if not alignment_fasta.exists():
            logger.error(f"Alignment file not found: {alignment_fasta}")
            return False
        
        # Build IQTree command
        cmd = [
            iqtree_cmd,
            '-s', str(alignment_fasta),
            '-pre', str(output_prefix),
            '-m', model,
            '-T', str(num_threads)
        ]
        
        logger.info(f"Running IQTree: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )
            
            if result.returncode != 0:
                logger.error(f"IQTree failed with error: {result.stderr}")
                return False
            
            # Check if tree file was created
            tree_file = Path(str(output_prefix) + ".treefile")
            if not tree_file.exists():
                logger.error("IQTree did not produce a tree file")
                return False
            
            logger.info(f"IQTree completed successfully. Tree: {tree_file}")
            return True
            
        except subprocess.TimeoutExpired:
            logger.error("IQTree timed out")
            return False
        except Exception as e:
            logger.error(f"Error running IQTree: {e}")
            return False
    
    @staticmethod
    def relabel_tree_with_taxonomy(
        tree_file: Path,
        blast_hits: List[BlastHit],
        output_tree: Path,
        label_field: str = "subject_rank_name"
    ) -> bool:
        """
        Relabel tree nodes with taxonomic names.
        
        This reads a Newick tree file and replaces sequence IDs with taxonomic names
        from the BLAST hits.
        
        Args:
            tree_file: Path to input tree file (Newick format)
            blast_hits: List of BlastHit objects with taxonomic information
            output_tree: Path to output tree file with relabeled nodes
            label_field: Field from BlastHit to use for labels (default: "subject_rank_name")
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Create mapping of sequence accession to taxonomic name
            # Trees will have accessions (e.g., MZ387488.1) not full IDs
            label_map = {}
            for hit in blast_hits:
                label = getattr(hit, label_field, None)
                if label:
                    # Clean up the label for tree format (Newick format constraints)
                    # Replace spaces, colons, parentheses, commas, and semicolons
                    clean_label = (label.replace(' ', '_')
                                  .replace(':', '_')
                                  .replace('(', '_')
                                  .replace(')', '_')
                                  .replace(',', '_')
                                  .replace(';', '_'))
                    # Use accession for mapping since that's what appears in trees
                    accession = hit.get_accession()
                    label_map[accession] = clean_label
                    # Also handle potential header prefixes from trimmed files
                    if hit.subject_rank_name:
                        header_with_tax = f"{accession} {hit.subject_rank_name}"
                        label_map[header_with_tax] = clean_label
            
            # Read the tree file
            with open(tree_file, 'r') as f:
                tree_string = f.read()
            
            # Replace sequence IDs with taxonomic names
            for seq_id, tax_name in label_map.items():
                # Handle different possible formats in the tree
                tree_string = tree_string.replace(f"({seq_id}:", f"({tax_name}:")
                tree_string = tree_string.replace(f",{seq_id}:", f",{tax_name}:")
                tree_string = tree_string.replace(f" {seq_id}:", f" {tax_name}:")
            
            # Write the relabeled tree
            with open(output_tree, 'w') as f:
                f.write(tree_string)
            
            logger.info(f"Tree relabeled with {len(label_map)} taxonomic names")
            logger.info(f"Relabeled tree saved to: {output_tree}")
            return True
            
        except Exception as e:
            logger.error(f"Error relabeling tree: {e}")
            return False


def process_query_alignment_and_tree(
    query_id: str,
    query_dir: Path,
    blast_hits: List[BlastHit],
    query_fasta: Path,
    num_threads: int = 1
) -> Dict[str, Optional[Path]]:
    """
    Complete pipeline for a single query: trim, align, and build tree.
    
    Args:
        query_id: Query sequence identifier
        query_dir: Directory containing query-specific files
        blast_hits: List of BlastHit objects for this query (with taxonomy info)
        query_fasta: Path to original query FASTA file
        num_threads: Number of threads to use
        
    Returns:
        Dictionary with paths to generated files:
            - 'trimmed_fasta': Trimmed sequences
            - 'alignment': Aligned sequences
            - 'tree': Phylogenetic tree
            - 'labeled_tree': Tree with taxonomic labels
    """
    results = {
        'trimmed_fasta': None,
        'alignment': None,
        'tree': None,
        'labeled_tree': None
    }
    
    # File paths
    safe_query_id = query_id.replace('|', '_').replace('/', '_')
    representatives_fasta = query_dir / f"{safe_query_id}_representatives.fasta"
    trimmed_fasta = query_dir / f"{safe_query_id}_trimmed.fasta"
    alignment_fasta = query_dir / f"{safe_query_id}_aligned.fasta"
    tree_prefix = query_dir / f"{safe_query_id}_tree"
    tree_file = Path(str(tree_prefix) + ".treefile")
    labeled_tree = query_dir / f"{safe_query_id}_tree_labeled.treefile"
    
    # Step 1: Read query sequence and add it to the representatives file
    try:
        query_sequences = FastaReader.read_fasta(query_fasta)
        if query_id not in query_sequences:
            logger.error(f"Query {query_id} not found in {query_fasta}")
            return results
        
        # Read representatives and combine with query
        if not representatives_fasta.exists():
            logger.error(f"Representatives file not found: {representatives_fasta}")
            return results
        
        # Create combined FASTA with query + representatives
        combined_fasta = query_dir / f"{safe_query_id}_with_query.fasta"
        with open(combined_fasta, 'w') as out:
            # Write query first
            query_seq = query_sequences[query_id]
            out.write(f">{query_id}\n")
            for i in range(0, len(query_seq), 60):
                out.write(query_seq[i:i+60] + "\n")
            
            # Append representatives
            with open(representatives_fasta, 'r') as rep:
                out.write(rep.read())
        
        logger.info(f"Combined query with representatives: {combined_fasta}")
        
    except Exception as e:
        logger.error(f"Error preparing sequences: {e}")
        return results
    
    # Step 2: Trim sequences based on BLAST coordinates
    logger.info(f"Trimming sequences for {query_id}...")
    if SequenceTrimmer.trim_sequences_from_blast_hits(
        fasta_path=combined_fasta,
        blast_hits=blast_hits,
        output_fasta=trimmed_fasta,
        query_id=query_id
    ):
        results['trimmed_fasta'] = trimmed_fasta
        logger.info(f"Trimmed sequences saved to: {trimmed_fasta}")
    else:
        logger.error(f"Failed to trim sequences for {query_id}")
        return results
    
    # Step 3: Align sequences with MAFFT
    logger.info(f"Aligning sequences for {query_id}...")
    if MAFFTAligner.align_sequences(
        input_fasta=trimmed_fasta,
        output_fasta=alignment_fasta,
        auto_orient=True,
        num_threads=num_threads
    ):
        results['alignment'] = alignment_fasta
        logger.info(f"Alignment saved to: {alignment_fasta}")
    else:
        logger.error(f"Failed to align sequences for {query_id}")
        return results
    
    # Step 4: Build phylogenetic tree with IQTree
    logger.info(f"Building phylogenetic tree for {query_id}...")
    if IQTreeBuilder.build_tree(
        alignment_fasta=alignment_fasta,
        output_prefix=tree_prefix,
        num_threads=num_threads
    ):
        results['tree'] = tree_file
        logger.info(f"Tree saved to: {tree_file}")
    else:
        logger.error(f"Failed to build tree for {query_id}")
        return results
    
    # Step 5: Relabel tree with taxonomic names
    logger.info(f"Relabeling tree with taxonomic names for {query_id}...")
    if IQTreeBuilder.relabel_tree_with_taxonomy(
        tree_file=tree_file,
        blast_hits=blast_hits,
        output_tree=labeled_tree
    ):
        results['labeled_tree'] = labeled_tree
        logger.info(f"Labeled tree saved to: {labeled_tree}")
    else:
        logger.warning(f"Failed to relabel tree for {query_id}, but unlabeled tree is available")
    
    return results
