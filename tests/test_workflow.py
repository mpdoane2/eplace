"""
Integration tests for the complete BLAST workflow.
"""

import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
import pytest
from collections import defaultdict

from eplace_lib.blast_analysis import run_blast_search, FastaReader, BlastHit
from eplace_lib.taxonomy import process_blast_results_for_taxonomy


class TestBlastWorkflow:
    """Test cases for the complete BLAST workflow."""
    
    @patch('eplace_lib.blast_analysis.BlastRunner.run_blastn')
    @patch('eplace_lib.blast_analysis.BlastRunner.parse_blast_results')
    @patch('eplace_lib.taxonomy.SequenceExtractor.extract_sequences')
    def test_complete_workflow(self, mock_extract, mock_parse, mock_run):
        """Test the complete workflow from FASTA to representative sequences."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create test query FASTA
            query_fasta = tmppath / "query.fasta"
            query_content = """>seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
"""
            query_fasta.write_text(query_content)
            
            # Mock BLAST results
            mock_run.return_value = True
            mock_hits = [
                BlastHit(
                    query_id='seq1', subject_id='gi|123|ref|NC_001',
                    percent_identity=95.0, alignment_length=40,
                    query_length=40, subject_length=1000,
                    query_start=1, query_end=40,
                    subject_start=100, subject_end=139,
                    evalue=1e-50, bit_score=250.0,
                    query_coverage=100.0
                ),
                BlastHit(
                    query_id='seq1', subject_id='gi|456|ref|NC_002',
                    percent_identity=92.0, alignment_length=38,
                    query_length=40, subject_length=900,
                    query_start=1, query_end=38,
                    subject_start=50, subject_end=87,
                    evalue=1e-45, bit_score=240.0,
                    query_coverage=95.0
                ),
                BlastHit(
                    query_id='seq2', subject_id='gi|789|ref|NC_003',
                    percent_identity=90.0, alignment_length=40,
                    query_length=40, subject_length=800,
                    query_start=1, query_end=40,
                    subject_start=200, subject_end=239,
                    evalue=1e-40, bit_score=230.0,
                    query_coverage=100.0
                ),
            ]
            mock_parse.return_value = mock_hits
            mock_extract.return_value = True
            
            # Step 1: Run BLAST search
            blast_output = tmppath / "blast_results.txt"
            success, filtered_hits = run_blast_search(
                query_fasta=query_fasta,
                output_file=blast_output,
                min_identity=90.0,
                min_coverage=80.0
            )
            
            assert success is True
            assert len(filtered_hits) == 3
            
            # Step 2: Process results for taxonomy
            output_dir = tmppath / "output"
            results = process_blast_results_for_taxonomy(
                blast_hits=filtered_hits,
                output_dir=output_dir,
                rank='species'
            )
            
            assert len(results) == 2
            assert 'seq1' in results
            assert 'seq2' in results
    
    def test_read_fasta_workflow(self):
        """Test reading a FASTA file as part of workflow."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create test FASTA
            query_fasta = tmppath / "query.fasta"
            query_content = """>seq1 Test sequence 1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2 Test sequence 2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>seq3 Test sequence 3
TTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGG
"""
            query_fasta.write_text(query_content)
            
            # Read sequences
            sequences = FastaReader.read_fasta(query_fasta)
            
            assert len(sequences) == 3
            assert 'seq1' in sequences
            assert 'seq2' in sequences
            assert 'seq3' in sequences
            assert len(sequences['seq1']) == 40
            assert len(sequences['seq2']) == 40
            assert len(sequences['seq3']) == 40
    
    @patch('eplace_lib.blast_analysis.BlastRunner.run_blastn')
    @patch('eplace_lib.blast_analysis.BlastRunner.parse_blast_results')
    def test_workflow_with_no_hits(self, mock_parse, mock_run):
        """Test workflow when BLAST returns no hits."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create test query FASTA
            query_fasta = tmppath / "query.fasta"
            query_content = """>seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""
            query_fasta.write_text(query_content)
            
            # Mock BLAST results with no hits
            mock_run.return_value = True
            mock_parse.return_value = []
            
            # Run BLAST search
            blast_output = tmppath / "blast_results.txt"
            success, filtered_hits = run_blast_search(
                query_fasta=query_fasta,
                output_file=blast_output
            )
            
            assert success is True
            assert len(filtered_hits) == 0
    
    @patch('eplace_lib.blast_analysis.BlastRunner.run_blastn')
    @patch('eplace_lib.blast_analysis.BlastRunner.parse_blast_results')
    @patch('eplace_lib.taxonomy.SequenceExtractor.extract_sequences')
    def test_workflow_different_ranks(self, mock_extract, mock_parse, mock_run):
        """Test workflow with different taxonomic ranks."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create test query FASTA
            query_fasta = tmppath / "query.fasta"
            query_content = """>seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""
            query_fasta.write_text(query_content)
            
            # Mock BLAST results
            mock_run.return_value = True
            mock_hits = [
                BlastHit(
                    query_id='seq1', subject_id='gi|123|ref|NC_001',
                    percent_identity=95.0, alignment_length=40,
                    query_length=40, subject_length=1000,
                    query_start=1, query_end=40,
                    subject_start=100, subject_end=139,
                    evalue=1e-50, bit_score=250.0,
                    query_coverage=100.0
                ),
            ]
            mock_parse.return_value = mock_hits
            mock_extract.return_value = True
            
            # Run BLAST search
            blast_output = tmppath / "blast_results.txt"
            success, filtered_hits = run_blast_search(
                query_fasta=query_fasta,
                output_file=blast_output
            )
            
            assert success is True
            
            # Test with different ranks
            for rank in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
                output_dir = tmppath / f"output_{rank}"
                results = process_blast_results_for_taxonomy(
                    blast_hits=filtered_hits,
                    output_dir=output_dir,
                    rank=rank
                )
                
                assert len(results) == 1
                assert 'seq1' in results
