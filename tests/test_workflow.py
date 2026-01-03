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
CAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTGGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCAGATGTGCCCAGATGGGATTAGCTTGTTGGTGAGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGTGTTGTGGTTAATAACCGCAGCAATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGGCACG
>seq2
AGAAGGCTTTGGCTTCTGATAGTCATGGACTCACTAGGCTGCTGAGGAAGATCAATAATACCTACTGGAATCAGTCATGAGAAGTCAAGCATGGAAATTGTGAATTGTGTGTGTGGCCAGACCAGTACCTCCAAGTGTTCAGAAGATGTGTGACCAGACAAAACACAGTAAATGCTGCCCAGCAAAAGGCAATCAATGCTGCCCACCACAGCAGAACCAGTGCTGCCAGTCAAAAGGCAATCAATGCTGCCCACCAAAACAGAACCAGTGCTGCCAGCCAAAAGGCAGTCAATGCTGCCCACCAAAACACAATCACTGCTGCCAGCCAAAACCCCCATGCTGCATTCAGGCCAGGTGCTGTGGTTTGGAGACCAAGCCTGAAGTCTCACCCCTTAACATGGAGTCTGAGCCCAACTCACC
"""
            query_fasta.write_text(query_content)
            
            # Mock BLAST results
            mock_run.return_value = True
            mock_hits = [
                BlastHit(
                    query_id='seq1', subject_id='gi|156763568|gb|EU014687.1|',
                    percent_identity=100.000, 
                    alignment_length=540,
                    query_length=540,
                    subject_length=1432,
                    query_start=1, 
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=998,
                    query_coverage=100,
                    subject_taxid="149539",
                    subject_taxids="149539",
                    subject_rank_tid="590",
                    subject_rank_name="Salmonella"
                ),
                BlastHit(
                    query_id='seq2', subject_id='gi|34190046|gb|BC014593.2|',
                    percent_identity=100.000, 
                    alignment_length=420,
                    query_length=420,
                    subject_length=784,
                    query_start=1, 
                    query_end=420,
                    subject_start=1,
                    subject_end=420,
                    evalue=0.0,
                    bit_score=776,
                    query_coverage=100,
                    subject_taxid="9606",
                    subject_taxids="9606",
                    subject_rank_tid="9605",
                    subject_rank_name="Homo"
                ),
                BlastHit(
                    query_id='seq3', subject_id='gi|2694387494|ref|XM_055113774.3|',
                    percent_identity=91.304, 
                    alignment_length=115,
                    query_length=420,
                    subject_length=1626,
                    query_start=218, 
                    query_end=331,
                    subject_start=173,
                    subject_end=286,
                    evalue=3.56e-33,
                    bit_score=156,
                    query_coverage=27.3809523809524,
                    subject_taxid="9597",
                    subject_taxids="9597",
                    subject_rank_tid="9596",
                    subject_rank_name="Pan"
                )
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
            assert len(filtered_hits) == 2
            
            # Step 2: Process results for taxonomy
            output_dir = tmppath / "output"
            results = process_blast_results_for_taxonomy(
                blast_hits=filtered_hits,
                output_dir=output_dir,
                rank='genus'
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
                    query_id='seq1', subject_id='gi|156763568|gb|EU014687.1|',
                    percent_identity=100.000, 
                    alignment_length=540,
                    query_length=540,
                    subject_length=1432,
                    query_start=1, 
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=998,
                    query_coverage=100,
                    subject_taxid="149539",
                    subject_taxids="149539",
                    subject_rank_tid="590",
                    subject_rank_name="Salmonella"
                ),
                BlastHit(
                    query_id='seq2', subject_id='gi|34190046|gb|BC014593.2|',
                    percent_identity=100.000, 
                    alignment_length=420,
                    query_length=420,
                    subject_length=784,
                    query_start=1, 
                    query_end=420,
                    subject_start=1,
                    subject_end=420,
                    evalue=0.0,
                    bit_score=776,
                    query_coverage=100,
                    subject_taxid="9606",
                    subject_taxids="9606",
                    subject_rank_tid="9605",
                    subject_rank_name="Homo"
                ),
                BlastHit(
                    query_id='seq3', subject_id='gi|2694387494|ref|XM_055113774.3|',
                    percent_identity=91.304, 
                    alignment_length=115,
                    query_length=420,
                    subject_length=1626,
                    query_start=218, 
                    query_end=331,
                    subject_start=173,
                    subject_end=286,
                    evalue=3.56e-33,
                    bit_score=156,
                    query_coverage=27.3809523809524,
                    subject_taxid="9597",
                    subject_taxids="9597",
                    subject_rank_tid="9596",
                    subject_rank_name="Pan"
                )
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
            for rank in ['phylum', 'class', 'order', 'family', 'genus']:
                output_dir = tmppath / f"output_{rank}"
                results = process_blast_results_for_taxonomy(
                    blast_hits=filtered_hits,
                    output_dir=output_dir,
                    rank=rank
                )
                
                assert len(results) == 2
                assert 'seq1' in results
