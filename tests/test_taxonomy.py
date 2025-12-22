"""
Tests for taxonomy extraction functionality.
"""

import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
import pytest

from eplace_lib.taxonomy import (
    TaxonomicInfo,
    TaxonomyExtractor,
    SequenceExtractor,
    process_blast_results_for_taxonomy
)
from eplace_lib.blast_analysis import BlastHit


class TestTaxonomyExtractor:
    """Test cases for TaxonomyExtractor class."""
    
    def test_parse_sequence_id_simple(self):
        """Test parsing a simple sequence ID."""
        seq_id = "NC_001234.5"
        tax_info = TaxonomyExtractor.parse_sequence_id(seq_id)
        
        assert tax_info.sequence_id == seq_id
    
    def test_parse_sequence_id_with_gi(self):
        """Test parsing sequence ID with GI number."""
        seq_id = "gi|123456|ref|NC_001234.5|"
        tax_info = TaxonomyExtractor.parse_sequence_id(seq_id)
        
        assert tax_info.sequence_id == seq_id
        assert tax_info.taxid == "123456"
    
    def test_extract_taxonomy_from_hits(self):
        """Test extracting taxonomy from BLAST hits."""
        hits = [
            BlastHit(
                query_id='seq1', subject_id='gi|123|ref|NC_001',
                percent_identity=95.0, alignment_length=200,
                query_length=500, subject_length=1000,
                query_start=1, query_end=450,
                subject_start=100, subject_end=549,
                evalue=1e-50, bit_score=250.0,
                query_coverage=90.0
            ),
            BlastHit(
                query_id='seq1', subject_id='gi|456|ref|NC_002',
                percent_identity=92.0, alignment_length=180,
                query_length=500, subject_length=900,
                query_start=1, query_end=400,
                subject_start=50, subject_end=449,
                evalue=1e-45, bit_score=240.0,
                query_coverage=80.0
            ),
        ]
        
        extractor = TaxonomyExtractor()
        taxonomy_info = extractor.extract_taxonomy_from_hits(hits)
        
        assert len(taxonomy_info) == 2
        assert 'gi|123|ref|NC_001' in taxonomy_info
        assert 'gi|456|ref|NC_002' in taxonomy_info
    
    def test_group_hits_by_query(self):
        """Test grouping BLAST hits by query."""
        hits = [
            BlastHit(
                query_id='seq1', subject_id='subj1',
                percent_identity=95.0, alignment_length=200,
                query_length=500, subject_length=1000,
                query_start=1, query_end=450,
                subject_start=100, subject_end=549,
                evalue=1e-50, bit_score=250.0,
                query_coverage=90.0
            ),
            BlastHit(
                query_id='seq1', subject_id='subj2',
                percent_identity=92.0, alignment_length=180,
                query_length=500, subject_length=900,
                query_start=1, query_end=400,
                subject_start=50, subject_end=449,
                evalue=1e-45, bit_score=240.0,
                query_coverage=80.0
            ),
            BlastHit(
                query_id='seq2', subject_id='subj3',
                percent_identity=88.0, alignment_length=150,
                query_length=400, subject_length=800,
                query_start=1, query_end=350,
                subject_start=200, subject_end=549,
                evalue=1e-40, bit_score=230.0,
                query_coverage=87.5
            ),
        ]
        
        extractor = TaxonomyExtractor()
        grouped = extractor.group_hits_by_query(hits)
        
        assert len(grouped) == 2
        assert len(grouped['seq1']) == 2
        assert len(grouped['seq2']) == 1
    
    def test_select_representatives_by_rank_valid(self):
        """Test selecting representatives by valid rank."""
        hits = [
            BlastHit(
                query_id='seq1', subject_id='gi|123|ref|NC_001',
                percent_identity=95.0, alignment_length=200,
                query_length=500, subject_length=1000,
                query_start=1, query_end=450,
                subject_start=100, subject_end=549,
                evalue=1e-50, bit_score=250.0,
                query_coverage=90.0
            ),
            BlastHit(
                query_id='seq1', subject_id='gi|456|ref|NC_002',
                percent_identity=92.0, alignment_length=180,
                query_length=500, subject_length=900,
                query_start=1, query_end=400,
                subject_start=50, subject_end=449,
                evalue=1e-45, bit_score=240.0,
                query_coverage=80.0
            ),
            BlastHit(
                query_id='seq1', subject_id='gi|789|ref|NC_003',
                percent_identity=90.0, alignment_length=170,
                query_length=500, subject_length=850,
                query_start=1, query_end=420,
                subject_start=60, subject_end=479,
                evalue=1e-40, bit_score=230.0,
                query_coverage=84.0
            ),
        ]
        
        extractor = TaxonomyExtractor()
        representatives = extractor.select_representatives_by_rank(
            hits=hits,
            rank='species'
        )
        
        # Should return at least one representative
        assert len(representatives) >= 1
        # Best hit should be included (highest bit score)
        assert representatives[0].bit_score == 250.0
    
    def test_select_representatives_by_rank_invalid(self):
        """Test selecting representatives with invalid rank."""
        hits = [
            BlastHit(
                query_id='seq1', subject_id='subj1',
                percent_identity=95.0, alignment_length=200,
                query_length=500, subject_length=1000,
                query_start=1, query_end=450,
                subject_start=100, subject_end=549,
                evalue=1e-50, bit_score=250.0,
                query_coverage=90.0
            ),
        ]
        
        extractor = TaxonomyExtractor()
        with pytest.raises(ValueError, match="Invalid rank"):
            extractor.select_representatives_by_rank(
                hits=hits,
                rank='invalid_rank'
            )
    
    def test_select_representatives_empty_hits(self):
        """Test selecting representatives from empty hits list."""
        extractor = TaxonomyExtractor()
        representatives = extractor.select_representatives_by_rank(
            hits=[],
            rank='species'
        )
        
        assert len(representatives) == 0


class TestSequenceExtractor:
    """Test cases for SequenceExtractor class."""
    
    def test_init_default(self):
        """Test SequenceExtractor initialization with defaults."""
        with patch.dict('os.environ', {}, clear=True):
            extractor = SequenceExtractor()
            assert extractor.blastdb_path == Path.home() / "blastdb"
    
    def test_init_with_env_var(self):
        """Test SequenceExtractor initialization with BLASTDB env var."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.dict('os.environ', {'BLASTDB': tmpdir}):
                extractor = SequenceExtractor()
                assert extractor.blastdb_path == Path(tmpdir)
    
    def test_init_with_path(self):
        """Test SequenceExtractor initialization with explicit path."""
        test_path = Path("/test/path")
        extractor = SequenceExtractor(blastdb_path=test_path)
        assert extractor.blastdb_path == test_path
    
    @patch('subprocess.run')
    def test_check_blastdbcmd_available_true(self, mock_run):
        """Test checking blastdbcmd availability when available."""
        mock_run.return_value = MagicMock(returncode=0)
        
        extractor = SequenceExtractor()
        assert extractor.check_blastdbcmd_available() is True
    
    @patch('subprocess.run')
    def test_check_blastdbcmd_available_false(self, mock_run):
        """Test checking blastdbcmd availability when not available."""
        mock_run.side_effect = FileNotFoundError()
        
        extractor = SequenceExtractor()
        assert extractor.check_blastdbcmd_available() is False
    
    @patch('eplace_lib.taxonomy.SequenceExtractor.check_blastdbcmd_available')
    @patch('subprocess.run')
    def test_extract_sequences_success(self, mock_run, mock_check):
        """Test extracting sequences successfully."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_fasta = tmppath / "output.fasta"
            
            mock_check.return_value = True
            mock_run.return_value = MagicMock(returncode=0, stderr="")
            
            extractor = SequenceExtractor(blastdb_path=tmppath)
            success = extractor.extract_sequences(
                sequence_ids=['seq1', 'seq2'],
                output_fasta=output_fasta
            )
            
            assert success is True
            mock_run.assert_called_once()
    
    @patch('eplace_lib.taxonomy.SequenceExtractor.check_blastdbcmd_available')
    def test_extract_sequences_not_available(self, mock_check):
        """Test extracting sequences when blastdbcmd not available."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_fasta = tmppath / "output.fasta"
            
            mock_check.return_value = False
            
            extractor = SequenceExtractor()
            with pytest.raises(RuntimeError, match="blastdbcmd is not available"):
                extractor.extract_sequences(
                    sequence_ids=['seq1', 'seq2'],
                    output_fasta=output_fasta
                )
    
    @patch('eplace_lib.taxonomy.SequenceExtractor.check_blastdbcmd_available')
    def test_extract_sequences_empty_list(self, mock_check):
        """Test extracting sequences with empty ID list."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_fasta = tmppath / "output.fasta"
            
            mock_check.return_value = True
            
            extractor = SequenceExtractor()
            success = extractor.extract_sequences(
                sequence_ids=[],
                output_fasta=output_fasta
            )
            
            assert success is False
    
    @patch('eplace_lib.taxonomy.SequenceExtractor.extract_sequences')
    def test_extract_representatives_for_query(self, mock_extract):
        """Test extracting representatives for a single query."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            hits = [
                BlastHit(
                    query_id='seq1', subject_id='subj1',
                    percent_identity=95.0, alignment_length=200,
                    query_length=500, subject_length=1000,
                    query_start=1, query_end=450,
                    subject_start=100, subject_end=549,
                    evalue=1e-50, bit_score=250.0,
                    query_coverage=90.0
                ),
            ]
            
            mock_extract.return_value = True
            
            extractor = SequenceExtractor()
            output_fasta = extractor.extract_representatives_for_query(
                query_id='seq1',
                representative_hits=hits,
                output_dir=tmppath
            )
            
            assert output_fasta is not None
            assert output_fasta.parent.exists()
    
    def test_extract_representatives_for_query_empty_hits(self):
        """Test extracting representatives with no hits."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            extractor = SequenceExtractor()
            output_fasta = extractor.extract_representatives_for_query(
                query_id='seq1',
                representative_hits=[],
                output_dir=tmppath
            )
            
            assert output_fasta is None


class TestProcessBlastResultsForTaxonomy:
    """Test cases for process_blast_results_for_taxonomy function."""
    
    @patch('eplace_lib.taxonomy.SequenceExtractor.extract_representatives_for_query')
    def test_process_blast_results_for_taxonomy(self, mock_extract):
        """Test processing BLAST results for taxonomy."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            hits = [
                BlastHit(
                    query_id='seq1', subject_id='subj1',
                    percent_identity=95.0, alignment_length=200,
                    query_length=500, subject_length=1000,
                    query_start=1, query_end=450,
                    subject_start=100, subject_end=549,
                    evalue=1e-50, bit_score=250.0,
                    query_coverage=90.0
                ),
                BlastHit(
                    query_id='seq2', subject_id='subj2',
                    percent_identity=92.0, alignment_length=180,
                    query_length=500, subject_length=900,
                    query_start=1, query_end=400,
                    subject_start=50, subject_end=449,
                    evalue=1e-45, bit_score=240.0,
                    query_coverage=80.0
                ),
            ]
            
            mock_extract.return_value = tmppath / "output.fasta"
            
            results = process_blast_results_for_taxonomy(
                blast_hits=hits,
                output_dir=tmppath,
                rank='species'
            )
            
            assert len(results) == 2
            assert 'seq1' in results
            assert 'seq2' in results
    
    @patch('eplace_lib.taxonomy.SequenceExtractor.extract_representatives_for_query')
    def test_process_blast_results_different_ranks(self, mock_extract):
        """Test processing BLAST results with different taxonomic ranks."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            hits = [
                BlastHit(
                    query_id='seq1', subject_id='subj1',
                    percent_identity=95.0, alignment_length=200,
                    query_length=500, subject_length=1000,
                    query_start=1, query_end=450,
                    subject_start=100, subject_end=549,
                    evalue=1e-50, bit_score=250.0,
                    query_coverage=90.0
                ),
            ]
            
            mock_extract.return_value = tmppath / "output.fasta"
            
            for rank in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
                results = process_blast_results_for_taxonomy(
                    blast_hits=hits,
                    output_dir=tmppath,
                    rank=rank
                )
                
                assert len(results) == 1
                assert 'seq1' in results
