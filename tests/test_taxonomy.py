"""
Tests for taxonomy extraction functionality.
"""

import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
import pytest

from eplace_lib.taxonomy import (
    TaxonomyExtractor,
    SequenceExtractor,
    process_blast_results_for_taxonomy
)
from eplace_lib.blast_analysis import BlastHit


class TestTaxonomyExtractor:
    """Test cases for TaxonomyExtractor class."""
    
    def setup_method(self):
        self.taxonomy_extractor = TaxonomyExtractor("genus")
        self.hits = [
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

    def test_parse_taxids(self):
        """Test parsing a simple sequence ID."""
        seq_id = "NC_001234.5"
        tax_info = self.taxonomy_extractor.parse_taxids([9606, 590])

        homo_tuple = ('9605', 'Homo')
        salm_tuple = ('590', 'Salmonella')
        assert tax_info["9606"] == homo_tuple
        assert tax_info["590"] == salm_tuple
    
    def test_group_hits_by_query(self):
        """Test grouping BLAST hits by query."""
        hits = self.hits
        
        extractor = self.taxonomy_extractor
        grouped = extractor.group_hits_by_query(hits)
        
        assert len(grouped) == 3
        assert len(grouped['seq1']) == 1
        assert len(grouped['seq2']) == 1
    
    def test_select_representatives_by_rank_valid(self):
        """Test selecting representatives by valid rank."""
        hits = self.hits
        
        extractor = self.taxonomy_extractor
        representatives = extractor.select_representatives_by_rank(
            hits=hits
        )
        
        # Should return at least one representative
        assert len(representatives) >= 1
        # Best hit should be included (highest bit score)
        assert representatives[0].bit_score == 998
    
    def test_select_representatives_empty_hits(self):
        """Test selecting representatives from empty hits list."""
        extractor = self.taxonomy_extractor
        representatives = extractor.select_representatives_by_rank(
            hits=[],
        )
        
        assert len(representatives) == 0


class TestSequenceExtractor:
    """Test cases for SequenceExtractor class."""
    
    def setup_method(self):
        self.hits = [
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
            
            hits = self.hits
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
            
            mock_extract.return_value = tmppath / "output.fasta"
            
            results = process_blast_results_for_taxonomy(
                blast_hits=hits,
                output_dir=tmppath,
                rank='genus'
            )
            
            assert len(results) == 3
            assert 'seq1' in results
            assert 'seq2' in results
    
