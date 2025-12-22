"""
Tests for BLAST analysis functionality.
"""

import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
import pytest

from eplace_lib.blast_analysis import (
    FastaReader,
    BlastRunner,
    BlastHit,
    run_blast_search
)


class TestFastaReader:
    """Test cases for FastaReader class."""
    
    def test_read_fasta_simple(self):
        """Test reading a simple FASTA file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            fasta_file = tmppath / "test.fasta"
            
            # Create test FASTA file
            fasta_content = """>seq1
ATCGATCG
>seq2
GCTAGCTA
"""
            fasta_file.write_text(fasta_content)
            
            # Read FASTA
            sequences = FastaReader.read_fasta(fasta_file)
            
            assert len(sequences) == 2
            assert sequences['seq1'] == 'ATCGATCG'
            assert sequences['seq2'] == 'GCTAGCTA'
    
    def test_read_fasta_multiline(self):
        """Test reading FASTA with multiline sequences."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            fasta_file = tmppath / "test.fasta"
            
            fasta_content = """>seq1 description here
ATCGATCG
GCTAGCTA
TTAACCGG
>seq2
AAAA
TTTT
"""
            fasta_file.write_text(fasta_content)
            
            sequences = FastaReader.read_fasta(fasta_file)
            
            assert len(sequences) == 2
            assert sequences['seq1'] == 'ATCGATCGGCTAGCTATTAACCGG'
            assert sequences['seq2'] == 'AAAATTTT'
    
    def test_read_fasta_nonexistent(self):
        """Test reading non-existent FASTA file."""
        with pytest.raises(FileNotFoundError):
            FastaReader.read_fasta(Path("/nonexistent/file.fasta"))
    
    def test_read_fasta_empty(self):
        """Test reading empty FASTA file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            fasta_file = tmppath / "empty.fasta"
            fasta_file.write_text("")
            
            with pytest.raises(ValueError, match="No sequences found"):
                FastaReader.read_fasta(fasta_file)
    
    def test_get_sequence_lengths(self):
        """Test getting sequence lengths."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            fasta_file = tmppath / "test.fasta"
            
            fasta_content = """>seq1
ATCGATCG
>seq2
GCTAGCTAGCTA
"""
            fasta_file.write_text(fasta_content)
            
            lengths = FastaReader.get_sequence_lengths(fasta_file)
            
            assert lengths['seq1'] == 8
            assert lengths['seq2'] == 12


class TestBlastRunner:
    """Test cases for BlastRunner class."""
    
    def test_init_default(self):
        """Test BlastRunner initialization with defaults."""
        with patch.dict('os.environ', {}, clear=True):
            runner = BlastRunner()
            assert runner.blastdb_path == Path.home() / "blastdb"
    
    def test_init_with_env_var(self):
        """Test BlastRunner initialization with BLASTDB env var."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.dict('os.environ', {'BLASTDB': tmpdir}):
                runner = BlastRunner()
                assert runner.blastdb_path == Path(tmpdir)
    
    def test_init_with_path(self):
        """Test BlastRunner initialization with explicit path."""
        test_path = Path("/test/path")
        runner = BlastRunner(blastdb_path=test_path)
        assert runner.blastdb_path == test_path
    
    @patch('subprocess.run')
    def test_check_blastn_available_true(self, mock_run):
        """Test checking blastn availability when available."""
        mock_run.return_value = MagicMock(returncode=0)
        
        runner = BlastRunner()
        assert runner.check_blastn_available() is True
    
    @patch('subprocess.run')
    def test_check_blastn_available_false(self, mock_run):
        """Test checking blastn availability when not available."""
        mock_run.side_effect = FileNotFoundError()
        
        runner = BlastRunner()
        assert runner.check_blastn_available() is False
    
    def test_parse_blast_results(self):
        """Test parsing BLAST tabular output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            blast_output = tmppath / "blast_results.txt"
            
            # Create mock BLAST output
            blast_content = """seq1\tgi|123|ref|NC_001\t95.5\t200\t500\t1000\t1\t200\t100\t299\t1e-50\t250.0
seq1\tgi|456|ref|NC_002\t92.0\t180\t500\t900\t1\t180\t50\t229\t1e-45\t240.0
seq2\tgi|789|ref|NC_003\t88.5\t150\t400\t800\t1\t150\t200\t349\t1e-40\t230.0
"""
            blast_output.write_text(blast_content)
            
            runner = BlastRunner()
            hits = runner.parse_blast_results(blast_output)
            
            assert len(hits) == 3
            assert hits[0].query_id == 'seq1'
            assert hits[0].subject_id == 'gi|123|ref|NC_001'
            assert hits[0].percent_identity == 95.5
            assert hits[0].alignment_length == 200
            assert hits[0].query_length == 500
            assert hits[0].evalue == 1e-50
            assert hits[0].bit_score == 250.0
            assert hits[0].query_coverage == pytest.approx(40.0, abs=0.1)
    
    def test_parse_blast_results_invalid_format(self):
        """Test parsing BLAST output with invalid format."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            blast_output = tmppath / "blast_results.txt"
            
            # Create invalid BLAST output (too few fields)
            blast_content = "seq1\tgi|123\t95.5\n"
            blast_output.write_text(blast_content)
            
            runner = BlastRunner()
            with pytest.raises(ValueError, match="Invalid BLAST output format"):
                runner.parse_blast_results(blast_output)
    
    def test_filter_blast_hits(self):
        """Test filtering BLAST hits."""
        hits = [
            BlastHit(
                query_id='seq1', subject_id='subj1',
                percent_identity=95.0, alignment_length=200,
                query_length=500, subject_length=1000,
                query_start=1, query_end=450,  # 90% coverage
                subject_start=100, subject_end=549,
                evalue=1e-50, bit_score=250.0,
                query_coverage=90.0
            ),
            BlastHit(
                query_id='seq1', subject_id='subj2',
                percent_identity=85.0, alignment_length=150,  # Below identity threshold
                query_length=500, subject_length=900,
                query_start=1, query_end=400,
                subject_start=50, subject_end=449,
                evalue=1e-40, bit_score=230.0,
                query_coverage=80.0
            ),
            BlastHit(
                query_id='seq1', subject_id='subj3',
                percent_identity=92.0, alignment_length=100,
                query_length=500, subject_length=800,
                query_start=1, query_end=100,  # Only 20% coverage
                subject_start=200, subject_end=299,
                evalue=1e-30, bit_score=220.0,
                query_coverage=20.0
            ),
        ]
        
        runner = BlastRunner()
        filtered = runner.filter_blast_hits(
            hits,
            min_identity=90.0,
            min_coverage=80.0
        )
        
        assert len(filtered) == 1
        assert filtered[0].subject_id == 'subj1'
    
    @patch('eplace_lib.blast_analysis.BlastRunner.check_blastn_available')
    @patch('subprocess.run')
    def test_run_blastn_success(self, mock_run, mock_check):
        """Test running blastn successfully."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            query_file = tmppath / "query.fasta"
            output_file = tmppath / "output.txt"
            
            # Create query file
            query_file.write_text(">seq1\nATCG\n")
            
            # Mock blastn availability and execution
            mock_check.return_value = True
            mock_run.return_value = MagicMock(returncode=0, stderr="")
            
            runner = BlastRunner(blastdb_path=tmppath)
            success = runner.run_blastn(query_file, output_file)
            
            assert success is True
            mock_run.assert_called_once()
    
    def test_run_blastn_no_query_file(self):
        """Test running blastn with non-existent query file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            query_file = tmppath / "nonexistent.fasta"
            output_file = tmppath / "output.txt"
            
            runner = BlastRunner()
            with pytest.raises(FileNotFoundError):
                runner.run_blastn(query_file, output_file)
    
    @patch('eplace_lib.blast_analysis.BlastRunner.check_blastn_available')
    def test_run_blastn_not_available(self, mock_check):
        """Test running blastn when not available."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            query_file = tmppath / "query.fasta"
            output_file = tmppath / "output.txt"
            
            query_file.write_text(">seq1\nATCG\n")
            mock_check.return_value = False
            
            runner = BlastRunner()
            with pytest.raises(RuntimeError, match="blastn is not available"):
                runner.run_blastn(query_file, output_file)


class TestRunBlastSearch:
    """Test cases for run_blast_search convenience function."""
    
    @patch('eplace_lib.blast_analysis.BlastRunner.run_blastn')
    @patch('eplace_lib.blast_analysis.BlastRunner.parse_blast_results')
    @patch('eplace_lib.blast_analysis.BlastRunner.filter_blast_hits')
    def test_run_blast_search_success(self, mock_filter, mock_parse, mock_run):
        """Test run_blast_search convenience function."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            query_file = tmppath / "query.fasta"
            output_file = tmppath / "output.txt"
            
            query_file.write_text(">seq1\nATCG\n")
            
            # Mock successful execution
            mock_run.return_value = True
            mock_parse.return_value = [MagicMock()]
            mock_filter.return_value = [MagicMock()]
            
            success, hits = run_blast_search(
                query_fasta=query_file,
                output_file=output_file
            )
            
            assert success is True
            assert len(hits) == 1
    
    @patch('eplace_lib.blast_analysis.BlastRunner.run_blastn')
    def test_run_blast_search_blast_fails(self, mock_run):
        """Test run_blast_search when BLAST fails."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            query_file = tmppath / "query.fasta"
            output_file = tmppath / "output.txt"
            
            query_file.write_text(">seq1\nATCG\n")
            
            # Mock failed execution
            mock_run.return_value = False
            
            success, hits = run_blast_search(
                query_fasta=query_file,
                output_file=output_file
            )
            
            assert success is False
            assert len(hits) == 0
