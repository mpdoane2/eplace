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
    run_blast_search,
    MMseqs2Runner,
    run_mmseqs_search,
    normalize_sequence_id,
    _extract_accession_from_pipe_id,
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
    
    def setup_method(self):
        self.salmonella_taxonomy = {
            'phylum': ('1224', 'Pseudomonadota'),
            'class': ('1236', 'Gammaproteobacteria'),
            'order': ('91347', 'Enterobacterales'),
            'family': ('543', 'Enterobacteriaceae'),
            'genus': ('590', 'Salmonella')
        }
        self.human_taxonomy = {
            'phylum': ('7711', 'Chordata'),
            'class': ('40674', 'Mammalia'),
            'order': ('9443', 'Primates'),
            'family': ('9604', 'Hominidae'),
            'genus': ('9605', 'Homo'),
            'species': ('9606', 'Homo sapiens')
        }
        self.pan_taxonomy = {
            'phylum': ('7711', 'Chordata'),
            'class': ('40674', 'Mammalia'),
            'order': ('9443', 'Primates'),
            'family': ('9604', 'Hominidae'),
            'genus': ('9596', 'Pan'),
            'species': ('9597', 'Pan paniscus')
        }
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
                subject_taxid="590",
                subject_taxids="590",
                subject_taxonomy=self.salmonella_taxonomy
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
                subject_taxids="9606;9605",
                subject_taxonomy=self.human_taxonomy
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
                subject_taxids="9597;9596;9604",
                subject_taxonomy=self.pan_taxonomy
            )
        ]

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
            blast_content = """seq1\tgi|123|ref|NC_001\t95.5\t200\t500\t1000\t1\t200\t100\t299\t1e-50\t250.0\t9606\t7711;40674;9443;9604;9605;9606
seq1\tgi|456|ref|NC_002\t92.0\t180\t500\t900\t1\t180\t50\t229\t1e-45\t240.0\t590\t1224;1236;91347;543;590
seq2\tgi|789|ref|NC_003\t88.5\t150\t400\t800\t1\t150\t200\t349\t1e-40\t230.0\t590\t1224;1236;91347;543;590
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
            assert hits[0].subject_taxid == "9606"
            assert hits[0].subject_taxids == "7711;40674;9443;9604;9605;9606"
    
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
        hits = self.hits
        
        runner = BlastRunner()
        filtered = runner.filter_blast_hits(
            hits,
            min_identity=90.0,
            min_coverage=80.0
        )
        
        assert len(filtered) == 2
        assert filtered[0].subject_id == 'gi|156763568|gb|EU014687.1|'
    
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


class TestBlastHit:
    """Test cases for BlastHit class."""
    
    def test_get_accession_gi_format(self):
        """Test extracting accession from gi|...|gb|...| format."""
        hit = BlastHit(
            query_id="test",
            subject_id="gi|2273658778|gb|MZ387488.1|",
            percent_identity=100.0,
            alignment_length=100,
            query_length=100,
            subject_length=100,
            query_start=1,
            query_end=100,
            subject_start=1,
            subject_end=100,
            evalue=0.0,
            bit_score=100,
            query_coverage=100.0,
            subject_taxid="12345",
            subject_taxids="12345"
        )
        assert hit.get_accession() == "MZ387488.1"
    
    def test_get_accession_ref_format(self):
        """Test extracting accession from ref|...| format."""
        hit = BlastHit(
            query_id="test",
            subject_id="ref|NZ_CP123456.1|",
            percent_identity=100.0,
            alignment_length=100,
            query_length=100,
            subject_length=100,
            query_start=1,
            query_end=100,
            subject_start=1,
            subject_end=100,
            evalue=0.0,
            bit_score=100,
            query_coverage=100.0,
            subject_taxid="12345",
            subject_taxids="12345"
        )
        assert hit.get_accession() == "NZ_CP123456.1"
    
    def test_get_accession_simple_format(self):
        """Test accession when already in simple format."""
        hit = BlastHit(
            query_id="test",
            subject_id="MZ387488.1",
            percent_identity=100.0,
            alignment_length=100,
            query_length=100,
            subject_length=100,
            query_start=1,
            query_end=100,
            subject_start=1,
            subject_end=100,
            evalue=0.0,
            bit_score=100,
            query_coverage=100.0,
            subject_taxid="12345",
            subject_taxids="12345"
        )
        assert hit.get_accession() == "MZ387488.1"
    
    def test_get_accession_gnl_format(self):
        """Test extracting identifier from gnl|database|identifier format."""
        hit = BlastHit(
            query_id="test",
            subject_id="gnl|BL_ORD_ID|12345",
            percent_identity=100.0,
            alignment_length=100,
            query_length=100,
            subject_length=100,
            query_start=1,
            query_end=100,
            subject_start=1,
            subject_end=100,
            evalue=0.0,
            bit_score=100,
            query_coverage=100.0,
            subject_taxid="12345",
            subject_taxids="12345"
        )
        assert hit.get_accession() == "12345"


class TestNormalizeSequenceId:
    """Test cases for normalize_sequence_id helper."""

    def test_gi_gb_format(self):
        """gi|...|gb|ACC| normalizes to ACC."""
        assert normalize_sequence_id("gi|336317909|gb|HQ641676.1|") == "HQ641676.1"

    def test_ref_format(self):
        """ref|ACC| normalizes to ACC."""
        assert normalize_sequence_id("ref|NZ_CP123456.1|") == "NZ_CP123456.1"

    def test_gb_format(self):
        """gb|ACC| normalizes to ACC."""
        assert normalize_sequence_id("gb|MZ387488.1|") == "MZ387488.1"

    def test_plain_accession(self):
        """A plain accession is returned unchanged."""
        assert normalize_sequence_id("HQ641676.1") == "HQ641676.1"

    def test_fasta_header_with_description(self):
        """>ACC description ... normalizes to ACC."""
        assert normalize_sequence_id(">HQ641676.1 Genypterus capensis voucher SAIAB 79629") == "HQ641676.1"

    def test_fasta_header_gi_with_description(self):
        """>gi|...|gb|ACC| description normalizes to ACC."""
        assert normalize_sequence_id(">gi|336317909|gb|HQ641676.1| Genypterus capensis") == "HQ641676.1"

    def test_mafft_leading_r_prefix(self):
        """_R_ACC (MAFFT reverse-complement marker) normalizes to ACC."""
        assert normalize_sequence_id("_R_HQ641676.1") == "HQ641676.1"

    def test_mafft_trailing_r_suffix(self):
        """ACC_R_ (MAFFT reverse-complement marker) normalizes to ACC."""
        assert normalize_sequence_id("HQ641676.1_R_") == "HQ641676.1"

    def test_empty_string(self):
        """Empty string is returned as-is."""
        assert normalize_sequence_id("") == ""

    def test_plain_accession_gi_gb_round_trip(self):
        """gi|...|gb|ACC| and plain ACC normalize to the same value."""
        assert normalize_sequence_id("gi|336317909|gb|HQ641676.1|") == normalize_sequence_id("HQ641676.1")


class TestExtractAccessionFromPipeId:
    """Test cases for the _extract_accession_from_pipe_id private helper."""

    def test_gi_gb_format(self):
        """gi|...|gb|ACC| extracts ACC."""
        assert _extract_accession_from_pipe_id("gi|2273658778|gb|MZ387488.1|") == "MZ387488.1"

    def test_ref_format(self):
        """ref|ACC| extracts ACC."""
        assert _extract_accession_from_pipe_id("ref|NZ_CP123456.1|") == "NZ_CP123456.1"

    def test_gb_format(self):
        """gb|ACC| extracts ACC."""
        assert _extract_accession_from_pipe_id("gb|MZ387488.1|") == "MZ387488.1"

    def test_gnl_format(self):
        """gnl|database|identifier extracts identifier."""
        assert _extract_accession_from_pipe_id("gnl|BL_ORD_ID|12345") == "12345"

    def test_plain_accession_unchanged(self):
        """Plain accession with no pipes is returned unchanged."""
        assert _extract_accession_from_pipe_id("MZ387488.1") == "MZ387488.1"

    def test_custom_pipe_id_no_fallback(self):
        """Custom pipe-delimited IDs with no known prefix are returned unchanged (no fallback)."""
        assert _extract_accession_from_pipe_id("sampleA|42") == "sampleA|42"

    def test_custom_pipe_ids_not_conflated(self):
        """Two custom IDs sharing a trailing segment do NOT normalize to the same value."""
        a = _extract_accession_from_pipe_id("sampleA|42")
        b = _extract_accession_from_pipe_id("sampleB|42")
        assert a != b


class TestMMseqs2Runner:
    """Test cases for MMseqs2Runner class."""

    def test_init_default(self):
        """Test MMseqs2Runner initialization with defaults."""
        with patch.dict('os.environ', {}, clear=True):
            runner = MMseqs2Runner()
            assert runner.db_path == Path.home() / "mmseqs2db"

    def test_init_with_env_var(self):
        """Test MMseqs2Runner initialization with MMSEQS_DB_DIR env var."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.dict('os.environ', {'MMSEQS_DB_DIR': tmpdir}):
                runner = MMseqs2Runner()
                assert runner.db_path == Path(tmpdir)

    def test_init_with_legacy_env_var(self):
        """Test MMseqs2Runner initialization with legacy MMSEQS2DB env var."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.dict('os.environ', {'MMSEQS2DB': tmpdir}, clear=True):
                runner = MMseqs2Runner()
                assert runner.db_path == Path(tmpdir)

    def test_init_prefers_mmseqs_db_dir_over_legacy(self):
        """Test MMSEQS_DB_DIR takes precedence over MMSEQS2DB when both are set."""
        with tempfile.TemporaryDirectory() as preferred, tempfile.TemporaryDirectory() as legacy:
            with patch.dict(
                'os.environ',
                {'MMSEQS_DB_DIR': preferred, 'MMSEQS2DB': legacy},
                clear=True
            ):
                runner = MMseqs2Runner()
                assert runner.db_path == Path(preferred)

    def test_init_with_path(self):
        """Test MMseqs2Runner initialization with explicit path."""
        test_path = Path("/test/mmseqs2db")
        runner = MMseqs2Runner(db_path=test_path)
        assert runner.db_path == test_path

    @patch('subprocess.run')
    def test_check_mmseqs_available_true(self, mock_run):
        """Test checking mmseqs availability when available."""
        mock_run.return_value = MagicMock(returncode=0)
        runner = MMseqs2Runner()
        assert runner.check_mmseqs_available() is True

    @patch('subprocess.run')
    def test_check_mmseqs_available_false(self, mock_run):
        """Test checking mmseqs availability when not available."""
        mock_run.side_effect = FileNotFoundError()
        runner = MMseqs2Runner()
        assert runner.check_mmseqs_available() is False

    def test_parse_mmseqs_results(self):
        """Test parsing MMseqs2 tabular output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            mmseqs_output = tmppath / "mmseqs_results.txt"

            # Create mock MMseqs2 output (14 tab-separated columns)
            content = (
                "seq1\tNC_001\t95.5\t200\t500\t1000\t1\t200\t100\t299\t1e-50\t250.0\t9606\tEukaryota\n"
                "seq1\tNC_002\t92.0\t180\t500\t900\t1\t180\t50\t229\t1e-45\t240.0\t590\tBacteria\n"
                "seq2\tNC_003\t88.5\t150\t400\t800\t1\t150\t200\t349\t1e-40\t230.0\t0\t\n"
            )
            mmseqs_output.write_text(content)

            runner = MMseqs2Runner()
            hits = runner.parse_mmseqs_results(mmseqs_output)

            assert len(hits) == 3
            assert hits[0].query_id == 'seq1'
            assert hits[0].subject_id == 'NC_001'
            assert hits[0].percent_identity == 95.5
            assert hits[0].alignment_length == 200
            assert hits[0].query_length == 500
            assert hits[0].evalue == 1e-50
            assert hits[0].bit_score == 250.0
            assert hits[0].query_coverage == pytest.approx(40.0, abs=0.1)
            assert hits[0].subject_taxid == "9606"
            assert hits[0].subject_taxids == "9606"
            # taxid "0" should remain "0"
            assert hits[2].subject_taxid == "0"

    def test_parse_mmseqs_results_no_taxonomy(self):
        """Test parsing MMseqs2 output without taxonomy columns."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            mmseqs_output = tmppath / "mmseqs_results.txt"

            # Only 12 columns (no taxid / taxlineage)
            content = (
                "seq1\tNC_001\t95.5\t200\t500\t1000\t1\t200\t100\t299\t1e-50\t250.0\n"
            )
            mmseqs_output.write_text(content)

            runner = MMseqs2Runner()
            hits = runner.parse_mmseqs_results(mmseqs_output)

            assert len(hits) == 1
            assert hits[0].subject_taxid == "0"
            assert hits[0].subject_taxids == "0"

    def test_parse_mmseqs_results_na_taxid(self):
        """Test that N/A taxid values are normalized to '0'."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            mmseqs_output = tmppath / "mmseqs_results.txt"

            content = (
                "seq1\tNC_001\t95.5\t200\t500\t1000\t1\t200\t100\t299\t1e-50\t250.0\tN/A\t\n"
            )
            mmseqs_output.write_text(content)

            runner = MMseqs2Runner()
            hits = runner.parse_mmseqs_results(mmseqs_output)

            assert hits[0].subject_taxid == "0"

    def test_parse_mmseqs_results_invalid_format(self):
        """Test parsing MMseqs2 output with too few fields."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            mmseqs_output = tmppath / "mmseqs_results.txt"

            mmseqs_output.write_text("seq1\tNC_001\t95.5\n")

            runner = MMseqs2Runner()
            with pytest.raises(ValueError, match="Invalid MMseqs2 output format"):
                runner.parse_mmseqs_results(mmseqs_output)

    def test_parse_mmseqs_results_nonexistent(self):
        """Test parsing a non-existent output file."""
        runner = MMseqs2Runner()
        with pytest.raises(FileNotFoundError):
            runner.parse_mmseqs_results(Path("/nonexistent/file.txt"))

    def test_filter_hits(self):
        """Test filtering MMseqs2 hits."""
        hits = [
            BlastHit(
                query_id='seq1', subject_id='NC_001',
                percent_identity=95.0, alignment_length=400,
                query_length=400, subject_length=400,
                query_start=1, query_end=400,
                subject_start=1, subject_end=400,
                evalue=0.0, bit_score=800,
                query_coverage=100.0,
                subject_taxid="9606", subject_taxids="9606"
            ),
            BlastHit(
                query_id='seq1', subject_id='NC_002',
                percent_identity=85.0,  # below threshold
                alignment_length=200,
                query_length=400, subject_length=400,
                query_start=1, query_end=200,
                subject_start=1, subject_end=200,
                evalue=1e-10, bit_score=300,
                query_coverage=50.0,
                subject_taxid="590", subject_taxids="590"
            ),
        ]

        runner = MMseqs2Runner()
        filtered = runner.filter_hits(hits, min_identity=90.0, min_coverage=80.0)

        assert len(filtered) == 1
        assert filtered[0].subject_id == 'NC_001'

    @patch('eplace_lib.blast_analysis.MMseqs2Runner.check_mmseqs_available')
    @patch('subprocess.run')
    def test_run_easy_search_success(self, mock_run, mock_check):
        """Test running mmseqs easy-search successfully."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            query_file = tmppath / "query.fasta"
            output_file = tmppath / "output.txt"

            query_file.write_text(">seq1\nATCG\n")

            mock_check.return_value = True
            mock_run.return_value = MagicMock(returncode=0, stderr="")

            runner = MMseqs2Runner(db_path=tmppath)
            success = runner.run_easy_search(query_file, output_file)

            assert success is True
            mock_run.assert_called_once()
            call_args = mock_run.call_args[0][0]
            assert call_args[0] == 'mmseqs'
            assert call_args[1] == 'easy-search'

    def test_run_easy_search_no_query_file(self):
        """Test running mmseqs with a non-existent query file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            runner = MMseqs2Runner()
            with pytest.raises(FileNotFoundError):
                runner.run_easy_search(
                    tmppath / "nonexistent.fasta",
                    tmppath / "output.txt"
                )

    @patch('eplace_lib.blast_analysis.MMseqs2Runner.check_mmseqs_available')
    def test_run_easy_search_not_available(self, mock_check):
        """Test running mmseqs when it is not installed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            query_file = tmppath / "query.fasta"
            query_file.write_text(">seq1\nATCG\n")

            mock_check.return_value = False

            runner = MMseqs2Runner()
            with pytest.raises(RuntimeError, match="mmseqs is not available"):
                runner.run_easy_search(query_file, tmppath / "output.txt")


class TestRunMMseqsSearch:
    """Test cases for the run_mmseqs_search convenience function."""

    @patch('eplace_lib.blast_analysis.MMseqs2Runner.run_easy_search')
    @patch('eplace_lib.blast_analysis.MMseqs2Runner.parse_mmseqs_results')
    @patch('eplace_lib.blast_analysis.MMseqs2Runner.filter_hits')
    def test_run_mmseqs_search_success(self, mock_filter, mock_parse, mock_run):
        """Test run_mmseqs_search convenience function."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            query_file = tmppath / "query.fasta"
            output_file = tmppath / "output.txt"

            query_file.write_text(">seq1\nATCG\n")

            mock_run.return_value = True
            mock_parse.return_value = [MagicMock()]
            mock_filter.return_value = [MagicMock()]

            success, hits = run_mmseqs_search(
                query_fasta=query_file,
                output_file=output_file
            )

            assert success is True
            assert len(hits) == 1

    @patch('eplace_lib.blast_analysis.MMseqs2Runner.run_easy_search')
    def test_run_mmseqs_search_fails(self, mock_run):
        """Test run_mmseqs_search when the search fails."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            query_file = tmppath / "query.fasta"
            output_file = tmppath / "output.txt"

            query_file.write_text(">seq1\nATCG\n")
            mock_run.return_value = False

            success, hits = run_mmseqs_search(
                query_fasta=query_file,
                output_file=output_file
            )

            assert success is False
            assert len(hits) == 0

    def test_run_mmseqs_search_skip_existing(self):
        """Test that run_mmseqs_search skips the search when output already exists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            query_file = tmppath / "query.fasta"
            output_file = tmppath / "output.txt"

            query_file.write_text(">seq1\nATCG\n")
            # Write a valid (but minimal) output file so parsing succeeds
            output_file.write_text(
                "seq1\tNC_001\t95.5\t200\t500\t1000\t1\t200\t100\t299\t1e-50\t250.0\t9606\t\n"
            )

            with patch('eplace_lib.blast_analysis.MMseqs2Runner.run_easy_search') as mock_run:
                success, hits = run_mmseqs_search(
                    query_fasta=query_file,
                    output_file=output_file,
                    skip_existing=True
                )
                mock_run.assert_not_called()
                assert success is True

    def test_run_mmseqs_search_invalid_sensitivity_too_low(self):
        """Test that sensitivity below 1.0 raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            query_file = tmppath / "query.fasta"
            output_file = tmppath / "output.txt"
            query_file.write_text(">seq1\nATCG\n")

            with pytest.raises(ValueError, match="sensitivity must be between"):
                run_mmseqs_search(
                    query_fasta=query_file,
                    output_file=output_file,
                    sensitivity=0.5
                )

    def test_run_mmseqs_search_invalid_sensitivity_too_high(self):
        """Test that sensitivity above 7.5 raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            query_file = tmppath / "query.fasta"
            output_file = tmppath / "output.txt"
            query_file.write_text(">seq1\nATCG\n")

            with pytest.raises(ValueError, match="sensitivity must be between"):
                run_mmseqs_search(
                    query_fasta=query_file,
                    output_file=output_file,
                    sensitivity=8.0
                )
