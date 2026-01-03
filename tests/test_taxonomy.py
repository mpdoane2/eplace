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
    process_blast_results_for_taxonomy,
    rewrite_blast_hits
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
        tax_info, phyla_info = self.taxonomy_extractor.parse_taxids([9606, 590])

        homo_tuple = ('9605', 'Homo')
        salm_tuple = ('590', 'Salmonella')
        homo_phylum = ('7711', 'Chordata')
        salm_phylum = ('1224', 'Pseudomonadota')
        assert tax_info["9606"] == homo_tuple
        assert tax_info["590"] == salm_tuple
        assert phyla_info['9606'] == homo_phylum
        assert phyla_info['590'] == salm_phylum
    
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
    
    @patch('eplace_lib.taxonomy.TaxonomyExtractor.parse_taxids')
    @patch('eplace_lib.taxonomy.SequenceExtractor.extract_representatives_for_query')
    def test_process_blast_results_for_taxonomy(self, mock_extract, mock_parse_taxids):
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
            
            # Mock parse_taxids to return taxonomy and phylum info
            mock_parse_taxids.return_value = (
                {
                    "149539": ("590", "Salmonella"),
                    "9606": ("9605", "Homo"),
                    "9597": ("9596", "Pan")
                },
                {
                    "149539": ("1224", "Pseudomonadota"),
                    "9606": ("7711", "Chordata"),
                    "9597": ("7711", "Chordata")
                }
            )
            
            mock_extract.return_value = tmppath / "output.fasta"
            
            results = process_blast_results_for_taxonomy(
                blast_hits=hits,
                output_dir=tmppath,
                rank='genus'
            )
            
            assert len(results) == 3
            assert 'seq1' in results
            assert 'seq2' in results
            
            # Verify phylum information is correctly set on blast hits
            # seq1 has taxid 149539 (Salmonella) which should map to Pseudomonadota phylum
            assert hits[0].subject_phylum_tid == '1224'
            assert hits[0].subject_phylum_name == 'Pseudomonadota'
            
            # seq2 has taxid 9606 (Homo sapiens) which should map to Chordata phylum
            assert hits[1].subject_phylum_tid == '7711'
            assert hits[1].subject_phylum_name == 'Chordata'
            
            # seq3 has taxid 9597 (Pan) which should also map to Chordata phylum
            assert hits[2].subject_phylum_tid == '7711'
            assert hits[2].subject_phylum_name == 'Chordata'

class TestRewriteBlastHits:
    """Test cases for rewrite_blast_hits function."""
    
    # Expected field names for blast hit output
    EXPECTED_FIELDS = [
        "query_id", "subject_id", "percent_identity", "alignment_length",
        "query_length", "subject_length", "query_start", "query_end",
        "subject_start", "subject_end", "evalue", "bit_score",
        "query_coverage", "subject_taxid", "subject_taxids",
        "subject_rank_tid", "subject_rank_name",
        "subject_phylum_tid", "subject_phylum_name"
    ]
    
    def test_rewrite_blast_hits_with_complete_annotations(self):
        """Test writing blast hits with all fields populated."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "blast_hits.tsv"
            
            hits = [
                BlastHit(
                    query_id='seq1',
                    subject_id='gi|156763568|gb|EU014687.1|',
                    percent_identity=100.0,
                    alignment_length=540,
                    query_length=540,
                    subject_length=1432,
                    query_start=1,
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=998.0,
                    query_coverage=100.0,
                    subject_taxid="149539",
                    subject_taxids="149539",
                    subject_rank_tid="590",
                    subject_rank_name="Salmonella",
                    subject_phylum_tid="1224",
                    subject_phylum_name="Pseudomonadota"
                )
            ]
            
            result = rewrite_blast_hits(hits, output_file, header=True)
            
            assert result is True
            assert output_file.exists()
            
            # Read and verify content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            assert len(lines) == 2  # header + 1 data line
            
            # Check header
            header = lines[0].strip().split('\t')
            assert header == self.EXPECTED_FIELDS
            
            # Check data line
            data = lines[1].strip().split('\t')
            assert data[0] == "seq1"
            assert data[1] == "gi|156763568|gb|EU014687.1|"
            assert data[15] == "590"
            assert data[16] == "Salmonella"
            assert data[17] == "1224"
            assert data[18] == "Pseudomonadota"
    
    def test_rewrite_blast_hits_with_none_values(self):
        """Test writing blast hits with None values in optional fields."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "blast_hits.tsv"
            
            hits = [
                BlastHit(
                    query_id='seq1',
                    subject_id='gi|156763568|gb|EU014687.1|',
                    percent_identity=100.0,
                    alignment_length=540,
                    query_length=540,
                    subject_length=1432,
                    query_start=1,
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=998.0,
                    query_coverage=100.0,
                    subject_taxid="149539",
                    subject_taxids="149539",
                    subject_rank_tid=None,
                    subject_rank_name=None,
                    subject_phylum_tid=None,
                    subject_phylum_name=None
                )
            ]
            
            result = rewrite_blast_hits(hits, output_file, header=True)
            
            assert result is True
            assert output_file.exists()
            
            # Read and verify content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            assert len(lines) == 2  # header + 1 data line
            
            # Check data line - None values should be empty strings
            # Remove only newline to preserve tab-delimited format
            data = lines[1].rstrip('\n').split('\t')
            assert len(data) == 19  # Verify all fields present
            assert data[0] == "seq1"
            assert data[15] == ""  # subject_rank_tid
            assert data[16] == ""  # subject_rank_name
            assert data[17] == ""  # subject_phylum_tid
            assert data[18] == ""  # subject_phylum_name
    
    def test_rewrite_blast_hits_without_header(self):
        """Test writing blast hits without header."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "blast_hits.tsv"
            
            hits = [
                BlastHit(
                    query_id='seq1',
                    subject_id='gi|156763568|gb|EU014687.1|',
                    percent_identity=100.0,
                    alignment_length=540,
                    query_length=540,
                    subject_length=1432,
                    query_start=1,
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=998.0,
                    query_coverage=100.0,
                    subject_taxid="149539",
                    subject_taxids="149539",
                    subject_rank_tid="590",
                    subject_rank_name="Salmonella",
                    subject_phylum_tid="1224",
                    subject_phylum_name="Pseudomonadota"
                )
            ]
            
            result = rewrite_blast_hits(hits, output_file, header=False)
            
            assert result is True
            assert output_file.exists()
            
            # Read and verify content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            # No header, just data line
            assert len(lines) == 1
            
            # Check data line
            data = lines[0].strip().split('\t')
            assert data[0] == "seq1"
            assert data[1] == "gi|156763568|gb|EU014687.1|"
    
    def test_rewrite_blast_hits_empty_list(self):
        """Test writing blast hits with empty list."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "blast_hits.tsv"
            
            hits = []
            
            result = rewrite_blast_hits(hits, output_file, header=True)
            
            assert result is True
            assert output_file.exists()
            
            # Read and verify content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            # Should only have header, no data lines
            assert len(lines) == 1
            
            # Check header exists with all expected fields
            header = lines[0].strip().split('\t')
            assert header == self.EXPECTED_FIELDS
    
    def test_rewrite_blast_hits_multiple_records(self):
        """Test writing multiple blast hits."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            output_file = tmppath / "blast_hits.tsv"
            
            hits = [
                BlastHit(
                    query_id='seq1',
                    subject_id='gi|156763568|gb|EU014687.1|',
                    percent_identity=100.0,
                    alignment_length=540,
                    query_length=540,
                    subject_length=1432,
                    query_start=1,
                    query_end=540,
                    subject_start=1,
                    subject_end=540,
                    evalue=0.0,
                    bit_score=998.0,
                    query_coverage=100.0,
                    subject_taxid="149539",
                    subject_taxids="149539",
                    subject_rank_tid="590",
                    subject_rank_name="Salmonella",
                    subject_phylum_tid="1224",
                    subject_phylum_name="Pseudomonadota"
                ),
                BlastHit(
                    query_id='seq2',
                    subject_id='gi|34190046|gb|BC014593.2|',
                    percent_identity=95.5,
                    alignment_length=420,
                    query_length=420,
                    subject_length=784,
                    query_start=1,
                    query_end=420,
                    subject_start=1,
                    subject_end=420,
                    evalue=1e-100,
                    bit_score=776.0,
                    query_coverage=100.0,
                    subject_taxid="9606",
                    subject_taxids="9606",
                    subject_rank_tid="9605",
                    subject_rank_name="Homo",
                    subject_phylum_tid="7711",
                    subject_phylum_name="Chordata"
                )
            ]
            
            result = rewrite_blast_hits(hits, output_file, header=True)
            
            assert result is True
            assert output_file.exists()
            
            # Read and verify content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            assert len(lines) == 3  # header + 2 data lines
            
            # Check first data line
            data1 = lines[1].strip().split('\t')
            assert data1[0] == "seq1"
            assert data1[17] == "1224"
            assert data1[18] == "Pseudomonadota"
            
            # Check second data line
            data2 = lines[2].strip().split('\t')
            assert data2[0] == "seq2"
            assert data2[17] == "7711"
            assert data2[18] == "Chordata"

