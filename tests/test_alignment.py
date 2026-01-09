"""
Tests for sequence alignment and phylogenetic tree building module.
"""

import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

from eplace_lib.blast_analysis import BlastHit
from eplace_lib.alignment import (
    SequenceTrimmer,
    MAFFTAligner,
    IQTreeBuilder,
    process_query_alignment_and_tree
)


class TestSequenceTrimmer:
    """Test cases for SequenceTrimmer class."""

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

    def test_trim_sequence_by_coordinates_forward(self):
        """Test trimming a sequence with forward strand coordinates."""
        sequence = "ATCGATCGATCGATCGATCG"  # 20bp
        # Trim positions 5-15 (1-indexed)
        trimmed = SequenceTrimmer.trim_sequence_by_coordinates(sequence, 5, 15)
        # Expected: positions 5-15 inclusive = "ATCGATCGATC" (11bp)
        # Position 5 in 1-indexed = index 4 in 0-indexed = 'A'
        assert trimmed == "ATCGATCGATC"
        assert len(trimmed) == 11
    
    def test_trim_sequence_by_coordinates_reverse(self):
        """Test trimming a sequence with reverse strand coordinates."""
        sequence = "ATCGATCGATCGATCGATCG"  # 20bp
        # Reverse strand: start > end
        trimmed = SequenceTrimmer.trim_sequence_by_coordinates(sequence, 15, 5)
        # Should handle reverse coordinates
        assert len(trimmed) == 11
    
    def test_trim_sequence_by_coordinates_boundaries(self):
        """Test trimming at sequence boundaries."""
        sequence = "ATCGATCGATCG"  # 12bp
        # Full sequence
        trimmed = SequenceTrimmer.trim_sequence_by_coordinates(sequence, 1, 12)
        assert trimmed == sequence
        
        # Beyond boundaries (should be clamped)
        trimmed = SequenceTrimmer.trim_sequence_by_coordinates(sequence, 1, 100)
        assert trimmed == sequence
    
    def test_trim_sequences_from_blast_hits(self):
        """Test trimming sequences based on BLAST hits."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create test FASTA with query and subject sequences
            fasta_path = tmppath / "test.fasta"
            fasta_content = """>query1
ATCGATCGATCGATCGATCGATCGATCGATCGATCG
>subject1
NNNNNATCGATCGATCGATCGATCGATCGNNNNNNNNNN
>subject2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
"""
            fasta_path.write_text(fasta_content)
            
            # Create BLAST hits
            blast_hits = [
                BlastHit(
                    query_id='query1',
                    subject_id='subject1',
                    percent_identity=95.0,
                    alignment_length=30,
                    query_length=36,
                    subject_length=40,
                    query_start=1,
                    query_end=30,
                    subject_start=6,
                    subject_end=35,
                    evalue=1e-10,
                    bit_score=100,
                    query_coverage=83.3,
                    subject_taxid="590",
                    subject_taxids="590",
                    subject_taxonomy=self.salmonella_taxonomy
                ),
                BlastHit(
                    query_id='query1',
                    subject_id='subject2',
                    percent_identity=90.0,
                    alignment_length=20,
                    query_length=36,
                    subject_length=30,
                    query_start=5,
                    query_end=24,
                    subject_start=1,
                    subject_end=20,
                    evalue=1e-8,
                    bit_score=80,
                    query_coverage=55.6,
                    subject_taxid="9606",
                    subject_taxids="9606;9605",
                    subject_taxonomy=self.human_taxonomy
                )
            ]
            
            # Trim sequences
            output_fasta = tmppath / "trimmed.fasta"
            success = SequenceTrimmer.trim_sequences_from_blast_hits(
                fasta_path=fasta_path,
                blast_hits=blast_hits,
                taxonomic_rank='genus',
                output_fasta=output_fasta,
                query_id='query1'
            )
            
            assert success is True
            assert output_fasta.exists()
            
            # Read and verify output
            content = output_fasta.read_text()
            assert '>query1' in content
            assert '>subject1 Salmonella' in content
            assert '>subject2 Homo' in content
    
    def test_trim_sequences_with_gi_format_blast_hits(self):
        """Test trimming when BLAST hits have gi|...|gb|...| format but FASTA has accessions."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create test FASTA with accession-only IDs (like blastdbcmd output)
            fasta_path = tmppath / "test.fasta"
            fasta_content = """>query1
ATCGATCGATCGATCGATCGATCGATCGATCGATCG
>MZ387488.1
NNNNNATCGATCGATCGATCGATCGATCGNNNNNNNNNN
>NZ_CP123456.1
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
"""
            fasta_path.write_text(fasta_content)
            
            # Create BLAST hits with full gi|...|gb|...| format
            blast_hits = [
                BlastHit(
                    query_id='query1',
                    subject_id='gi|2273658778|gb|MZ387488.1|',  # Full format
                    percent_identity=95.0,
                    alignment_length=30,
                    query_length=36,
                    subject_length=40,
                    query_start=1,
                    query_end=30,
                    subject_start=6,
                    subject_end=35,
                    evalue=1e-10,
                    bit_score=100,
                    query_coverage=83.3,
                    subject_taxid="590",
                    subject_taxids="590",
                    subject_taxonomy=self.salmonella_taxonomy
                ),
                BlastHit(
                    query_id='query1',
                    subject_id='ref|NZ_CP123456.1|',  # ref format
                    percent_identity=90.0,
                    alignment_length=20,
                    query_length=36,
                    subject_length=30,
                    query_start=5,
                    query_end=24,
                    subject_start=1,
                    subject_end=20,
                    evalue=1e-8,
                    bit_score=80,
                    query_coverage=55.6,
                    subject_taxid="590",
                    subject_taxids="590",
                    subject_taxonomy=self.human_taxonomy
                )
            ]
            
            # Trim sequences - should match accessions correctly
            output_fasta = tmppath / "trimmed.fasta"
            success = SequenceTrimmer.trim_sequences_from_blast_hits(
                fasta_path=fasta_path,
                blast_hits=blast_hits,
                output_fasta=output_fasta,
                taxonomic_rank='genus',
                query_id='query1'
            )
            
            assert success is True
            assert output_fasta.exists()
            
            # Read and verify output
            content = output_fasta.read_text()
            assert '>query1' in content
            assert '>MZ387488.1 Salmonella' in content
            assert '>NZ_CP123456.1 Homo' in content


class TestMAFFTAligner:
    """Test cases for MAFFTAligner class."""
    
    def test_check_mafft_available(self):
        """Test checking if MAFFT is available."""
        # This will return True or False depending on system
        result = MAFFTAligner.check_mafft_available()
        assert isinstance(result, bool)
    
    @patch('subprocess.run')
    def test_align_sequences_success(self, mock_run):
        """Test successful MAFFT alignment."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create test input
            input_fasta = tmppath / "input.fasta"
            input_fasta.write_text(">seq1\nATCG\n>seq2\nATCG\n")
            
            output_fasta = tmppath / "output.fasta"
            
            # Mock successful MAFFT run
            mock_run.return_value = MagicMock(returncode=0, stderr="")
            
            # Mock check_mafft_available to return True
            with patch.object(MAFFTAligner, 'check_mafft_available', return_value=True):
                result = MAFFTAligner.align_sequences(
                    input_fasta=input_fasta,
                    output_fasta=output_fasta
                )
            
            # Verify MAFFT was called
            assert mock_run.called
            # Check that output file path was used
            call_args = mock_run.call_args
            # The output file is passed as stdout in the actual implementation
    
    def test_align_sequences_mafft_not_available(self):
        """Test alignment fails when MAFFT is not available."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            input_fasta = tmppath / "input.fasta"
            input_fasta.write_text(">seq1\nATCG\n")
            
            output_fasta = tmppath / "output.fasta"
            
            # Mock MAFFT not available
            with patch.object(MAFFTAligner, 'check_mafft_available', return_value=False):
                result = MAFFTAligner.align_sequences(
                    input_fasta=input_fasta,
                    output_fasta=output_fasta
                )
            
            assert result is False
    
    def test_align_sequences_input_not_found(self):
        """Test alignment fails when input file doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            input_fasta = tmppath / "nonexistent.fasta"
            output_fasta = tmppath / "output.fasta"
            
            with patch.object(MAFFTAligner, 'check_mafft_available', return_value=True):
                result = MAFFTAligner.align_sequences(
                    input_fasta=input_fasta,
                    output_fasta=output_fasta
                )
            
            assert result is False


class TestIQTreeBuilder:
    """Test cases for IQTreeBuilder class."""

    def setup_method(self):
        self.salmonella_taxonomy = {
            'phylum': ('1224', 'Pseudomonadota'),
            'class': ('1236', 'Gammaproteobacteria'),
            'order': ('91347', 'Enterobacterales'),
            'family': ('543', 'Enterobacteriaceae'),
            'genus': ('590', 'Salmonella'),
            'species': ('28901', 'Salmonella enterica')
        }
        self.ecoli_taxonomy = {
            'phylum': ('1224', 'Pseudomonadota'),
            'class': ('1236', 'Gammaproteobacteria'),
            'order': ('91347', 'Enterobacterales'),
            'family': ('543', 'Enterobacteriaceae'),
            'genus': ('561', 'Escherichia'),
            'species': ('562', 'Escherichia coli')
        }

    def test_check_iqtree_available(self):
        """Test checking if IQTree is available."""
        available, cmd = IQTreeBuilder.check_iqtree_available()
        assert isinstance(available, bool)
        if available:
            assert cmd in ['iqtree', 'iqtree2', 'iqtree3']
        else:
            assert cmd is None
    
    @patch('subprocess.run')
    def test_build_tree_success(self, mock_run):
        """Test successful tree building."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create test alignment
            alignment_fasta = tmppath / "alignment.fasta"
            alignment_fasta.write_text(">seq1\nATCG\n>seq2\nATCG\n")
            
            output_prefix = tmppath / "tree"
            tree_file = Path(str(output_prefix) + ".treefile")
            
            # Mock successful IQTree run
            mock_run.return_value = MagicMock(returncode=0, stderr="")
            
            # Create the expected tree file
            tree_file.write_text("(seq1:0.1,seq2:0.1);")
            
            # Mock check_iqtree_available
            with patch.object(IQTreeBuilder, 'check_iqtree_available', return_value=(True, 'iqtree2')):
                result = IQTreeBuilder.build_tree(
                    alignment_fasta=alignment_fasta,
                    output_prefix=output_prefix
                )
            
            assert result is True
            assert mock_run.called
    
    def test_build_tree_iqtree_not_available(self):
        """Test tree building fails when IQTree is not available."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            alignment_fasta = tmppath / "alignment.fasta"
            alignment_fasta.write_text(">seq1\nATCG\n")
            
            output_prefix = tmppath / "tree"
            
            # Mock IQTree not available
            with patch.object(IQTreeBuilder, 'check_iqtree_available', return_value=(False, None)):
                result = IQTreeBuilder.build_tree(
                    alignment_fasta=alignment_fasta,
                    output_prefix=output_prefix
                )
            
            assert result is False
    
    def test_relabel_tree_with_taxonomy(self):
        """Test relabeling tree with taxonomic names."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create test tree
            tree_file = tmppath / "tree.treefile"
            tree_content = "(subject1:0.1,subject2:0.2,query1:0.0);"
            tree_file.write_text(tree_content)
            
            # Create BLAST hits with taxonomy
            blast_hits = [
                BlastHit(
                    query_id='query1',
                    subject_id='subject1',
                    percent_identity=95.0,
                    alignment_length=30,
                    query_length=36,
                    subject_length=40,
                    query_start=1,
                    query_end=30,
                    subject_start=6,
                    subject_end=35,
                    evalue=1e-10,
                    bit_score=100,
                    query_coverage=83.3,
                    subject_taxid="12345",
                    subject_taxids="12345",
                    subject_taxonomy=self.salmonella_taxonomy
                ),
                BlastHit(
                    query_id='query1',
                    subject_id='subject2',
                    percent_identity=90.0,
                    alignment_length=20,
                    query_length=36,
                    subject_length=30,
                    query_start=5,
                    query_end=24,
                    subject_start=1,
                    subject_end=20,
                    evalue=1e-8,
                    bit_score=80,
                    query_coverage=55.6,
                    subject_taxid="67890",
                    subject_taxids="67890",
                    subject_taxonomy=self.ecoli_taxonomy
                )
            ]
            
            output_tree = tmppath / "tree_labeled.treefile"
            
            # Relabel tree
            success = IQTreeBuilder.relabel_tree_with_taxonomy(
                tree_file=tree_file,
                blast_hits=blast_hits,
                output_tree=output_tree,
                taxonomic_rank='species'
            )
            
            assert success is True
            assert output_tree.exists()
            
            # Verify labels were replaced
            content = output_tree.read_text()
            assert 'Salmonella' in content
            assert 'Escherichia_coli' in content  # Spaces replaced with underscores
    
    def test_relabel_tree_with_reversed_sequences(self):
        """Test relabeling tree with reversed sequences (MAFFT _R_ prefix)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create test tree with _R_ prefix on subject2 (reversed by MAFFT)
            tree_file = tmppath / "tree.treefile"
            tree_content = "(subject1:0.1,_R_subject2:0.2,query1:0.0);"
            tree_file.write_text(tree_content)
            
            # Create BLAST hits with taxonomy
            blast_hits = [
                BlastHit(
                    query_id='query1',
                    subject_id='subject1',
                    percent_identity=95.0,
                    alignment_length=30,
                    query_length=36,
                    subject_length=40,
                    query_start=1,
                    query_end=30,
                    subject_start=6,
                    subject_end=35,
                    evalue=1e-10,
                    bit_score=100,
                    query_coverage=83.3,
                    subject_taxid="12345",
                    subject_taxids="12345",
                    subject_taxonomy=self.salmonella_taxonomy
                ),
                BlastHit(
                    query_id='query1',
                    subject_id='subject2',
                    percent_identity=90.0,
                    alignment_length=20,
                    query_length=36,
                    subject_length=30,
                    query_start=5,
                    query_end=24,
                    subject_start=1,
                    subject_end=20,
                    evalue=1e-8,
                    bit_score=80,
                    query_coverage=55.6,
                    subject_taxid="67890",
                    subject_taxids="67890",
                    subject_taxonomy=self.ecoli_taxonomy
                )
            ]
            
            output_tree = tmppath / "tree_labeled.treefile"
            
            # Relabel tree
            success = IQTreeBuilder.relabel_tree_with_taxonomy(
                tree_file=tree_file,
                blast_hits=blast_hits,
                output_tree=output_tree,
                taxonomic_rank='species'
            )
            
            assert success is True
            assert output_tree.exists()
            
            # Verify labels were replaced
            content = output_tree.read_text()
            # subject1 should be normal
            assert 'Salmonella_enterica:0.1' in content
            # subject2 should have _R suffix (reversed)
            assert 'Escherichia_coli_R:0.2' in content
            # query should remain unchanged
            assert 'query1:0.0' in content



class TestProcessQueryAlignmentAndTree:
    """Test cases for the complete pipeline."""

    def setup_method(self):
        self.salmonella_taxonomy = {
            'phylum': ('1224', 'Pseudomonadota'),
            'class': ('1236', 'Gammaproteobacteria'),
            'order': ('91347', 'Enterobacterales'),
            'family': ('543', 'Enterobacteriaceae'),
            'genus': ('590', 'Salmonella')
        }

    @patch('eplace_lib.alignment.MAFFTAligner.align_sequences')
    @patch('eplace_lib.alignment.IQTreeBuilder.build_tree')
    @patch('eplace_lib.alignment.IQTreeBuilder.relabel_tree_with_taxonomy')
    def test_process_query_alignment_and_tree(
        self,
        mock_relabel,
        mock_build_tree,
        mock_align
    ):
        """Test the complete pipeline for a single query."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create query directory and files
            query_id = 'test_query'
            query_dir = tmppath / query_id
            query_dir.mkdir()
            
            # Create representatives file
            representatives_fasta = query_dir / f"{query_id}_representatives.fasta"
            representatives_fasta.write_text(">subject1\nATCGATCG\n")
            
            # Create query FASTA
            query_fasta = tmppath / "query.fasta"
            query_fasta.write_text(f">{query_id}\nATCGATCG\n")
            
            # Create BLAST hits
            blast_hits = [
                BlastHit(
                    query_id=query_id,
                    subject_id='subject1',
                    percent_identity=95.0,
                    alignment_length=8,
                    query_length=8,
                    subject_length=8,
                    query_start=1,
                    query_end=8,
                    subject_start=1,
                    subject_end=8,
                    evalue=1e-10,
                    bit_score=100,
                    query_coverage=100.0,
                    subject_taxid="590",
                    subject_taxids="590",
                    subject_taxonomy=self.salmonella_taxonomy
                )
            ]
            
            # Mock successful operations
            mock_align.return_value = True
            mock_build_tree.return_value = True
            mock_relabel.return_value = True
            
            # Create expected tree file
            tree_file = query_dir / f"{query_id}_tree.treefile"
            tree_file.write_text("(subject1:0.1,test_query:0.0);")
            
            # Run pipeline
            results = process_query_alignment_and_tree(
                query_id=query_id,
                query_dir=query_dir,
                blast_hits=blast_hits,
                query_fasta=query_fasta,
                taxonomic_rank='genus',
                num_threads=1
            )
            
            # Verify results
            assert results['trimmed_fasta'] is not None
            assert results['alignment'] is not None
            assert results['tree'] is not None
            
            # Verify trimmed file was created
            assert results['trimmed_fasta'].exists()
