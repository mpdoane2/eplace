"""
Tests for NCBI database download functionality.
"""

import os
import tempfile
import hashlib
import tarfile
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
import pytest

from eplace_lib.ncbi_download import (
    NCBIDownloader,
    MMseqsDownloader,
    setup_ncbi_database,
    setup_mmseqs_database,
    setup_mmseqs_taxonomy,
    check_available_memory_gb
)


class TestNCBIDownloader:
    """Test cases for NCBIDownloader class."""
    
    def test_init(self):
        """Test NCBIDownloader initialization."""
        downloader = NCBIDownloader()
        assert downloader.blastdb_dir is None
        assert downloader.NCBI_FTP_BASE == "https://ftp.ncbi.nlm.nih.gov/blast/db/"
        assert downloader.CORE_NT_PREFIX == "core_nt"
    
    def test_get_blastdb_directory_with_env_var(self):
        """Test get_blastdb_directory when BLASTDB env var is set."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.dict(os.environ, {'BLASTDB': tmpdir}):
                downloader = NCBIDownloader()
                result = downloader.get_blastdb_directory()
                assert result == Path(tmpdir)
                assert result.exists()
    
    def test_get_blastdb_directory_without_env_var(self):
        """Test get_blastdb_directory when BLASTDB env var is not set."""
        with patch.dict(os.environ, {}, clear=True):
            # Remove BLASTDB if it exists
            os.environ.pop('BLASTDB', None)
            downloader = NCBIDownloader()
            result = downloader.get_blastdb_directory()
            expected = Path.home() / "blastdb"
            assert result == expected
    
    def test_get_blastdb_directory_with_invalid_env_var(self):
        """Test get_blastdb_directory when BLASTDB points to non-existent dir."""
        with patch.dict(os.environ, {'BLASTDB': '/nonexistent/path'}):
            downloader = NCBIDownloader()
            result = downloader.get_blastdb_directory()
            expected = Path.home() / "blastdb"
            assert result == expected
    
    def test_check_database_exists_true(self):
        """Test check_database_exists when database files exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            # Create a mock database file
            (tmppath / "core_nt.00.nhr").touch()
            
            downloader = NCBIDownloader()
            result = downloader.check_database_exists(tmppath)
            assert result is True
    
    def test_check_database_exists_false(self):
        """Test check_database_exists when no database files exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            downloader = NCBIDownloader()
            result = downloader.check_database_exists(tmppath)
            assert result is False
    
    def test_check_database_exists_only_tarballs(self):
        """Test check_database_exists when only tarballs exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            # Create only tarball, not extracted files
            (tmppath / "core_nt.00.tar.gz").touch()
            
            downloader = NCBIDownloader()
            result = downloader.check_database_exists(tmppath)
            assert result is False
    
    def test_verify_md5_success(self):
        """Test MD5 verification with correct checksum."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create a test file
            test_file = tmppath / "test.txt"
            test_content = b"Hello, World!"
            test_file.write_bytes(test_content)
            
            # Calculate MD5
            md5_hash = hashlib.md5(test_content).hexdigest()
            
            # Create MD5 file
            md5_file = tmppath / "test.txt.md5"
            md5_file.write_text(f"{md5_hash}  test.txt\n")
            
            downloader = NCBIDownloader()
            result = downloader.verify_md5(test_file, md5_file)
            assert result is True
    
    def test_verify_md5_failure(self):
        """Test MD5 verification with incorrect checksum."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create a test file
            test_file = tmppath / "test.txt"
            test_file.write_bytes(b"Hello, World!")
            
            # Create MD5 file with wrong hash (valid format but wrong value)
            md5_file = tmppath / "test.txt.md5"
            md5_file.write_text("00000000000000000000000000000000  test.txt\n")
            
            downloader = NCBIDownloader()
            result = downloader.verify_md5(test_file, md5_file)
            assert result is False
    
    def test_verify_md5_invalid_format(self):
        """Test MD5 verification with invalid file format."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create a test file
            test_file = tmppath / "test.txt"
            test_file.write_bytes(b"Hello, World!")
            
            downloader = NCBIDownloader()
            
            # Test empty MD5 file
            md5_file = tmppath / "empty.md5"
            md5_file.write_text("")
            with pytest.raises(ValueError, match="MD5 file is empty"):
                downloader.verify_md5(test_file, md5_file)
            
            # Test invalid hash format (too short)
            md5_file.write_text("tooshort")
            with pytest.raises(ValueError, match="Invalid MD5 hash format"):
                downloader.verify_md5(test_file, md5_file)
            
            # Test invalid characters
            md5_file.write_text("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
            with pytest.raises(ValueError, match="Invalid MD5 hash format"):
                downloader.verify_md5(test_file, md5_file)
    
    def test_extract_tarball(self):
        """Test tarball extraction."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create a test file to add to tarball
            test_file = tmppath / "test.txt"
            test_file.write_text("Test content")
            
            # Create tarball
            tarball_path = tmppath / "test.tar.gz"
            with tarfile.open(tarball_path, 'w:gz') as tar:
                tar.add(test_file, arcname="test.txt")
            
            # Remove original file
            test_file.unlink()
            
            # Extract
            extract_dir = tmppath / "extracted"
            extract_dir.mkdir()
            
            downloader = NCBIDownloader()
            downloader.extract_tarball(tarball_path, extract_dir)
            
            # Verify extraction
            extracted_file = extract_dir / "test.txt"
            assert extracted_file.exists()
            assert extracted_file.read_text() == "Test content"
    
    def test_extract_tarball_with_unsafe_path(self):
        """Test tarball extraction rejects unsafe paths."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            # Create a test file
            test_file = tmppath / "test.txt"
            test_file.write_text("Test content")
            
            # Create tarball with path traversal
            tarball_path = tmppath / "unsafe.tar.gz"
            with tarfile.open(tarball_path, 'w:gz') as tar:
                # Use a path that escapes the extraction directory
                tar.add(test_file, arcname="../evil.txt")
            
            # Try to extract
            extract_dir = tmppath / "extracted"
            extract_dir.mkdir()
            
            downloader = NCBIDownloader()
            # Should raise ValueError for unsafe path
            with pytest.raises(ValueError, match="Unsafe path in tarball"):
                downloader.extract_tarball(tarball_path, extract_dir)
    
    def test_download_file_with_unsafe_filename(self):
        """Test download_file rejects unsafe filenames."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            
            downloader = NCBIDownloader()
            
            # Test path traversal attempts
            with pytest.raises(ValueError, match="Invalid filename"):
                downloader.download_file("../../../etc/passwd", tmppath)
            
            with pytest.raises(ValueError, match="Invalid filename"):
                downloader.download_file("/etc/passwd", tmppath)
            
            with pytest.raises(ValueError, match="Invalid filename"):
                downloader.download_file("..\\..\\windows\\system32", tmppath)
    
    @patch('eplace_lib.ncbi_download.urlopen')
    def test_get_available_files(self, mock_urlopen):
        """Test getting available files from FTP server."""
        # Mock HTML response
        mock_html = '''
        <html>
        <a href="core_nt.00.tar.gz">core_nt.00.tar.gz</a>
        <a href="core_nt.00.tar.gz.md5">core_nt.00.tar.gz.md5</a>
        <a href="core_nt.01.tar.gz">core_nt.01.tar.gz</a>
        <a href="core_nt.01.tar.gz.md5">core_nt.01.tar.gz.md5</a>
        <a href="other_file.txt">other_file.txt</a>
        </html>
        '''
        
        mock_response = MagicMock()
        mock_response.read.return_value = mock_html.encode('utf-8')
        mock_response.__enter__.return_value = mock_response
        mock_response.__exit__.return_value = None
        mock_urlopen.return_value = mock_response
        
        downloader = NCBIDownloader()
        files = downloader.get_available_files()
        
        assert 'core_nt.00.tar.gz' in files
        assert 'core_nt.00.tar.gz.md5' in files
        assert 'core_nt.01.tar.gz' in files
        assert 'core_nt.01.tar.gz.md5' in files
        assert 'other_file.txt' not in files


class TestSetupNCBIDatabase:
    """Test cases for setup_ncbi_database convenience function."""
    
    @patch('eplace_lib.ncbi_download.NCBIDownloader.download_and_setup_database')
    def test_setup_ncbi_database(self, mock_setup):
        """Test setup_ncbi_database convenience function."""
        mock_setup.return_value = (True, "Success")
        
        success, message = setup_ncbi_database()
        
        assert success is True
        assert message == "Success"
        mock_setup.assert_called_once_with(False, True)
    
    @patch('eplace_lib.ncbi_download.NCBIDownloader.download_and_setup_database')
    def test_setup_ncbi_database_with_force(self, mock_setup):
        """Test setup_ncbi_database with force_download."""
        mock_setup.return_value = (True, "Success")
        
        success, message = setup_ncbi_database(force_download=True)
        
        assert success is True
        mock_setup.assert_called_once_with(True, True)
    
    @patch('eplace_lib.ncbi_download.NCBIDownloader.download_and_setup_database')
    def test_setup_ncbi_database_with_verbose(self, mock_setup):
        """Test setup_ncbi_database with verbose parameter."""
        mock_setup.return_value = (True, "Success")
        
        success, message = setup_ncbi_database(verbose=False)
        
        assert success is True
        mock_setup.assert_called_once_with(False, False)


class TestMMseqsDownloader:
    """Test cases for MMseqsDownloader class and helpers."""

    def test_get_mmseqsdb_directory_with_preferred_env(self):
        """MMSEQS_DB_DIR should be used when present."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.dict(os.environ, {'MMSEQS_DB_DIR': tmpdir}, clear=True):
                downloader = MMseqsDownloader()
                assert downloader.get_mmseqsdb_directory() == Path(tmpdir)

    def test_get_mmseqsdb_directory_with_legacy_env(self):
        """MMSEQS2DB should be used as a fallback when MMSEQS_DB_DIR is absent."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.dict(os.environ, {'MMSEQS2DB': tmpdir}, clear=True):
                downloader = MMseqsDownloader()
                assert downloader.get_mmseqsdb_directory() == Path(tmpdir)

    @patch('eplace_lib.ncbi_download.MMseqsDownloader.download_nt_database')
    def test_setup_mmseqs_database(self, mock_download):
        """setup_mmseqs_database delegates to MMseqsDownloader."""
        mock_path = Path("/tmp/mock_mmseqs/NT.20260515/NT")
        mock_download.return_value = (True, "Success", mock_path)

        success, message, db_path = setup_mmseqs_database(force_download=True, threads=4, db_dir=Path("/tmp/mock_mmseqs"))

        assert success is True
        assert message == "Success"
        assert db_path == mock_path
        mock_download.assert_called_once_with(force_download=True, threads=4)

    @patch('eplace_lib.ncbi_download.MMseqsDownloader.add_taxonomy_to_database')
    def test_setup_mmseqs_taxonomy(self, mock_taxonomy):
        """setup_mmseqs_taxonomy delegates to MMseqsDownloader."""
        mock_taxonomy.return_value = (True, "Taxonomy ready")
        mmseqs_db = Path("/tmp/mmseqs/NT.20260515/NT")
        ncbi_taxonomy = Path("/tmp/ncbi_taxonomy")

        success, message = setup_mmseqs_taxonomy(
            mmseqs_db=mmseqs_db,
            ncbi_taxonomy=ncbi_taxonomy,
            threads=8,
            acc2taxid_dir=Path("/tmp/acc2taxid"),
            taxonomy_workdir=Path("/tmp/workdir"),
            db_dir=Path("/tmp/mmseqs")
        )

        assert success is True
        assert message == "Taxonomy ready"
        mock_taxonomy.assert_called_once_with(
            mmseqs_db=mmseqs_db,
            ncbi_taxonomy=ncbi_taxonomy,
            threads=8,
            acc2taxid_dir=Path("/tmp/acc2taxid"),
            taxonomy_workdir=Path("/tmp/workdir")
        )

    @patch('eplace_lib.ncbi_download.get_total_memory_gb')
    def test_check_available_memory_gb(self, mock_total_memory):
        """check_available_memory_gb should compare against requirement."""
        mock_total_memory.return_value = 96.0
        ok, total = check_available_memory_gb(64.0)
        assert ok is True
        assert total == 96.0
