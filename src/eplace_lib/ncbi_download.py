"""
NCBI database download module.

This module provides functionality for downloading and managing NCBI BLAST databases,
specifically the core nucleotide (nt) database.
"""

import os
import hashlib
import tarfile
from pathlib import Path
from typing import Optional, List, Tuple
from urllib.request import urlopen, urlretrieve
from urllib.error import URLError


class NCBIDownloader:
    """
    A class for managing NCBI BLAST database downloads.
    
    This class handles checking for existing databases, downloading from NCBI FTP,
    verifying checksums, and extracting database files.
    """
    
    NCBI_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/blast/db/"
    CORE_NT_PREFIX = "core_nt"
    
    def __init__(self):
        """Initialize the NCBIDownloader."""
        self.blastdb_dir = None
    
    def get_blastdb_directory(self) -> Path:
        """
        Get or determine the BLASTDB directory.
        
        Checks if the BLASTDB environment variable is set. If it exists and points
        to a valid directory, uses that. Otherwise, creates and returns a path to
        ~/blastdb.
        
        Returns:
            Path object pointing to the BLASTDB directory
        """
        if self.blastdb_dir is not None:
            return self.blastdb_dir
        
        # Check BLASTDB environment variable
        blastdb_env = os.environ.get('BLASTDB')
        
        if blastdb_env:
            blastdb_path = Path(blastdb_env)
            if blastdb_path.exists() and blastdb_path.is_dir():
                self.blastdb_dir = blastdb_path
                return self.blastdb_dir
        
        # Use ~/blastdb as default
        home_dir = Path.home()
        blastdb_path = home_dir / "blastdb"
        
        # Create directory if it doesn't exist
        blastdb_path.mkdir(parents=True, exist_ok=True)
        
        self.blastdb_dir = blastdb_path
        return self.blastdb_dir
    
    def check_database_exists(self, db_dir: Optional[Path] = None) -> bool:
        """
        Check if NCBI core_nt database files exist in the specified directory.
        
        Args:
            db_dir: Directory to check. If None, uses the default BLASTDB directory.
            
        Returns:
            True if at least one core_nt database file exists, False otherwise
        """
        if db_dir is None:
            db_dir = self.get_blastdb_directory()
        
        # Look for core_nt.* files
        core_nt_files = list(db_dir.glob(f"{self.CORE_NT_PREFIX}.*"))
        
        # Check for actual database files (not just tarballs)
        db_extensions = ['.nhr', '.nin', '.nsq', '.ndb', '.not', '.ntf', '.nto']
        for file in core_nt_files:
            if any(str(file).endswith(ext) for ext in db_extensions):
                return True
        
        return False
    
    def get_available_files(self) -> List[str]:
        """
        Get list of available core_nt files from NCBI FTP server.
        
        Returns:
            List of filenames matching core_nt pattern
            
        Raises:
            URLError: If unable to connect to FTP server
        """
        try:
            with urlopen(self.NCBI_FTP_BASE) as response:
                html_content = response.read().decode('utf-8')
            
            # Parse HTML to find core_nt files
            files = []
            for line in html_content.split('\n'):
                if f'href="{self.CORE_NT_PREFIX}.' in line:
                    # Extract filename from href
                    start = line.find(f'href="{self.CORE_NT_PREFIX}.')
                    if start != -1:
                        start += 6  # Length of 'href="'
                        end = line.find('"', start)
                        if end != -1:
                            filename = line[start:end]
                            # Only include tar.gz files and their md5 files
                            if filename.endswith('.tar.gz') or filename.endswith('.tar.gz.md5'):
                                files.append(filename)
            
            return sorted(set(files))
        except URLError as e:
            raise URLError(f"Failed to fetch file list from NCBI FTP: {e}")
    
    def download_file(self, filename: str, dest_dir: Path, show_progress: bool = True) -> Path:
        """
        Download a file from NCBI FTP server.
        
        Args:
            filename: Name of the file to download
            dest_dir: Destination directory
            show_progress: Whether to show download progress (not implemented yet)
            
        Returns:
            Path to the downloaded file
            
        Raises:
            URLError: If download fails
        """
        url = self.NCBI_FTP_BASE + filename
        dest_path = dest_dir / filename
        
        try:
            urlretrieve(url, dest_path)
            return dest_path
        except URLError as e:
            raise URLError(f"Failed to download {filename}: {e}")
    
    def verify_md5(self, file_path: Path, md5_file_path: Path) -> bool:
        """
        Verify the MD5 checksum of a file.
        
        Args:
            file_path: Path to the file to verify
            md5_file_path: Path to the MD5 checksum file
            
        Returns:
            True if checksum matches, False otherwise
        """
        # Read expected MD5 from file
        with open(md5_file_path, 'r') as f:
            md5_content = f.read().strip()
            # MD5 files typically have format: "checksum filename"
            expected_md5 = md5_content.split()[0]
        
        # Calculate actual MD5
        md5_hash = hashlib.md5()
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                md5_hash.update(chunk)
        actual_md5 = md5_hash.hexdigest()
        
        return actual_md5 == expected_md5
    
    def extract_tarball(self, tarball_path: Path, dest_dir: Path) -> None:
        """
        Extract a tar.gz file to the specified directory.
        
        Args:
            tarball_path: Path to the tar.gz file
            dest_dir: Destination directory for extraction
            
        Raises:
            tarfile.TarError: If extraction fails
        """
        with tarfile.open(tarball_path, 'r:gz') as tar:
            tar.extractall(path=dest_dir)
    
    def download_and_setup_database(self, force_download: bool = False) -> Tuple[bool, str]:
        """
        Main function to download and setup the NCBI core_nt database.
        
        This function:
        1. Determines the BLASTDB directory
        2. Checks if database already exists (unless force_download is True)
        3. Downloads all core_nt.* files from NCBI FTP
        4. Verifies MD5 checksums
        5. Extracts the database files
        
        Args:
            force_download: If True, downloads even if database exists
            
        Returns:
            Tuple of (success: bool, message: str)
        """
        db_dir = self.get_blastdb_directory()
        
        # Check if database already exists
        if not force_download and self.check_database_exists(db_dir):
            return True, f"Database already exists in {db_dir}"
        
        try:
            # Get list of available files
            files = self.get_available_files()
            
            if not files:
                return False, "No core_nt files found on NCBI FTP server"
            
            # Separate tar.gz and md5 files
            tar_files = [f for f in files if f.endswith('.tar.gz') and not f.endswith('.md5')]
            md5_files = [f for f in files if f.endswith('.tar.gz.md5')]
            
            # Download all files
            downloaded_tarballs = []
            for tar_file in tar_files:
                print(f"Downloading {tar_file}...")
                tar_path = self.download_file(tar_file, db_dir)
                
                # Download corresponding MD5 file
                md5_file = f"{tar_file}.md5"
                if md5_file in md5_files:
                    print(f"Downloading {md5_file}...")
                    md5_path = self.download_file(md5_file, db_dir)
                    
                    # Verify MD5
                    print(f"Verifying checksum for {tar_file}...")
                    if not self.verify_md5(tar_path, md5_path):
                        return False, f"MD5 checksum verification failed for {tar_file}"
                    
                    print(f"Checksum verified for {tar_file}")
                else:
                    print(f"Warning: No MD5 file found for {tar_file}")
                
                downloaded_tarballs.append(tar_path)
            
            # Extract all tarballs
            for tar_path in downloaded_tarballs:
                print(f"Extracting {tar_path.name}...")
                self.extract_tarball(tar_path, db_dir)
                print(f"Extracted {tar_path.name}")
            
            return True, f"Successfully downloaded and extracted database to {db_dir}"
            
        except URLError as e:
            return False, f"Network error: {e}"
        except Exception as e:
            return False, f"Error during download/extraction: {e}"


def setup_ncbi_database(force_download: bool = False) -> Tuple[bool, str]:
    """
    Convenience function to setup the NCBI core_nt database.
    
    Args:
        force_download: If True, downloads even if database exists
        
    Returns:
        Tuple of (success: bool, message: str)
    """
    downloader = NCBIDownloader()
    return downloader.download_and_setup_database(force_download)
