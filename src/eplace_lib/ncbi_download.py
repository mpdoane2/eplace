"""
NCBI database download module.

This module provides functionality for downloading and managing NCBI BLAST databases,
specifically the core nucleotide (nt) database.
"""

import os
import csv
import gzip
import hashlib
import shutil
import subprocess
import tarfile
from datetime import datetime, timezone
import tempfile
import logging
from pathlib import Path
from typing import Optional, List, Tuple
from urllib.request import urlopen, urlretrieve
from urllib.error import URLError

# Configure module logger
logger = logging.getLogger(__name__)


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
            ValueError: If filename contains path traversal sequences
        """
        # Validate filename to prevent path traversal attacks
        if '..' in filename or filename.startswith('/') or '\\' in filename:
            raise ValueError(f"Invalid filename: {filename}")
        
        url = self.NCBI_FTP_BASE + filename
        dest_path = dest_dir / filename
        
        # Ensure the destination path is within dest_dir (Python 3.9+)
        try:
            dest_path.relative_to(dest_dir)
        except ValueError:
            raise ValueError(f"Invalid destination path: {dest_path}")
        
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
            
        Raises:
            ValueError: If MD5 file format is invalid
        """
        # Read expected MD5 from file
        with open(md5_file_path, 'r') as f:
            md5_content = f.read().strip()
            
        if not md5_content:
            raise ValueError(f"MD5 file is empty: {md5_file_path}")
        
        # MD5 files typically have format: "checksum filename"
        parts = md5_content.split()
        if not parts:
            raise ValueError(f"Invalid MD5 file format: {md5_file_path}")
        
        expected_md5 = parts[0]
        
        # Validate MD5 hash format (32 hex characters)
        if len(expected_md5) != 32 or not all(c in '0123456789abcdefABCDEF' for c in expected_md5):
            raise ValueError(f"Invalid MD5 hash format: {expected_md5}")
        
        # Calculate actual MD5
        md5_hash = hashlib.md5()
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                md5_hash.update(chunk)
        actual_md5 = md5_hash.hexdigest()
        
        return actual_md5.lower() == expected_md5.lower()
    
    def extract_tarball(self, tarball_path: Path, dest_dir: Path) -> None:
        """
        Extract a tar.gz file to the specified directory.
        
        Args:
            tarball_path: Path to the tar.gz file
            dest_dir: Destination directory for extraction
            
        Raises:
            tarfile.TarError: If extraction fails
            ValueError: If tarball contains unsafe paths
        """
        dest_dir_resolved = dest_dir.resolve()
        
        with tarfile.open(tarball_path, 'r:gz') as tar:
            # Validate all member paths before extraction to prevent path traversal
            for member in tar.getmembers():
                member_path = (dest_dir / member.name).resolve()
                # Use relative_to() for safe path validation
                try:
                    member_path.relative_to(dest_dir_resolved)
                except ValueError:
                    raise ValueError(f"Unsafe path in tarball: {member.name}")
            
            # Safe to extract after validation
            tar.extractall(path=dest_dir)
    
    def download_and_setup_database(self, force_download: bool = False, verbose: bool = True) -> Tuple[bool, str]:
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
            verbose: If True, logs progress information (default: True)
            
        Returns:
            Tuple of (success: bool, message: str)
        """
        # Set up logging level
        if verbose:
            logging.basicConfig(level=logging.INFO, format='%(message)s')
        
        db_dir = self.get_blastdb_directory()
        
        # Check if database already exists
        if not force_download and self.check_database_exists(db_dir):
            msg = f"Database already exists in {db_dir}"
            logger.info(msg)
            return True, msg
        
        try:
            # Get list of available files
            logger.info("Fetching list of available files from NCBI FTP...")
            files = self.get_available_files()
            
            if not files:
                return False, "No core_nt files found on NCBI FTP server"
            
            # Separate tar.gz and md5 files
            tar_files = [f for f in files if f.endswith('.tar.gz') and not f.endswith('.md5')]
            md5_files = [f for f in files if f.endswith('.tar.gz.md5')]
            
            # Download all files
            downloaded_tarballs = []
            for tar_file in tar_files:
                logger.info(f"Downloading {tar_file}...")
                tar_path = self.download_file(tar_file, db_dir)
                
                # Download corresponding MD5 file
                md5_file = f"{tar_file}.md5"
                if md5_file in md5_files:
                    logger.info(f"Downloading {md5_file}...")
                    md5_path = self.download_file(md5_file, db_dir)
                    
                    # Verify MD5
                    logger.info(f"Verifying checksum for {tar_file}...")
                    if not self.verify_md5(tar_path, md5_path):
                        return False, f"MD5 checksum verification failed for {tar_file}"
                    
                    logger.info(f"Checksum verified for {tar_file}")
                else:
                    logger.warning(f"Warning: No MD5 file found for {tar_file}")
                
                downloaded_tarballs.append(tar_path)
            
            # Extract all tarballs
            for tar_path in downloaded_tarballs:
                logger.info(f"Extracting {tar_path.name}...")
                self.extract_tarball(tar_path, db_dir)
                logger.info(f"Extracted {tar_path.name}")
            
            msg = f"Successfully downloaded and extracted database to {db_dir}"
            logger.info(msg)
            return True, msg
            
        except URLError as e:
            return False, f"Network error: {e}"
        except Exception as e:
            return False, f"Error during download/extraction: {e}"


def setup_ncbi_database(force_download: bool = False, verbose: bool = True) -> Tuple[bool, str]:
    """
    Convenience function to setup the NCBI core_nt database.
    
    Args:
        force_download: If True, downloads even if database exists
        verbose: If True, logs progress information (default: True)
        
    Returns:
        Tuple of (success: bool, message: str)
    """
    downloader = NCBIDownloader()
    return downloader.download_and_setup_database(force_download, verbose)


def get_total_memory_gb() -> float:
    """Get total system memory in GiB.

    On Linux this first reads ``/proc/meminfo`` (``MemTotal``). If that is not
    available, it falls back to POSIX ``os.sysconf``. Returns ``0.0`` when both
    strategies fail.
    """
    try:
        with open("/proc/meminfo", "r") as meminfo:
            for line in meminfo:
                if line.startswith("MemTotal:"):
                    parts = line.split()
                    if len(parts) >= 2:
                        total_kb = int(parts[1])
                        return total_kb / (1024 ** 2)
    except (OSError, ValueError):
        pass

    try:
        page_size = os.sysconf("SC_PAGE_SIZE")
        num_pages = os.sysconf("SC_PHYS_PAGES")
        return (page_size * num_pages) / (1024 ** 3)
    except (AttributeError, ValueError, OSError):
        return 0.0


def check_available_memory_gb(required_gb: float) -> Tuple[bool, float]:
    """Check whether total system memory meets a required threshold in GiB."""
    total_gb = get_total_memory_gb()
    return total_gb >= required_gb, total_gb


class MMseqsDownloader:
    """Download and configure MMseqs2 NT databases and taxonomy sidecar files.

    Directory resolution for MMseqs2 databases prefers ``$MMSEQS_DB_DIR``, then
    ``$MMSEQS2DB`` (legacy), then ``~/mmseqs2db``.

    Workflow:
    1. Download NT with ``mmseqs databases``.
    2. Optionally fetch accession2taxid files and build mapping TSV.
    3. Attach taxonomy sidecars with ``mmseqs createtaxdb``.
    """

    MMSEQS_DB_DIR_ENV = "MMSEQS_DB_DIR"
    LEGACY_MMSEQS_DB_DIR_ENV = "MMSEQS2DB"
    NUCLEOTIDE_DB_NAME = "NT"
    ACC2TAXID_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/"

    def __init__(self, db_dir: Optional[Path] = None):
        """Initialize the MMseqs downloader."""
        self.db_dir = db_dir

    def get_mmseqsdb_directory(self) -> Path:
        """Get or determine the MMseqs2 database directory.

        Resolution order:
        1. Explicit path passed at initialization.
        2. ``$MMSEQS_DB_DIR``.
        3. ``$MMSEQS2DB`` (legacy fallback).
        4. ``~/mmseqs2db``.
        """
        if self.db_dir is None:
            mmseqs_env = (
                os.environ.get(self.MMSEQS_DB_DIR_ENV)
                or os.environ.get(self.LEGACY_MMSEQS_DB_DIR_ENV)
            )
            if mmseqs_env:
                self.db_dir = Path(mmseqs_env)
            else:
                self.db_dir = Path.home() / "mmseqs2db"

        self.db_dir.mkdir(parents=True, exist_ok=True)
        return self.db_dir

    @staticmethod
    def _run_command(cmd: list[str], error_context: str) -> Tuple[bool, str]:
        """Run a command and return success plus an error message on failure."""
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
        except FileNotFoundError:
            executable = cmd[0] if cmd else "command"
            return False, f"{executable} was not found. Please install it and ensure it is in PATH."
        except subprocess.SubprocessError as e:
            return False, f"{error_context}: {e}"

        if result.returncode != 0:
            stderr = result.stderr.strip() if result.stderr else "unknown error"
            return False, f"{error_context}: {stderr}"

        return True, ""

    def download_nt_database(
        self,
        force_download: bool = False,
        threads: int = 1
    ) -> Tuple[bool, str, Optional[Path]]:
        """Download MMseqs2 NT database using ``mmseqs databases``."""
        mmseqs_dir = self.get_mmseqsdb_directory()
        date_stamp = datetime.now(timezone.utc).strftime("%Y%m%d")
        target_dir = mmseqs_dir / f"{self.NUCLEOTIDE_DB_NAME}.{date_stamp}"
        target_db = target_dir / self.NUCLEOTIDE_DB_NAME

        if target_db.exists() and not force_download:
            return True, f"MMseqs2 NT database already exists at {target_db}", target_db

        if force_download and target_dir.exists():
            shutil.rmtree(target_dir)

        target_dir.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory(prefix="mmseqs_tmp_", dir=str(mmseqs_dir)) as tmp_dir:
            cmd = [
                "mmseqs",
                "databases",
                "--threads",
                str(threads),
                self.NUCLEOTIDE_DB_NAME,
                str(target_db),
                tmp_dir
            ]
            success, error = self._run_command(cmd, "MMseqs2 database download failed")
            if not success:
                return False, error, None

        return True, f"MMseqs2 NT database downloaded to {target_db}", target_db

    @staticmethod
    def _validate_mmseqs_database(mmseqs_db: Path) -> Tuple[bool, str]:
        """Validate that MMseqs2 NT database files exist."""
        required_files = [
            mmseqs_db,
            mmseqs_db.with_suffix(".index"),
            mmseqs_db.with_suffix(".dbtype"),
            mmseqs_db.with_suffix(".lookup"),
        ]
        missing = [str(path) for path in required_files if not path.exists()]
        if missing:
            return False, f"Missing required MMseqs2 database files: {', '.join(missing)}"
        return True, ""

    @staticmethod
    def _validate_taxonomy_dump(ncbi_taxonomy: Path) -> Tuple[bool, str]:
        """Validate taxonomy dump files required by ``mmseqs createtaxdb``."""
        required_files = ["nodes.dmp", "names.dmp", "merged.dmp"]
        missing = [f for f in required_files if not (ncbi_taxonomy / f).is_file()]
        if missing:
            return False, f"Missing required NCBI taxonomy files in {ncbi_taxonomy}: {', '.join(missing)}"
        return True, ""

    def _resolve_acc2taxid_dir(self, ncbi_taxonomy: Path, acc2taxid_dir: Optional[Path]) -> Path:
        """Resolve accession2taxid directory from arg, env, or taxonomy path.

        Resolution order:
        1. Explicit ``acc2taxid_dir`` argument.
        2. ``$ACC2TAXID_DIR`` environment variable.
        3. ``<ncbi_taxonomy>/accession2taxid``.
        """
        if acc2taxid_dir is not None:
            return acc2taxid_dir

        env_value = os.environ.get("ACC2TAXID_DIR")
        if env_value:
            return Path(env_value)

        return ncbi_taxonomy / "accession2taxid"

    def _ensure_accession2taxid_files(self, acc2taxid_dir: Path) -> Tuple[bool, str]:
        """Download required accession2taxid files if they are missing."""
        acc2taxid_dir.mkdir(parents=True, exist_ok=True)
        required_files = [
            "nucl_gb.accession2taxid.gz",
            "nucl_wgs.accession2taxid.gz",
        ]

        for filename in required_files:
            target = acc2taxid_dir / filename
            if target.exists() and target.stat().st_size > 0:
                continue
            try:
                urlretrieve(self.ACC2TAXID_BASE + filename, target)
            except URLError as e:
                return False, f"Failed to download {filename}: {e}"
        return True, ""

    @staticmethod
    def _build_tax_mapping_file(acc2taxid_dir: Path, mapping_file: Path) -> Tuple[bool, str]:
        """Create accession-to-taxid mapping file for MMseqs2 taxonomy.

        The output is a two-column TSV: ``accession<TAB>taxid``. Both versioned
        (``accession.version``) and unversioned accessions are emitted, then
        deduplicated via ``sort -u``.
        """
        unsorted_mapping = mapping_file.with_suffix(".unsorted.tsv")
        try:
            with open(unsorted_mapping, "w") as out:
                for gz_file in sorted(acc2taxid_dir.glob("nucl_*.accession2taxid.gz")):
                    with gzip.open(gz_file, "rt", newline="") as infile:
                        reader = csv.reader(infile, delimiter="\t")
                        for row in reader:
                            if len(row) < 3:
                                continue
                            if row[0] == "accession":
                                continue
                            accession, accession_version, taxid = row[0], row[1], row[2]
                            if not taxid.isdigit():
                                continue
                            # MMseqs2 lookups can use either versioned or
                            # unversioned accessions depending on how the
                            # database was created; emit both for coverage.
                            out.write(f"{accession_version}\t{taxid}\n")
                            out.write(f"{accession}\t{taxid}\n")

            result = subprocess.run(
                ["sort", "-u", str(unsorted_mapping), "-o", str(mapping_file)],
                capture_output=True,
                text=True
            )
            if result.returncode != 0:
                stderr = result.stderr.strip() if result.stderr else "unknown error"
                return False, f"Failed to sort mapping file: {stderr}"
        except (OSError, ValueError, gzip.BadGzipFile) as e:
            return False, f"Failed to create taxonomy mapping file: {e}"
        finally:
            if unsorted_mapping.exists():
                unsorted_mapping.unlink()

        return True, ""

    def add_taxonomy_to_database(
        self,
        mmseqs_db: Path,
        ncbi_taxonomy: Path,
        threads: int = 1,
        acc2taxid_dir: Optional[Path] = None,
        taxonomy_workdir: Optional[Path] = None
    ) -> Tuple[bool, str]:
        """Add NCBI taxonomy sidecar files to an MMseqs2 NT database.

        Args:
            mmseqs_db: Path to MMseqs2 NT database base file (e.g. ``.../NT``).
            ncbi_taxonomy: Directory with NCBI taxonomy dump files.
            threads: Number of threads for ``mmseqs createtaxdb``.
            acc2taxid_dir: Optional directory with accession2taxid files.
            taxonomy_workdir: Optional working directory for mapping files.

        Returns:
            Tuple ``(success, message)``.

        Side effects:
            Creates mapping files in ``taxonomy_workdir`` and writes MMseqs2
            taxonomy sidecar files adjacent to ``mmseqs_db``.
        """
        db_ok, db_error = self._validate_mmseqs_database(mmseqs_db)
        if not db_ok:
            return False, db_error

        taxonomy_ok, taxonomy_error = self._validate_taxonomy_dump(ncbi_taxonomy)
        if not taxonomy_ok:
            return False, taxonomy_error

        resolved_acc2taxid = self._resolve_acc2taxid_dir(ncbi_taxonomy, acc2taxid_dir)
        ensure_ok, ensure_error = self._ensure_accession2taxid_files(resolved_acc2taxid)
        if not ensure_ok:
            return False, ensure_error

        workdir = taxonomy_workdir if taxonomy_workdir is not None else (mmseqs_db.parent / "taxonomy_build")
        workdir.mkdir(parents=True, exist_ok=True)
        mapping_file = workdir / "nt.accession2taxid.mmseqs.tsv"
        map_ok, map_error = self._build_tax_mapping_file(resolved_acc2taxid, mapping_file)
        if not map_ok:
            return False, map_error

        with tempfile.TemporaryDirectory(prefix="mmseqs_tax_tmp_", dir=str(workdir)) as tmp_dir:
            cmd = [
                "mmseqs",
                "createtaxdb",
                str(mmseqs_db),
                tmp_dir,
                "--ncbi-tax-dump",
                str(ncbi_taxonomy),
                "--tax-mapping-file",
                str(mapping_file),
                "--threads",
                str(threads),
            ]
            success, error = self._run_command(cmd, "MMseqs2 taxonomy creation failed")
            if not success:
                return False, error

        return True, f"MMseqs2 taxonomy files created for database {mmseqs_db}"


def setup_mmseqs_database(
    force_download: bool = False,
    threads: int = 1,
    db_dir: Optional[Path] = None
) -> Tuple[bool, str, Optional[Path]]:
    """Convenience function to download MMseqs2 NT database."""
    downloader = MMseqsDownloader(db_dir=db_dir)
    return downloader.download_nt_database(force_download=force_download, threads=threads)


def setup_mmseqs_taxonomy(
    mmseqs_db: Path,
    ncbi_taxonomy: Path,
    threads: int = 1,
    acc2taxid_dir: Optional[Path] = None,
    taxonomy_workdir: Optional[Path] = None,
    db_dir: Optional[Path] = None
) -> Tuple[bool, str]:
    """Convenience function to add taxonomy to an MMseqs2 database."""
    downloader = MMseqsDownloader(db_dir=db_dir)
    return downloader.add_taxonomy_to_database(
        mmseqs_db=mmseqs_db,
        ncbi_taxonomy=ncbi_taxonomy,
        threads=threads,
        acc2taxid_dir=acc2taxid_dir,
        taxonomy_workdir=taxonomy_workdir
    )
