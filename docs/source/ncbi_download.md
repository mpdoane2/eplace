# NCBI Database Download Module

This module provides functionality for downloading and managing NCBI BLAST databases, specifically the core nucleotide (nt) database.

## Features

- Automatic detection of BLASTDB environment variable
- Fallback to `~/blastdb` directory if BLASTDB is not set
- Downloads core_nt.* files from NCBI FTP server
- Verifies MD5 checksums for data integrity
- Automatic extraction of tar.gz files
- Security features to prevent path traversal attacks
- Configurable logging output

## Quick Start

### Basic Usage

```python
from eplace_lib.ncbi_download import setup_ncbi_database

# Download and setup the database
success, message = setup_ncbi_database()
print(f"Success: {success}")
print(f"Message: {message}")
```

### Check if Database Exists

```python
from eplace_lib.ncbi_download import NCBIDownloader

downloader = NCBIDownloader()
db_dir = downloader.get_blastdb_directory()
exists = downloader.check_database_exists()

print(f"BLASTDB directory: {db_dir}")
print(f"Database exists: {exists}")
```

### Force Redownload

```python
from eplace_lib.ncbi_download import setup_ncbi_database

# Force redownload even if database exists
success, message = setup_ncbi_database(force_download=True)
```

### Disable Verbose Output

```python
from eplace_lib.ncbi_download import setup_ncbi_database

# Download without progress messages
success, message = setup_ncbi_database(verbose=False)
```

## BLASTDB Environment Variable

The module checks for the `BLASTDB` environment variable:

- If set and points to a valid directory, uses that location
- If not set or invalid, creates and uses `~/blastdb`

### Setting BLASTDB

```bash
# In bash
export BLASTDB=/path/to/your/blastdb

# Or in Python
import os
os.environ['BLASTDB'] = '/path/to/your/blastdb'
```

## Advanced Usage

```python
from eplace_lib.ncbi_download import NCBIDownloader

downloader = NCBIDownloader()

# Get database directory
db_dir = downloader.get_blastdb_directory()

# Get list of available files from NCBI
files = downloader.get_available_files()

# Download a specific file
downloader.download_file('core_nt.00.tar.gz', db_dir)

# Verify MD5 checksum
downloader.verify_md5(tar_path, md5_path)

# Extract tarball
downloader.extract_tarball(tar_path, db_dir)
```

## Security

The module includes security features to prevent path traversal attacks:

- Filename validation in `download_file()`
- Path validation before extraction in `extract_tarball()`
- Safe extraction that validates all member paths

## Requirements

Uses Python standard library modules only:
- `os`
- `hashlib`
- `tarfile`
- `logging`
- `pathlib`
- `urllib.request`

No external dependencies required.

## Examples

See `examples/download_ncbi_example.py` for comprehensive usage examples.

## Testing

Run the test suite:

```bash
pytest tests/test_ncbi_download.py -v
```

## Note

The NCBI core_nt database is large (hundreds of GB). Ensure you have sufficient disk space and bandwidth before downloading.
