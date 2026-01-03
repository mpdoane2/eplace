#!/usr/bin/env python
"""
Example script demonstrating how to use the NCBI database download functionality.

This script shows how to:
1. Check if the NCBI database exists
2. Download and setup the database
3. Use custom BLASTDB locations
"""

import os
from pathlib import Path
from eplace_lib.ncbi_download import NCBIDownloader, setup_ncbi_database


def example_basic_usage():
    """Basic usage example."""
    print("=" * 60)
    print("Example 1: Basic Usage")
    print("=" * 60)
    
    # Simple one-line setup
    # Note: This will actually download large files in a real scenario
    success, message = setup_ncbi_database()
    print(f"Success: {success}")
    print(f"Message: {message}")
    
    print("To download the database, call:")
    print("  success, message = setup_ncbi_database()")
    print("\nThis will:")
    print("  - Check BLASTDB environment variable")
    print("  - Use ~/blastdb if BLASTDB is not set")
    print("  - Download core_nt.* files from NCBI")
    print("  - Verify MD5 checksums")
    print("  - Extract the files")
    print()


def example_check_database():
    """Example of checking if database exists."""
    print("=" * 60)
    print("Example 2: Check if Database Exists")
    print("=" * 60)
    
    downloader = NCBIDownloader()
    db_dir = downloader.get_blastdb_directory()
    exists = downloader.check_database_exists()
    
    print(f"BLASTDB directory: {db_dir}")
    print(f"Database exists: {exists}")
    print()


def example_custom_location():
    """Example using custom BLASTDB location."""
    print("=" * 60)
    print("Example 3: Using Custom BLASTDB Location")
    print("=" * 60)
    
    # Set BLASTDB environment variable
    custom_path = Path.home() / "my_custom_blastdb"
    os.environ['BLASTDB'] = str(custom_path)
    
    downloader = NCBIDownloader()
    db_dir = downloader.get_blastdb_directory()
    
    print(f"Custom BLASTDB directory: {db_dir}")
    print(f"Directory created: {db_dir.exists()}")
    
    # Clean up - remove from environment
    del os.environ['BLASTDB']
    print()


def example_force_redownload():
    """Example of forcing a redownload."""
    print("=" * 60)
    print("Example 4: Force Redownload")
    print("=" * 60)
    
    print("To force redownload even if database exists:")
    print("  success, message = setup_ncbi_database(force_download=True)")
    print("\nThis is useful if:")
    print("  - Database files are corrupted")
    print("  - You want to update to the latest version")
    print()


def example_advanced_usage():
    """Example of advanced usage with the class."""
    print("=" * 60)
    print("Example 5: Advanced Usage")
    print("=" * 60)
    
    downloader = NCBIDownloader()
    
    # Get database directory
    db_dir = downloader.get_blastdb_directory()
    print(f"Database directory: {db_dir}")
    
    # Check if database exists
    exists = downloader.check_database_exists(db_dir)
    print(f"Database exists: {exists}")
    
    # In a real scenario, you would:
    # 1. Get available files
    # files = downloader.get_available_files()
    # 
    # 2. Download files
    # for file in files:
    #     if file.endswith('.tar.gz'):
    #         downloader.download_file(file, db_dir)
    #
    # 3. Verify checksums
    # downloader.verify_md5(tar_file, md5_file)
    #
    # 4. Extract files
    # downloader.extract_tarball(tar_file, db_dir)
    
    print("\nFor full control, use the NCBIDownloader class methods directly.")
    print()


if __name__ == "__main__":
    print("\nNCBI Database Download Examples")
    print("=" * 60)
    print("\nNOTE: Actual downloads are commented out to avoid")
    print("downloading large files during testing.")
    print("=" * 60)
    print()
    
    example_basic_usage()
    example_check_database()
    example_custom_location()
    example_force_redownload()
    example_advanced_usage()
    
    print("=" * 60)
    print("Examples completed!")
    print("=" * 60)
