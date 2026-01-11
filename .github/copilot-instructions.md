# ePLACE Copilot Instructions

## Repository Overview

**ePLACE** (environmental Phylogenetic Localisation and Clade Estimation) is a Python bioinformatics library for analyzing environmental DNA (eDNA) sequences through BLAST comparison and taxonomic classification. The repository is approximately 1.1MB with ~9,000 lines of Python code across 20 Python files.

### Key Information
- **Primary Language**: Python 3.8+ (tested with Python 3.12)
- **Package Type**: Python library with CLI interface (installed via pip)
- **Main Module**: `eplace_lib` (located in `src/eplace_lib/`)
- **Distribution**: PyPI package built with setuptools
- **External Dependencies**: BLAST+, TaxonKit, MAFFT, IQTree (installed separately via conda/system packages)
- **Python Dependencies**: pytaxonkit (available only via conda/mamba from bioconda channel)

## Critical Build & Environment Setup

### Installation Order (ALWAYS follow this sequence)

1. **Install External Bioinformatics Tools First** (via conda/mamba):
   ```bash
   mamba install -y bioconda::blast bioconda::pytaxonkit bioconda::iqtree bioconda::mafft
   ```
   **CRITICAL**: pytaxonkit is ONLY available via conda/mamba, NOT via pip. The package will fail to import if pytaxonkit is missing.

2. **Install Python Package**:
   ```bash
   # For development
   pip install -e ".[dev]"
   
   # Basic installation
   pip install .
   ```

3. **Verify Installation**:
   ```bash
   eplace --help  # Should display help without errors
   ```

### Known Installation Issues & Workarounds

**Issue**: `ModuleNotFoundError: No module named 'pytaxonkit'`
- **Cause**: pytaxonkit is not available via pip
- **Solution**: MUST install via conda/mamba: `mamba install -y bioconda::pytaxonkit`
- **Note**: There is a known issue with pytaxonkit installation on Python >=3.12 (see INSTALL.md). The code works with older pytaxonkit versions.

**Issue**: Tests fail to run due to import errors
- **Cause**: pytaxonkit not installed or NCBI taxonomy databases not set up
- **Solution**: Install pytaxonkit first, then set up NCBI taxonomy databases for TaxonKit as described in [taxonkit documentation](https://bioinf.shenwei.me/taxonkit/usage/#before-use)

## Build & Test Commands

### Building the Package
```bash
# Install build tool (if not already installed)
pip install build

# Build distribution packages (creates wheel and source dist)
python -m build
# Result: Creates dist/eplace-0.1.2-py3-none-any.whl and dist/eplace-0.1.2.tar.gz
# Time: ~10-15 seconds
```
**Status**: ✅ Build works without external dependencies installed

### Running Tests
```bash
# PREREQUISITE: pytaxonkit must be installed via conda/mamba
# Install test dependencies
pip install -e ".[dev]"

# Run all tests
pytest tests/ -v

# Run specific test modules
pytest tests/test_blast_analysis.py -v
pytest tests/test_taxonomy.py -v
pytest tests/test_workflow.py -v

# Run with coverage
pytest tests/ --cov=eplace_lib --cov-report=html
```
**Status**: ⚠️ Tests CANNOT run without pytaxonkit installed (conda/mamba required)
**Expected behavior**: All tests use mocking for external dependencies (BLAST, MAFFT, IQTree)

### Building Documentation
```bash
cd docs

# Install documentation dependencies (first time only)
pip install -r requirements.txt
# Requirements: sphinx>=7.0.0, sphinx-rtd-theme>=2.0.0, myst-parser>=2.0.0

# Build HTML documentation
make html
# Time: ~30-45 seconds
# Output: docs/build/html/index.html
# Expected: Succeeds with ~17 warnings about duplicate object descriptions (safe to ignore)
```
**Status**: ✅ Documentation builds successfully (warnings are expected and non-blocking)

### Linting
There are NO linting configurations (no .flake8, .pylintrc, ruff.toml, etc.) in the repository.
**Action**: Do NOT add linting unless specifically requested.

## Project Structure & Architecture

### Repository Layout
```
eplace/
├── .github/
│   └── workflows/
│       └── python-publish.yml     # PyPI publishing workflow (triggered on releases)
├── docs/
│   ├── source/                    # Sphinx documentation source files (.rst, .md)
│   ├── Makefile                   # Documentation build commands
│   └── requirements.txt           # Doc build dependencies
├── examples/
│   ├── blast_workflow_example.py  # Example: individual workflow
│   ├── grouped_workflow_example.py # Example: grouped workflow
│   ├── download_ncbi_example.py   # Example: database download
│   └── relabel_example.py         # Example: tree relabeling
├── slurm/
│   ├── eplace.slurm               # SLURM job script for individual workflow
│   └── eplace_groups.slurm        # SLURM job script for grouped workflow
├── src/eplace_lib/
│   ├── __init__.py                # Package exports and version
│   ├── blast_analysis.py          # BLAST operations (BlastRunner, BlastHit, FastaReader)
│   ├── ncbi_download.py           # NCBI database download (NCBIDownloader)
│   ├── sequences.py               # Sequence utilities (SequenceAnalyzer)
│   ├── taxonomy.py                # Taxonomy extraction (TaxonomyExtractor, SequenceExtractor)
│   ├── alignment.py               # MSA & phylogeny (SequenceTrimmer, MAFFTAligner, IQTreeBuilder)
│   └── cli.py                     # Command-line interface (1101 lines)
├── tests/
│   ├── test_blast_analysis.py     # BlastRunner, FastaReader tests
│   ├── test_ncbi_download.py      # Database download tests
│   ├── test_taxonomy.py           # Taxonomy extraction tests (48KB - largest test file)
│   ├── test_alignment.py          # Alignment and tree building tests (32KB)
│   ├── test_workflow.py           # Integration workflow tests
│   ├── test_cli.py                # CLI command tests
│   └── test_relabel_command.py    # Tree relabeling tests
├── pyproject.toml                 # Package configuration (setuptools-based)
├── requirements.txt               # Python dependencies (pytaxonkit only)
├── README.md                      # Comprehensive usage documentation
├── INSTALL.md                     # Installation instructions
└── citation.cff                   # Citation metadata
```

### Key Source Modules (src/eplace_lib/)

**blast_analysis.py** (481 lines)
- `BlastHit`: Dataclass representing BLAST hit with coverage calculations
- `FastaReader`: Read/validate FASTA files
- `BlastRunner`: Execute blastn searches with filtering by identity/coverage
- External command: `blastn` (from BLAST+ suite)

**ncbi_download.py** (322 lines)
- `NCBIDownloader`: Download and verify NCBI core_nt database
- Uses BLASTDB environment variable or ~/blastdb as default
- Handles MD5 checksum verification and tar extraction

**taxonomy.py** (657 lines)
- `TaxonomyExtractor`: Extract taxonomy info using pytaxonkit
- `SequenceExtractor`: Extract sequences from BLAST databases using blastdbcmd
- External commands: `taxonkit` (via pytaxonkit), `blastdbcmd`

**alignment.py** (1596 lines - largest module)
- `SequenceTrimmer`: Trim sequences to BLAST-aligned regions
- `MAFFTAligner`: Multiple sequence alignment using MAFFT
- `IQTreeBuilder`: Build phylogenetic trees using IQTree
- External commands: `mafft`, `iqtree`
- Functions for both individual and grouped workflows

**cli.py** (1101 lines)
- Three main subcommands: `download`, `blast`, `grouped`
- `blast`: Individual analysis (one tree per query)
- `grouped`: Group queries by taxonomy (one tree per taxonomic group)
- Extensive argument parsing and workflow orchestration

**sequences.py** (68 lines)
- `SequenceAnalyzer`: Basic sequence utilities (GC content, reverse complement)

### Workflow Modes

**Individual Workflow** (`eplace blast`):
- Processes each query sequence independently
- Creates one phylogenetic tree per query sequence
- Output: One directory per query with representatives, alignment, tree

**Grouped Workflow** (`eplace grouped`):
- Groups queries by taxonomic classification (e.g., class, order)
- Creates one tree per taxonomic group
- Removes redundant reference sequences within groups
- Output: Per-query directories + grouped directories with combined analyses

## Testing Strategy

- **Framework**: pytest with pytest-cov
- **Approach**: Heavy use of mocking (`unittest.mock.patch`) for external commands
- **Coverage**: Unit tests for all major modules + integration tests
- **Test Data**: Tests create temporary directories and mock BLAST/taxonomy results
- **No External Dependencies**: Tests do not require BLAST, MAFFT, or IQTree to be installed
- **Exception**: pytaxonkit import is required (conda/mamba installation needed)

## Configuration Files

- **pyproject.toml**: Package metadata, dependencies, entry points
  - `[project.scripts]`: Defines `eplace` command
  - `[project.optional-dependencies]`: Separates dev and docs dependencies
  - No linting or code quality tool configurations present
- **.readthedocs.yaml**: ReadTheDocs build configuration
  - Python 3.12, Sphinx documentation
- **citation.cff**: CITATION.cff format for academic citation

## CI/CD & Validation

### GitHub Actions
- **File**: `.github/workflows/python-publish.yml`
- **Trigger**: On release publication only
- **Steps**:
  1. Build distributions using `python -m build`
  2. Upload to PyPI using trusted publishing
- **No CI Tests**: There is NO automated testing workflow in GitHub Actions
- **No Pre-commit Hooks**: No pre-commit, black, isort, or other automation

### Manual Validation Steps
Since there's no CI, ALWAYS validate changes manually:
1. Install package in clean environment: `pip install -e ".[dev]"`
2. Run tests: `pytest tests/ -v` (requires pytaxonkit)
3. Build package: `python -m build`
4. Build docs: `cd docs && make html`
5. Test CLI: `eplace --help`, `eplace blast --help`, `eplace grouped --help`

## Common Development Patterns

### Code Style Conventions (observed from codebase)
- Google-style docstrings for all public functions/classes
- Type hints on function signatures (Python 3.8+ style)
- Dataclasses for structured data (e.g., BlastHit)
- Logging via Python's logging module (not print statements)
- PEP 8 compliance (though not enforced by tools)

### Error Handling
- Functions return tuples: `(success: bool, result/message)`
- Extensive logging at INFO and WARNING levels
- subprocess calls with explicit error handling and logging

### External Command Execution
- All external tools (blastn, mafft, iqtree, blastdbcmd) are called via `subprocess.run()`
- Commands check for tool availability before execution
- stdout/stderr are captured and logged

## Important Development Notes

1. **pytaxonkit is Required**: The package cannot be imported without pytaxonkit. ALWAYS install via conda/mamba before running tests or using the package.

2. **External Tools are Optional for Core Library**: The package will import successfully without BLAST+, MAFFT, or IQTree installed. These are only required at runtime when calling specific functions.

3. **Large Test Files**: `test_taxonomy.py` (48KB) and `test_alignment.py` (32KB) are the largest test files. They contain extensive taxonomy data and alignment test cases.

4. **No Makefile Commands**: There is no top-level Makefile. Build commands are direct: `python -m build`, `pytest tests/`.

5. **Documentation Warnings are Expected**: Sphinx build produces ~17 warnings about duplicate object descriptions. These are safe to ignore and do not block the build.

6. **SLURM Scripts**: The `slurm/` directory contains example job submission scripts for HPC environments. These are examples, not required for development.

7. **Version Number**: Package version is defined in `pyproject.toml` as `0.1.2`. The `__init__.py` has `__version__ = "0.1.0"` which is outdated but not used by packaging.

## Quick Reference

**Installing for development**:
```bash
mamba create -yn eplace 'python>=3.12'
mamba activate eplace
mamba install -y bioconda::blast bioconda::pytaxonkit bioconda::iqtree bioconda::mafft
pip install -e ".[dev]"
```

**Running the full validation suite**:
```bash
# Build package
python -m build

# Run tests (requires pytaxonkit)
pytest tests/ -v

# Build documentation
cd docs && make html && cd ..

# Test CLI
eplace --help
```

**Making changes**: 
- ALWAYS test imports: `python -c "from eplace_lib import *"`
- If adding new modules, update `__init__.py` exports
- If changing CLI, update both code and documentation
- Use logging instead of print statements for consistency

## Trust These Instructions

These instructions are based on direct inspection and validation of the repository. If you encounter discrepancies:
1. First, verify you followed the installation order (pytaxonkit via conda/mamba first)
2. Check that external tools are installed if needed at runtime
3. Only then perform additional exploration or report the discrepancy
