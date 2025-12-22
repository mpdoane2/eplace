# eplace

ePLACE: environmental Phylogenetic Localisation and Clade Estimation

A Python library for analyzing environmental DNA (eDNA) sequences through BLAST comparison and taxonomic classification.

## Features

- **NCBI Database Management**: Download and manage NCBI BLAST databases (core_nt)
- **FASTA File Processing**: Read and validate FASTA files
- **BLAST Search**: Run blastn searches with configurable parameters
- **Result Filtering**: Filter BLAST results by identity and coverage thresholds
- **Taxonomic Analysis**: Extract representative sequences per taxonomic rank
- **Sequence Extraction**: Retrieve sequences from BLAST databases

## Installation

```bash
# Clone the repository
git clone https://github.com/linsalrob/eplace.git
cd eplace

# Install in development mode
pip install -e .

# Or with development dependencies
pip install -e ".[dev]"
```

## Requirements

- Python 3.8 or higher
- BLAST+ tools (blastn, blastdbcmd) must be installed separately:
  ```bash
  # Ubuntu/Debian
  sudo apt-get install ncbi-blast+
  
  # macOS with Homebrew
  brew install blast
  ```

## Quick Start

### 1. Download NCBI Database

```python
from eplace_lib import setup_ncbi_database

# Download the core_nt database
success, message = setup_ncbi_database()
```

### 2. Run BLAST Analysis

```python
from pathlib import Path
from eplace_lib import run_blast_search, process_blast_results_for_taxonomy

# Run BLAST search with filtering
success, filtered_hits = run_blast_search(
    query_fasta=Path("query.fasta"),
    output_file=Path("blast_results.txt"),
    min_identity=90.0,    # 90% identity threshold
    min_coverage=80.0     # 80% query coverage threshold
)

# Extract representative sequences by taxonomic rank
results = process_blast_results_for_taxonomy(
    blast_hits=filtered_hits,
    output_dir=Path("output"),
    rank="species"  # Options: phylum, class, order, family, genus, species
)
```

### 3. Command-Line Interface

```bash
# Basic usage
python examples/blast_workflow_example.py query.fasta output_dir

# With custom parameters
python examples/blast_workflow_example.py query.fasta output_dir \
    --rank genus \
    --min-identity 95 \
    --min-coverage 85 \
    --num-threads 4

# Show help
python examples/blast_workflow_example.py --help
```

## Documentation

- [NCBI Database Download](docs/ncbi_download.md) - Downloading and managing BLAST databases
- [BLAST Workflow](docs/blast_workflow.md) - Complete guide to BLAST analysis and taxonomy extraction

## Examples

See the `examples/` directory for comprehensive examples:

- `download_ncbi_example.py` - Download and manage NCBI databases
- `blast_workflow_example.py` - Complete BLAST workflow with command-line interface

## Testing

Run the test suite:

```bash
# Run all tests
pytest tests/ -v

# Run specific test modules
pytest tests/test_blast_analysis.py -v
pytest tests/test_taxonomy.py -v
pytest tests/test_workflow.py -v

# Run with coverage
pytest tests/ --cov=eplace_lib --cov-report=html
```

## Project Structure

```
eplace/
├── src/
│   └── eplace_lib/
│       ├── __init__.py
│       ├── blast_analysis.py    # BLAST operations
│       ├── ncbi_download.py     # Database management
│       ├── sequences.py         # Sequence analysis utilities
│       └── taxonomy.py          # Taxonomy extraction
├── tests/
│   ├── test_blast_analysis.py
│   ├── test_ncbi_download.py
│   ├── test_taxonomy.py
│   └── test_workflow.py
├── examples/
│   ├── blast_workflow_example.py
│   └── download_ncbi_example.py
├── docs/
│   ├── blast_workflow.md
│   └── ncbi_download.md
└── pyproject.toml
```

## Workflow Overview

1. **Download Database**: Use `setup_ncbi_database()` to download NCBI core_nt database
2. **Prepare Query**: Create a FASTA file with your query sequences
3. **Run BLAST**: Use `run_blast_search()` to search against the database
4. **Filter Results**: Automatically filter by identity and coverage thresholds
5. **Extract Representatives**: Select representative sequences per taxonomic rank
6. **Output**: Get FASTA files with representative sequences (one per query)

## Output Structure

```
output_dir/
├── blast_results.txt              # Raw BLAST results
├── query1_id/
│   └── query1_id_representatives.fasta
├── query2_id/
│   └── query2_id_representatives.fasta
└── ...
```

## License

MIT License - See LICENSE file for details

## Authors

- Rob Edwards (raedwards@gmail.com)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use ePLACE in your research, please cite:

```
Edwards, R. (2024). ePLACE: environmental Phylogenetic Localisation and Clade Estimation.
GitHub repository: https://github.com/linsalrob/eplace
```

## Support

For issues, questions, or suggestions, please open an issue on GitHub:
https://github.com/linsalrob/eplace/issues
