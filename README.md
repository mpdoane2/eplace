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
- **Sequence Trimming**: Trim reference sequences to aligned regions based on BLAST coordinates
- **Multiple Sequence Alignment**: Align sequences using MAFFT with auto-orientation
- **Phylogenetic Trees**: Build and label phylogenetic trees using IQTree

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
- TaxonKit (for taxonomy lookup):
  ```bash
  # Download from: https://github.com/shenwei356/taxonkit/releases
  # Or install via conda:
  conda install -c bioconda taxonkit
  ```
- MAFFT (optional, for sequence alignment):
  ```bash
  # Ubuntu/Debian
  sudo apt-get install mafft
  
  # macOS with Homebrew
  brew install mafft
  
  # Or via conda:
  conda install -c bioconda mafft
  ```
- IQTree (optional, for phylogenetic tree building):
  ```bash
  # Ubuntu/Debian
  sudo apt-get install iqtree
  
  # macOS with Homebrew
  brew install iqtree
  
  # Or via conda:
  conda install -c bioconda iqtree
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
# Basic usage - BLAST search and representative sequence extraction
python examples/blast_workflow_example.py query.fasta output_dir

# With custom parameters
python examples/blast_workflow_example.py query.fasta output_dir \
    --rank genus \
    --min-identity 95 \
    --min-coverage 85 \
    --num-threads 4

# Skip alignment and tree building (BLAST and extraction only)
python examples/blast_workflow_example.py query.fasta output_dir \
    --skip-alignment

# Show help
python examples/blast_workflow_example.py --help
```

## Documentation

- [NCBI Database Download](docs/ncbi_download.md) - Downloading and managing BLAST databases
- [BLAST Workflow](docs/blast_workflow.md) - Complete guide to BLAST analysis and taxonomy extraction

## Examples

See the `examples/` directory for comprehensive examples:

- `download_ncbi_example.py` - Download and manage NCBI databases
- `blast_workflow_example.py` - Complete BLAST workflow with per-query analysis and trees
- `grouped_workflow_example.py` - Grouped workflow that combines queries by taxonomic rank

### Grouped Workflow

The grouped workflow (`grouped_workflow_example.py`) differs from the standard workflow by:
- Grouping all queries that match to the same taxonomic rank (e.g., class, order)
- Creating one FASTA file per group containing all queries and unique reference sequences
- Removing redundant reference sequences within each group
- Building one alignment and phylogenetic tree per group (instead of per query)

This is useful when you want to analyze multiple related queries together in a single phylogenetic context.

```bash
# Group queries by class (default)
python examples/grouped_workflow_example.py query.fasta output_dir

# Group by a different taxonomic rank
python examples/grouped_workflow_example.py query.fasta output_dir --group-rank order

# Specify both representative rank and grouping rank
python examples/grouped_workflow_example.py query.fasta output_dir --rank genus --group-rank family
```

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
6. **Trim Sequences**: Extract aligned regions from reference sequences based on BLAST coordinates
7. **Align Sequences**: Use MAFFT to align query with trimmed reference sequences (optional)
8. **Build Tree**: Build phylogenetic tree using IQTree with taxonomic labels (optional)
9. **Output**: Get FASTA files, alignments, and trees (one set per query)

### Grouped Workflow Overview

The grouped workflow adds an additional step:
1-5. Same as standard workflow through representative extraction
6. **Group by Rank**: Group all queries by specified taxonomic rank (e.g., class)
7. **Create Grouped FASTA**: Combine all queries and unique references for each group
8. **Trim Sequences**: Trim references to aligned regions
9. **Check Consistency**: Verify BLAST hits align to similar locations on references
10. **Align and Build Trees**: Create one alignment and tree per taxonomic group

## Output Structure

### Standard Workflow Output

```
output_dir/
├── blast_results.txt              # Raw BLAST results
├── blast_results_annotated.txt    # BLAST results with taxonomic annotations
├── query1_id/
│   ├── query1_id_representatives.fasta          # Representative sequences
│   ├── query1_id_with_query.fasta              # Query + representatives
│   ├── query1_id_trimmed.fasta                 # Trimmed to aligned regions
│   ├── query1_id_aligned.fasta                 # Multiple sequence alignment
│   ├── query1_id_tree.treefile                 # Phylogenetic tree
│   ├── query1_id_tree_labeled.treefile         # Tree with taxonomic labels
│   └── query1_id_tree.* (other IQTree files)
├── query2_id/
│   └── query2_id_representatives.fasta
└── ...
```

### Grouped Workflow Output

```
output_dir/
├── blast_results.txt              # Raw BLAST results
├── blast_results_annotated.txt    # BLAST results with taxonomic annotations
├── query1_id/                     # Per-query representative sequences (from step 5)
│   └── query1_id_representatives.fasta
├── query2_id/
│   └── query2_id_representatives.fasta
├── Taxonomic_Group_1/             # One directory per taxonomic group
│   ├── Taxonomic_Group_1_combined.fasta        # All queries + unique references
│   ├── Taxonomic_Group_1_trimmed.fasta         # Trimmed to aligned regions
│   ├── Taxonomic_Group_1_aligned.fasta         # Multiple sequence alignment
│   ├── Taxonomic_Group_1_tree.treefile         # Phylogenetic tree
│   ├── Taxonomic_Group_1_tree_labeled.treefile # Tree with taxonomic labels
│   └── Taxonomic_Group_1_tree.* (other IQTree files)
├── Taxonomic_Group_2/
│   └── ...
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
