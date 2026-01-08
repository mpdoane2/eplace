# ePLACE Command-Line Interface

The `eplace` command provides a unified interface to all ePLACE functionality through three subcommands.

## Installation

After installing the package with `pip install .` or `pip install -e .`, the `eplace` command will be available:

```bash
# Verify installation
eplace --help

# Check version
eplace --version
```

## Commands Overview

| Command | Description |
|---------|-------------|
| `eplace download` | Download NCBI BLAST database |
| `eplace blast` | Run individual BLAST workflow (one tree per query) |
| `eplace grouped` | Run grouped BLAST workflow (one tree per taxonomic group) |

## eplace download

Download and setup the NCBI core_nt BLAST database.

### Usage

```bash
eplace download [--force]
```

### Options

- `--force`: Force redownload even if database exists

### Examples

```bash
# Download database to default location ($BLASTDB or ~/blastdb)
eplace download

# Force redownload
eplace download --force
```

### Notes

- The download is large (several GB) and may take time
- Database will be stored in `$BLASTDB` if set, otherwise `~/blastdb`
- MD5 checksums are verified automatically

## eplace blast

Run BLAST search with individual taxonomy analysis. Creates one phylogenetic tree per query sequence.

### Usage

```bash
eplace blast QUERY_FASTA OUTPUT_DIR [OPTIONS]
```

### Required Arguments

- `QUERY_FASTA`: Path to query FASTA file containing sequences to search
- `OUTPUT_DIR`: Output directory for results (will be created if it doesn't exist)

### Optional Arguments

#### Taxonomy Options
- `--rank {phylum,class,order,family,genus,species}`: Taxonomic rank for representative selection (default: genus)
- `--tree-label-rank {phylum,class,order,family,genus,species}`: Taxonomic rank for tree labeling (default: genus)

#### Filtering Options
- `--min-identity FLOAT`: Minimum percent identity for BLAST hits (default: 90.0)
- `--min-coverage FLOAT`: Minimum query coverage percentage (default: 80.0)

#### Database Options
- `--database NAME`: BLAST database name (default: core_nt)
- `--blastdb-path PATH`: Path to BLAST database directory

#### Performance Options
- `--num-threads INT`: Number of threads for BLAST and alignment (default: 1)

#### Workflow Options
- `--overwrite-existing-blast`: Overwrite existing BLAST results
- `--skip-alignment`: Skip alignment and tree building steps
- `--output-classification PATH`: Path to output classification TSV file

### Examples

```bash
# Basic usage with default parameters
eplace blast query.fasta output_dir

# With custom parameters
eplace blast query.fasta output_dir \
    --rank genus \
    --min-identity 95 \
    --min-coverage 85 \
    --num-threads 4

# Skip alignment and tree building (BLAST only)
eplace blast query.fasta output_dir --skip-alignment

# Use custom BLAST database location
eplace blast query.fasta output_dir --blastdb-path /path/to/blastdb
```

### Output Structure

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
└── ...
```

## eplace grouped

Run BLAST search with grouped taxonomy analysis. Groups queries by taxonomic rank and creates one phylogenetic tree per group.

### Usage

```bash
eplace grouped QUERY_FASTA OUTPUT_DIR [OPTIONS]
```

### Required Arguments

- `QUERY_FASTA`: Path to query FASTA file containing sequences to search
- `OUTPUT_DIR`: Output directory for results (will be created if it doesn't exist)

### Optional Arguments

#### Taxonomy Options
- `--rank {phylum,class,order,family,genus,species}`: Taxonomic rank for representative selection (default: genus)
- `--group-rank {phylum,class,order,family,genus,species}`: Taxonomic rank for grouping sequences (default: class)
- `--tree-label-rank {phylum,class,order,family,genus,species}`: Taxonomic rank for tree labeling (default: genus)

#### Filtering Options
- `--min-identity FLOAT`: Minimum percent identity for BLAST hits (default: 90.0)
- `--min-coverage FLOAT`: Minimum query coverage percentage (default: 80.0)

#### Database Options
- `--database NAME`: BLAST database name (default: core_nt)
- `--blastdb-path PATH`: Path to BLAST database directory

#### Performance Options
- `--num-threads INT`: Number of threads for BLAST and alignment (default: 1)

#### Workflow Options
- `--overwrite-existing-blast`: Overwrite existing BLAST results
- `--skip-alignment`: Skip alignment and tree building steps
- `--alignment-tolerance INT`: Maximum coordinate difference for alignment consistency (default: 50)
- `--output-classification PATH`: Path to output classification TSV file

### Examples

```bash
# Basic usage (group by class, default)
eplace grouped query.fasta output_dir

# Group by different taxonomic rank
eplace grouped query.fasta output_dir --group-rank order

# Specify both representative and grouping ranks
eplace grouped query.fasta output_dir --rank genus --group-rank family

# Skip alignment and tree building
eplace grouped query.fasta output_dir --skip-alignment
```

### Output Structure

```
output_dir/
├── blast_results.txt              # Raw BLAST results
├── blast_results_annotated.txt    # BLAST results with taxonomic annotations
├── query1_id/                     # Per-query representative sequences
│   └── query1_id_representatives.fasta
├── Taxonomic_Group_1/             # One directory per taxonomic group
│   ├── Taxonomic_Group_1_combined.fasta        # All queries + unique references
│   ├── Taxonomic_Group_1_trimmed.fasta         # Trimmed to aligned regions
│   ├── Taxonomic_Group_1_aligned.fasta         # Multiple sequence alignment
│   ├── Taxonomic_Group_1_tree.treefile         # Phylogenetic tree
│   ├── Taxonomic_Group_1_tree_labeled.treefile # Tree with taxonomic labels
│   └── Taxonomic_Group_1_tree.* (other IQTree files)
└── ...
```

## Workflow Comparison

### Individual Workflow (eplace blast)

**Best for:** Analyzing each query sequence in its own phylogenetic context

**Process:**
1. Run BLAST search for all queries
2. Extract representative sequences for each query at specified rank
3. Create one directory per query
4. Build one alignment and tree per query

**Output:** Separate phylogenetic trees for each query sequence

### Grouped Workflow (eplace grouped)

**Best for:** Analyzing multiple related queries together in a single phylogenetic context

**Process:**
1. Run BLAST search for all queries
2. Extract representative sequences for each query
3. Group queries by specified taxonomic rank (e.g., class, order)
4. Combine all queries in a group with unique reference sequences
5. Build one alignment and tree per group

**Output:** Phylogenetic trees with multiple queries grouped by taxonomy

## Common Use Cases

### Quick BLAST search without trees

```bash
eplace blast query.fasta results --skip-alignment
```

### High stringency search

```bash
eplace blast query.fasta results \
    --min-identity 95 \
    --min-coverage 90
```

### Multi-threaded analysis

```bash
eplace blast query.fasta results --num-threads 8
```

### Group related sequences

```bash
eplace grouped query.fasta results \
    --rank species \
    --group-rank genus
```

## Troubleshooting

### Command not found

If `eplace` command is not found after installation:

```bash
# Check if it's installed
pip show eplace

# Reinstall
pip install --force-reinstall .

# Or add to PATH
export PATH="$HOME/.local/bin:$PATH"
```

### Dependencies missing

Some features require external tools:

- **BLAST+**: Required for all workflows
  ```bash
  sudo apt-get install ncbi-blast+  # Ubuntu/Debian
  brew install blast                 # macOS
  ```

- **TaxonKit**: Required for taxonomy lookups
  ```bash
  conda install -c bioconda taxonkit
  ```

- **MAFFT**: Required for alignment (unless --skip-alignment)
  ```bash
  sudo apt-get install mafft  # Ubuntu/Debian
  brew install mafft          # macOS
  ```

- **IQTree**: Required for tree building (unless --skip-alignment)
  ```bash
  sudo apt-get install iqtree  # Ubuntu/Debian
  brew install iqtree          # macOS
  ```

## See Also

- [README.md](../README.md) - Project overview and library API
- [INSTALL.md](../INSTALL.md) - Installation instructions
- [NCBI Database Download](ncbi_download.md) - Database management details
- [BLAST Workflow](blast_workflow.md) - Complete workflow guide
