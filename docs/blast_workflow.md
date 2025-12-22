# BLAST Sequence Comparison Module

This module provides functionality for running BLAST searches, filtering results, and extracting representative sequences by taxonomic rank.

## Features

- Read FASTA files with validation
- Run blastn searches against NCBI databases
- Parse and filter BLAST results by identity and coverage
- Select representative sequences per taxonomic rank
- Extract sequences from BLAST databases
- Save results to separate FASTA files per query

## Quick Start

### Complete Workflow

```python
from pathlib import Path
from eplace_lib.blast_analysis import run_blast_search
from eplace_lib.taxonomy import process_blast_results_for_taxonomy

# Step 1: Run BLAST search with filtering
success, filtered_hits = run_blast_search(
    query_fasta=Path("query.fasta"),
    output_file=Path("blast_results.txt"),
    min_identity=90.0,  # 90% identity threshold
    min_coverage=80.0,   # 80% query coverage threshold
    database="core_nt",
    num_threads=4
)

# Step 2: Extract representative sequences by taxonomic rank
results = process_blast_results_for_taxonomy(
    blast_hits=filtered_hits,
    output_dir=Path("output"),
    rank="species"  # Can be: phylum, class, order, family, genus, species
)

# Results contains mapping of query IDs to output FASTA files
for query_id, output_fasta in results.items():
    print(f"{query_id}: {output_fasta}")
```

### Read FASTA Files

```python
from eplace_lib.blast_analysis import FastaReader

# Read sequences from FASTA file
sequences = FastaReader.read_fasta(Path("input.fasta"))

# Get sequence lengths
lengths = FastaReader.get_sequence_lengths(Path("input.fasta"))

for seq_id, length in lengths.items():
    print(f"{seq_id}: {length} bp")
```

### Run BLAST Search Only

```python
from eplace_lib.blast_analysis import BlastRunner

runner = BlastRunner()

# Check if blastn is available
if runner.check_blastn_available():
    # Run BLAST
    success = runner.run_blastn(
        query_fasta=Path("query.fasta"),
        output_file=Path("blast_results.txt"),
        database="core_nt",
        num_threads=4,
        max_target_seqs=100,
        evalue=1e-5
    )
```

### Parse and Filter BLAST Results

```python
from eplace_lib.blast_analysis import BlastRunner

runner = BlastRunner()

# Parse BLAST tabular output
hits = runner.parse_blast_results(Path("blast_results.txt"))

# Filter hits by identity and coverage
filtered_hits = runner.filter_blast_hits(
    hits,
    min_identity=90.0,
    min_coverage=80.0
)

print(f"Filtered {len(hits)} hits to {len(filtered_hits)} hits")
```

### Select Representatives by Taxonomic Rank

```python
from eplace_lib.taxonomy import TaxonomyExtractor

extractor = TaxonomyExtractor()

# Group hits by query
grouped_hits = extractor.group_hits_by_query(blast_hits)

# Select representatives for each query
for query_id, query_hits in grouped_hits.items():
    representatives = extractor.select_representatives_by_rank(
        hits=query_hits,
        rank="species",  # Taxonomic rank
        max_per_rank=1   # Number per rank
    )
    print(f"{query_id}: {len(representatives)} representatives")
```

### Extract Sequences from Database

```python
from eplace_lib.taxonomy import SequenceExtractor

extractor = SequenceExtractor()

# Check if blastdbcmd is available
if extractor.check_blastdbcmd_available():
    # Extract sequences
    success = extractor.extract_sequences(
        sequence_ids=["NC_001234.5", "NC_005678.9"],
        output_fasta=Path("extracted.fasta"),
        database="core_nt"
    )
```

## Command-Line Usage

The package includes an example script for command-line usage:

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

### Command-Line Options

- `query_fasta`: Path to query FASTA file (required)
- `output_dir`: Directory for output files (required)
- `--rank`: Taxonomic rank (phylum, class, order, family, genus, species) [default: species]
- `--min-identity`: Minimum percent identity [default: 90.0]
- `--min-coverage`: Minimum query coverage percentage [default: 80.0]
- `--database`: BLAST database name [default: core_nt]
- `--blastdb-path`: Path to BLAST database directory [default: $BLASTDB or ~/blastdb]
- `--num-threads`: Number of threads for BLAST [default: 1]
- `--dry-run`: Display what would be done without running BLAST

## Output Structure

The workflow creates the following output structure:

```
output_dir/
├── blast_results.txt              # Raw BLAST results (tabular format)
├── query1_id/
│   └── query1_id_representatives.fasta
├── query2_id/
│   └── query2_id_representatives.fasta
└── ...
```

Each query sequence gets its own directory with a FASTA file containing representative sequences.

## Data Classes

### BlastHit

Represents a single BLAST hit with the following attributes:

- `query_id`: Query sequence identifier
- `subject_id`: Subject (database) sequence identifier
- `percent_identity`: Percentage of identical matches
- `alignment_length`: Length of alignment
- `query_length`: Length of query sequence
- `subject_length`: Length of subject sequence
- `query_start`: Start position in query
- `query_end`: End position in query
- `subject_start`: Start position in subject
- `subject_end`: End position in subject
- `evalue`: Expectation value
- `bit_score`: Bit score
- `query_coverage`: Percentage of query covered by alignment

### TaxonomicInfo

Represents taxonomic information for a sequence:

- `sequence_id`: Sequence identifier
- `taxid`: NCBI taxonomy ID
- `kingdom`: Kingdom name
- `phylum`: Phylum name
- `class_name`: Class name
- `order`: Order name
- `family`: Family name
- `genus`: Genus name
- `species`: Species name

## Requirements

### System Requirements

- BLAST+ tools must be installed:
  - `blastn`: For running BLAST searches
  - `blastdbcmd`: For extracting sequences from databases
- NCBI BLAST database (e.g., core_nt) must be downloaded

### Python Requirements

Uses Python standard library modules only:
- `os`, `subprocess`, `logging`, `pathlib`
- `typing`, `collections`, `dataclasses`

No external dependencies required.

## Testing

Run the test suite:

```bash
# Test BLAST analysis module
pytest tests/test_blast_analysis.py -v

# Test taxonomy module
pytest tests/test_taxonomy.py -v

# Test complete workflow
pytest tests/test_workflow.py -v

# Run all tests
pytest tests/ -v
```

## Examples

See `examples/blast_workflow_example.py` for a comprehensive example demonstrating the complete workflow.

## Performance Considerations

- BLAST searches can be slow for large query files or databases
- Use `--num-threads` to parallelize BLAST searches
- Consider splitting large query files into smaller batches
- Ensure sufficient disk space for BLAST output and extracted sequences

## Limitations

- Taxonomic information extraction is simplified in this version
- Full taxonomy resolution would require querying NCBI taxonomy database
- Representative selection uses sequence ID patterns as a proxy for taxonomic grouping
- For production use, consider integrating with NCBI Entrez or taxonomy databases

## Security

The module includes security features:

- Input validation for FASTA files
- Safe subprocess execution with timeouts
- Path validation to prevent directory traversal
- Sanitized filenames for output files

## Troubleshooting

### "blastn is not available"

Install BLAST+ tools:
```bash
# Ubuntu/Debian
sudo apt-get install ncbi-blast+

# macOS with Homebrew
brew install blast

# Or download from NCBI
# https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

### "blastdbcmd is not available"

blastdbcmd is included with BLAST+ tools. Install BLAST+ as shown above.

### "Database already exists"

If you want to use a different database location, set the `BLASTDB` environment variable:

```bash
export BLASTDB=/path/to/your/blastdb
```

### "No hits found"

Check your filtering parameters:
- Try lowering `--min-identity` threshold
- Try lowering `--min-coverage` threshold
- Ensure your query sequences are appropriate for the database

## Future Enhancements

Potential improvements for future versions:

- Integration with NCBI Entrez API for full taxonomy information
- Support for additional BLAST programs (blastp, blastx, etc.)
- Multiple sequence alignment of representative sequences
- Phylogenetic tree construction from representatives
- Web interface for the workflow
- Support for custom databases beyond core_nt
