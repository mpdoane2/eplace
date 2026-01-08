# Installation

1. Create a conda environment for ePLACE:

```
mamba create -yn eplace 'python>=3.12'
```

2. Install the required dependencies

```
mamba install -y bioconda::blast bioconda::pytaxonkit bioconda::iqtree bioconda::mafft
```

> Note 1.
> At the time of writing there is an [issue](https://github.com/bioforensics/pytaxonkit/issues/50) with conda not installing the
> most current version of pytaxonkit if you are using python >=3.12. This code works with older versions of pytaxonkit.


> Note 2.
> You will need to download and set up the NCBI taxonomy databases for pytaxonkit; see the [taxonkit documentation](https://bioinf.shenwei.me/taxonkit/usage/#before-use) for detailed instructions of which NCBI taxonomy files to download.

3. Get and install eplace

```bash
pip install git+https://github.com/linsalrob/eplace.git
```

After installation, the `eplace` command will be available in your environment:

```bash
# Verify installation
eplace --help

# Show version
eplace --version
```

## Usage

Once installed, you can use the `eplace` command with three subcommands:

- `eplace download` - Download NCBI BLAST database
- `eplace blast` - Run individual BLAST workflow (one tree per query)
- `eplace grouped` - Run grouped BLAST workflow (one tree per taxonomic group)

For detailed help on each command:
```bash
eplace download --help
eplace blast --help
eplace grouped --help
```

See the [README.md](README.md) for complete documentation and examples.



