# Installation

1. Create a conda environment for ePLACE:

```
mamba create -yn eplace 'python>=3.12'
```

2. Install the required dependencies

```
mamba install -y bioconda::blast bioconda::pytaxonkit
```

> Note 1.
> At the time of writing there is an [issue](https://github.com/bioforensics/pytaxonkit/issues/50) with conda not installing the
> most current version of pytaxonkit if you are using python >=3.12. This code works with older versions of pytaxonkit.


> Note 2.
> You will need to download and set up the NCBI taxonomy databases for pytaxonkit; see the [pytaxonkit documentation](https://github.com/bioforensics/pytaxonkit#setting-up-the-ncbi-taxonomy-database) for detailed instructions.

3. Get and install eplace

```
pip install git+https://github.com/linsalrob/eplace.git
```



