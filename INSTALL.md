# Installation

1. Create a conda environment for ePLACE:

```
mamba create -yn eplace 'python>3.12'
```

2. Install the required dependencies

```
mamba install -y bioconda::blast bioconda::pytaxonkit
```

> Note
> At the time of writing there is an [issue](https://github.com/bioforensics/pytaxonkit/issues/50) with conda not installing the
> most current version of pytaxonkit if you are using python >=3.12. This code works with older versions of pytaxonkit.

3. Get eplace and install

```
pip install git+https://github.com/linsalrob/eplace.git
```



