#!/bin/bash
#SBATCH --job-name=eplace_groups
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

eval "$(conda shell.bash hook)"
conda activate eplace

if [[ $# != 2 ]]; then echo "$0  <fasta file> <output dir> " >&2; exit 1; fi

# Note: This script does not build the combined tree (which can take many hours with large datasets).
# To build the combined tree, add: --combined-tree-label-rank genus (or another taxonomic rank)
eplace grouped --output-classification classifications.tsv --rank genus --group-rank class --tree-label-rank class --num-threads 32 $1 $2
