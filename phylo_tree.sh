#!/bin/bash -l

#SBATCH --job-name=Strepdrought-tree
#SBATCH -o Strepdrought-tree.%j.%N.out
#SBATCH -e Strepdrought-tree.err
#SBATCH -p bigmemht
#SBATCH --ntasks=16
#SBATCH --time=288:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aigwe@ucdavis.edu

#SBATCH --export==ALL

Rscript /home/aigwe/strepdrought/code/phylo_tree.R