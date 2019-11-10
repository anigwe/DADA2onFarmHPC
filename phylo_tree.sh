#!/bin/bash -l

#SBATCH --job-name=ctri-tree
#SBATCH -o ctri-tree.%j.%N.out
#SBATCH -e ctri-tree.err
#SBATCH -p bigmemht
#SBATCH --ntasks=16
#SBATCH --time=288:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aigwe@ucdavis.edu

#SBATCH --export==ALL

module load R
Rscript /home/aigwe/CTRI/code/ctri_phylo_tree.R
