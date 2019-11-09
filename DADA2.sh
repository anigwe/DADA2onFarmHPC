#!/bin/bash -l

#SBATCH --job-name=Phenol16S
#SBATCH -o Phenol16S.%j.%N.out
#SBATCH -e Phenol16S.err
#SBATCH -p bigmemht
#SBATCH --ntasks=16
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aigwe@ucdavis.edu
#SBATCH --export==ALL

module load R
Rscript /home/aigwe/phenol/code/DADA2.R
