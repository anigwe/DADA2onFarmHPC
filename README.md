# DADA2 on Farm HPC
Contains code used to run DADA2 and create a phylogenetic tree on the Farm (research and teaching cluster for the College of Agricultural and Environmental Sciences) and using the terminal.
- Download Cyberduck to easily interface with folders in the HPC.
  - SFTP (SSH file transfer protocol), server: agri.cse.ucdavis.edu, port: 2022

- DADA2.R contains the code used for the DADA2 pipeline
- DADA2.sh job script (include absolute file location of DADA2.R)
- sbatch ~/DADA2.sh in the terminal to begin process

Run R interactively until strepdrought_plotErrF.pdf and strepdrought_plotErrr.pdf is produced. Download those files and use them to determine where to trim forward and reverse reads. 

## DADA2 inputs
- Sequencing files formatted correctly for input into DADA2

## DADA2 outputs
### major - OTU (ESV) and tax tables to path for upload and analysis
- strepdrought_silva.rds
  - saveRDS(taxa, "strepdrought_silva.rds")
- strepdrought_otu_dada.rds
  - saveRDS(seqtab.nochim, "strepdrought_otu_dada.rds")
- strepdrought_Reads.rds
  - saveRDS(track, "strepdrought_Reads.rds")
### minor
- strepdrought_package-version.txt 
  - list of package versions
- strepdrought_list_files.txt
  -  list of forward and reverse reads
- strepdrought_plotQP_F.pdf
- strepdrought_plotQP_R.pdf
  - plot of quality scores of the forward/reverse reads for the first two samples
- strepdrought_filt_trim.txt
  - list of reads in and reads out for each
- strepdrought_plotErrF.pdf
- strepdrought_plotErrr.pdf
  - plot of forward/reverse error rates
- strepdrought_dereps.txt
  - text of R object describing dereplicated sequencing reads
- strepdrought_unique_seqs.txt
- strepdrought_unique_seqs2.txt
  - see the number of unique sequences from each sample - also reports total number of seqs per sample
- strepdrought_mergers.txt
- strepdrought_make_seq_ESV.txt
- strepdrought_seqtab.txt
- strepdrought_nochim.txt
- strepdrought_seqtab.nochim.txt
- strepdrought_pipeline.txt
- strepdrought_sessioninfo.txt
  - session info

## phangorn inputs (phylogenetic tree)
- seqtab.nochim (strepdrought_otu_dada.rds)

## phangorn outputs
- ctri_phylo_tree.R
  - Contains R code used for constructing the phylogenetic tree

- ctri-tree.23909765.%N.out
  - Standard output file with the job number in the name. Contains information related to job run time and memory usage.

- ctri-tree.23909765.bigmem10.out
  - Output file named with job number and the node the job landed on. Contains information related to the R environment.

- ctri-tree.err
  - Errors written to this file

- phylo_tree.sh
  - Script that was used to sbatch file containing R code

- seqtab.nochim_alignment.RDS
  - RDS of alignment step

- seqtab.nochim_fit.RDS
  - RDS of phangorn fit step

- seqtab.nochim_tree.RData
  - Saved disk image of R environment

- seqtab.nochim_tree.RDS
  - RDS of DECIPHER constructed tree
    - Tree can be accessed in phyloseq via the following commands:
      - fitGTR <- readRDS("seqtab.nochim_tree.RDS")
      - TREE = phy_tree(fitGTR$tree)
