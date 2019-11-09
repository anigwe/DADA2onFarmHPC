#set working directory to sernonserp_summer2016 dataset folder by defining the variable path
#path <- "/home/aigwe/strepdrought/input/VannetteSoil16S"
#taxa.tax <- "/home/aigwe/strepdrought/input/silva_train/silva_nr_v128_train_set.fa"
#taxa.species <- "/home/aigwe/strepdrought/input/silva_train/silva_species_assignment_v128.fa"

#load required libraries and determine package versions: dada2, phyloseq, vegan and ggplot2
sink('/home/aigwe/strepdrought/output/strepdrought_package-version.txt')
cat("Package versions of dada2, phyloseq, vegan and ggplot2", "\n")
library(dada2)
library(phyloseq)
library(vegan)
library(ggplot2)
dada2version <- packageVersion("dada2")
phyloseqversion<- packageVersion("phyloseq")
veganversion<- packageVersion("vegan")
ggplot2version<- packageVersion("ggplot2")
cat("The package version of dada2 is...")
print(dada2version)
cat("The package version of phyloseq is...")
print(phyloseqversion)
cat("The package version of vegan is ...")
print(veganversion)
cat("The package version of ggplot2 is...")
print(ggplot2version)
sink()

#Read in the files in forward and reverse in matched order
sink('/home/aigwe/strepdrought/output/strepdrought_list_files.txt')
print(list.files("/home/aigwe/strepdrought/input/VannetteSoil16S"))
fnFs <- sort(list.files("/home/aigwe/strepdrought/input/VannetteSoil16S", pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files("/home/aigwe/strepdrought/input/VannetteSoil16S", pattern="_R2_001.fastq", full.names = TRUE))
sink()

#extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#view quality scores of the forward reads for the first two samples
pdf('/home/aigwe/strepdrought/output/strepdrought_plotQP_F.pdf')
plotQualityProfile(fnFs[1:2])
dev.off()
#view quality scores of the forward reads for the first two samples
pdf('/home/aigwe/strepdrought/output/strepdrought_plotQP_R.pdf')
plotQualityProfile(fnRs[1:2])
dev.off()

#creates new variable for filtered (trimmed) reads, and assigns file names for the new filtered reads
filt_path <- file.path("/home/aigwe/strepdrought/input/VannetteSoil16S", 'filtered')
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#trimming first 10 base pairs (trimLeft = 10) because the first tend to be poor quality with Illumina, and cutting off forward reads at 280 and reverse reads at 200 based on QC graphs 
#trim all reads, set at 280 for the forward and 200 for the reverse
#note that the F/R reads must still overlap after trimming for merging to work
#this also creates "filtered' as a directory as part of the process
#remember head() just views the first several files from the directory
sink('/home/aigwe/strepdrought/output/strepdrought_filt_trim.txt')
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
	trimLeft=10, truncLen=c(280,200), 
	maxN=0, maxEE=c(2,2), truncQ=2, 
	rm.phix=TRUE, compress=TRUE, multithread=TRUE) 
#head(out)
#View(out)
print(out)
sink()

#Randomly choosing 8 samples to look at to learn errors ("1:8")
errF <- learnErrors(filtFs[1:60], multithread=TRUE)
errR <- learnErrors(filtRs[1:60], multithread=TRUE)

#plot the error rates
pdf('/home/aigwe/strepdrought/output/strepdrought_plotErrF.pdf')
plotErrors(errF, nominalQ = TRUE)
dev.off()
pdf('/home/aigwe/strepdrought/output/strepdrought_plotErrR.pdf')
plotErrors(errR, nominalQ = TRUE)
dev.off()

#DEREPLICATION - this step matches exact sequences from all samples to each other and records their abundance
#this step finds and determines the exact sequence variants
sink('/home/aigwe/strepdrought/output/strepdrought_dereps.txt')
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
print(derepFs)
print(derepRs)
sink()

#see the number of unique sequences from each sample - also reports total number of seqs per sample
sink('/home/aigwe/strepdrought/output/strepdrought_unique_seqs.txt')
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
sink()
sink('/home/aigwe/strepdrought/output/strepdrought_unique_seqs2.txt')
cat("The number of unique sequences in the Forward samples are...")
print(dadaFs)
cat("The number of unique sequences in the Reverse samples are...")
print(dadaRs)
sink()

#describe data from the forward reads of the first sample
#sink('/home/aigwe/serpnonserp_summer2016/output/unique_seqs3.txt')
#cat("Describe data from the forward reads of the first sample")
#dadaFs[[1]]
#sink()

#merge the trimmed paired reads.
#all reads should merge, assumming trimming did not remove overlaps
#non exact matching sequences will be discared
sink('/home/aigwe/strepdrought/output/strepdrought_mergers.txt')
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
mergers_head <- head(mergers[1])
print(mergers)
cat("This is the head of merger 1...")
print(mergers_head)
sink()

#make sequence ('OTU') (actually ESV or ASV) table
sink('/home/aigwe/strepdrought/output/strepdrought_make_seq_ESV.txt')
seqtab <- makeSequenceTable(mergers)
print(dim(seqtab))
seqtabseq <- (table(nchar(getSequences(seqtab))))
print(seqtabseq)
print(dim(seqtabseq))
sink()

sink('/home/aigwe/strepdrought/output/strepdrought_seqtab.txt')
print(seqtab)
sink()

#remove chimeras
sink('/home/aigwe/strepdrought/output/strepdrought_nochim.txt')
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#determine the percentage of total sequence reads that are NOT chimeras (i.e. .96 = 4% are chimeras of total sequences)
#the percentage of reads that are chimeras should not be too high. No ballpark but >30% seems crazy
print(dim(seqtab.nochim))
print(sum(seqtab.nochim)/sum(seqtab))
sink()

sink('/home/aigwe/strepdrought/output/strepdrought_seqtab.nochim.txt')
print(seqtab.nochim)
sink()


#review number reads that have made it through the pipeline - trimmed, filtered, chimeria removal, etc.
sink('/home/aigwe/strepdrought/output/strepdrought_pipeline.txt')
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
print(colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim"))
print(rownames(track) <- sample.names)
#head(track)
#View(track)
print(track)
sink()

#ASSIGN TAXONOMY - dada gives the option of using RDP, greengeens, or Silva as the database for assigning taxonomy
#IMPORTANT - go to this site https://benjjneb.github.io/dada2/training.html and download the training set
#and place it in the same working directory 'path' before starting the command. For now just use the silva set.
#NOTE: Provide the full/absolute path to the Tax file otherwise you will get an error.
taxa <- assignTaxonomy(seqtab.nochim, "/home/aigwe/strepdrought/input/silva_train/silva_nr_v128_train_set.fa", multithread=TRUE)

#SPECIES ASSIGNMENT - optional additional method to make species assignments - available for silva and RDP
#also be sure to download the appropriate files and add to the same directory as 'path'
#NOTE: no need to unzip either of the above two tax files before using.
taxa <- addSpecies(taxa, "/home/aigwe/strepdrought/input/silva_train/silva_species_assignment_v128.fa")

#view the assigned taxonomic table
#NOTE:This is the end of the DADA2 pipeline - the tax table can now be exported, merged with metadata and opened
#with another package (Phyloseq, vegan, etc) for analysis
sink('/home/aigwe/strepdrought/output/strepdrought_taxa_print.txt')
taxa.print <- taxa
rownames(taxa.print) <- NULL
#head(taxa.print)
#View(taxa.print)
print(taxa.print)
sink()

#Export OTU (ESV) and tax tables to path for upload and analysis
setwd("/home/aigwe/strepdrought/output/")
table.tax <- write.table(taxa, "strepdrought_silva.txt", sep="\t", quote=F, row.names=T, col.names=NA)

table.nochim <- write.table(seqtab.nochim, "strepdrought_otu_dada.txt", sep="\t", quote=F, row.names=T, col.names=NA)

table.track <- write.table(track, "strepdrought_Reads.txt", sep="\t", quote=F, row.names=T, col.names=NA)

saveRDS(taxa, "strepdrought_silva.rds")
saveRDS(seqtab.nochim, "strepdrought_otu_dada.rds")
saveRDS(track, "strepdrought_Reads.rds")

sink('/home/aigwe/strepdrought/output/strepdrought_sessioninfo.txt')
print(sessionInfo())
sink()