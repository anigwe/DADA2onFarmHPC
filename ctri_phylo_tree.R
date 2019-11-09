library(dada2)
library(phangorn)
library(DECIPHER)

seqtab <- readRDS('/home/aigwe/CTRI/seqtab.nochim.rds')

seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

saveRDS(alignment, "seqtab.nochim_alignment.RDS")

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

saveRDS(fit, "seqtab.nochim_fit.RDS")

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

saveRDS(fitGTR, "seqtab.nochim_tree.RDS")
save.image(file = "seqtab.nochim_tree.RData")