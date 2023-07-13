library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
setwd("D:/2023_16S_data/July_2023_1/Julianne")

#save.image(file='D:/2023_16S_data/July_2023_1/DADA2_workspace.R')
load('D:/2023_16S_data/July_2023_1/Julianne/DADA2_workspace.R')

# File parsing
pathF <- "D:/2023_16S_data/July_2023_1/Julianne"
pathR <- "D:/2023_16S_data/July_2023_1/Julianne"
filtpathF <- file.path(pathF, "filtered")
filtpathR <- file.path(pathR, "filtered")
fastqFs <- sort(list.files(pathF, pattern="_L001_R1_001.fastq.gz"))  # 
fastqRs <- sort(list.files(pathR, pattern="_L001_R2_001.fastq.gz"))   #
plotQualityProfile(fastqFs[1:2])
plotQualityProfile(fastqRs[1:2])
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
out<-filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(240,200), maxEE=c(2,2), truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
head(out)

# File parsing
filtpathF <- "D:/2023_16S_data/July_2023_1/Julianne/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "D:/2023_16S_data/July_2023_1/Julianne/filtered" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="_L001_R1_001.fastq.gz", full.names = TRUE)   # _L001_
filtRs <- list.files(filtpathR, pattern="_L001_R2_001.fastq.gz", full.names = TRUE)   # _L001_
sample.names <- sapply(strsplit(basename(filtFs), "_S"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_S"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE,randomize=TRUE)
save(errF, file="errF.Rdata")
plotErrors(errF, nominalQ = TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE,randomize=TRUE)
save(errR, file="errR.Rdata")
plotErrors(errR, nominalQ = TRUE)
# Dereplication and DADA2 inference
derepFs <- derepFastq(filtFs, verbose = TRUE)
names(derepFs) <- sample.names
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)   # this step
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepRs) <- sample.names
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)

#saving RDS files 
saveRDS(seqtab, "seqtab_IL10.rds")
#After saving the seqtab file, Redo for each run and change the seqtab1.rds file to seqtab2.rds, seqtab3.rds etc. 

#Merging multiple runs 
#seqtab_2020_02_14 <- readRDS(choose.files()) 
#seqtab_2020_10_15 <- readRDS(choose.files()) 
#seqtab_2020_12_21 <- readRDS(choose.files()) 

#st.all <- mergeSequenceTables(seqtab_2020_02_14, seqtab_2020_10_15, seqtab_2020_12_21, seqtab)

# Remove chimeras
#seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)


# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Jonathan/Documents/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "C:/Users/Jonathan/Documents/silva_species_assignment_v138.1.fa.gz", verbose=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#Export to be compatible with QIIME2
write.table(t(seqtab.nochim), "seqtab-nochim.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab.nochim, fout='rep-seqs.fna', ids=colnames(seqtab.nochim))

# Export data so that it can be converted into BIOM
taxa[is.na(taxa)] <- ""
taxonomy<-paste("k__",taxa[,1],"; ","p__",taxa[,2],"; ","c__",taxa[,3],"; ","o__",taxa[,4],"; ","f__",taxa[,5],"; ","g__",taxa[,6],"; ","s__",taxa[,7],sep="")
output<-cbind(t(seqtab.nochim), taxonomy)
write.table(output, "D:/2023_16S_data/July_2023_1/Julianne/Julianne_July_2023_ASV_count_table.txt", sep="\t", col.names=NA)
##### Need to modify .txt file by typing "#OTU" in the upper left box, can then import into QIIME


#creating a phylogenetic tree: process takes a long time 
library("DECIPHER")
library("phangorn")
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "NNI", control = pml.control(trace = 1))
# pml.control(trace = 1) may be useful to check if optim.pml is still doing something or is crashed
#The most time consuming part will be the stochastic rearrangements
fitGTR <- optim.pml(fitGTR, rearrangement = "stochastic", rfadfatchet.par = list(iter = 5L, maxit = 5L, prop = 1/3))
#This is only 5 times instead of (a maximum) of 100 rearrangements. Now you could extrapolate how long phangorn would run and try
fitGTR <- optim.pml(fitGTR, rearrangement = "stochastic")
# this could be parallelized at the cost of memory consumption and in the end optimize all the other parameters again (they should not change much)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "NNI", control = pml.control(trace = 1))
detach("package:phangorn", unload=TRUE)

#Combined data to phyloseq
library(phyloseq)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa),phy_tree(fitGTR$tree))
ps


#Exporting tree
library(ape)
tree1 = phy_tree(ps)
ape::write.tree(tree1, "D:/M3LD_biopsy_Final2020/M3LD_finalSept2020_tree.tree")
