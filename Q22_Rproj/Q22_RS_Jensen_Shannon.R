library(vegan) #v 2.6.2
library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/q22/") # CHANGE to the directory containing the fastq files
here::i_am("Q22_Rproj/Q22_RS_Jensen_Shannon.R")

metadata <- read.delim("Q22_Microbiome/starting_files/Q22_Metadata.tsv", header=TRUE)
metadata$SampleID <- gsub("-",".",metadata$SampleID)
metadata$SampleID <- paste0("X","", metadata$SampleID)
counts <- read.delim("Q22_Microbiome/starting_files/Q22_ASV.tsv", header = TRUE,row.names=1)

## Apply minimum sequencing depth threshold --
counts <- counts %>% select(-c(taxonomy))
summary(colSums(counts))
counts <- counts[colSums(counts) >= 3000]

## Split counts into colon and ileum -- 

# Ileum
ileum_meta <- metadata %>% filter(Site=="ILE", SampleID %in% names(counts))
row.names(ileum_meta) <- ileum_meta$SampleID
ileum <- ileum_meta$SampleID
ileum_counts <- counts %>% select(all_of(ileum))

# Mucosal Colon
colon_meta <- metadata %>% filter(Site=="COL", SampleID %in% names(counts))
row.names(colon_meta) <- colon_meta$SampleID
colon <- colon_meta$SampleID
colon_counts <- counts %>% select(all_of(colon))

## Prevalence filter datasets -- 
# Ileum
0.15*28 #4 samples
ileum_counts <- prevalence_filter(ileum_counts,4)

#  Colon 
0.15*28 #4 samples
colon_counts <- prevalence_filter(colon_counts,4)

## Calculate RS Jensen Shannon distance matrix -- 


colon.dist <- calculate_rsjensen(colon_counts)
ileum.dist <- calculate_rsjensen(ileum_counts)

## Principal Coordinates Analysis -- 

cols <- c("WT"="black", "KO"="red")


ileum_pcoa <- generate_pcoA_plots(distance_matrix=ileum.dist,
                                     counts = ileum_counts,
                                     metadata = ileum_meta,
                                     title="Long Term - Luminal Colon RS Jensen",
                                     colorvariable = Q22,
                                     colorvector = cols,
                                     wa_scores_filepath = "Q22_Microbiome/beta_diversity/ileum_Top_Taxa_RSJ_PcoA.csv")
ileum_pcoa

unique(ileum_meta$Litter)
cols <- viridis(13)
ileum_pcoa_Litter <- generate_pcoA_plots(distance_matrix=ileum.dist,
                                  counts = ileum_counts,
                                  metadata = ileum_meta,
                                  title="Long Term - Luminal Colon RS Jensen",
                                  colorvariable = Litter,
                                  colorvector = cols,
                                  wa_scores_filepath = "Q22_Microbiome/beta_diversity/ileum_Top_Taxa_RSJ_PcoA.csv")
ileum_pcoa_Litter

cols <- c("WT"="black", "KO"="red")
colon_pcoa <- generate_pcoA_plots(distance_matrix=colon.dist,
                                     counts = colon_counts,
                                     metadata = colon_meta,
                                     title="Long Term - Mucosal Colon RS Jensen",
                                     colorvariable = Q22,
                                     colorvector = cols,
                                     wa_scores_filepath = "Q22_Microbiome/beta_diversity/colon_Top_Taxa_PcoA.csv")
colon_pcoa 
plot_grid(slt_lc_pcoa,slt_mc_pcoa, labels=c("C","D"),label_size = 20)
plot_grid(nohet_slt_lc_pcoa,nohet_slt_mc_pcoa, labels=c("C","D"),label_size = 20)


## PERMANOVA

# Ileum
data.dist<-ileum.dist
metadata <- ileum_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Litter+ Sex + Q22, data=metadata, permutations=10000)
data.adonis$aov.tab

# Mucosal Colon
data.dist<-colon.dist
metadata <- colon_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sex + Q22, data=metadata, permutations=10000)
data.adonis$aov.tab

# Luminal Colon -- no HET 
data.dist<-lt_ileum.dist
metadata <- nohet_lt_ileum_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

# Mucosal Colon -- no HET
data.dist<-lt_muccol.dist
metadata <- nohet_lt_muccol_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab


## Repeat-Measures-Aware 
# Luminal Colon
site <- c("Site")
mouseID <- c("Sex","Genotype","MouseID")
order_vector <- c("Sex","Site","Genotype")
run_repeated_PERMANOVA(ileum.dist,
                       ileum_meta,
                       site,
                       mouseID,
                       order_vector)
site <- c("Site")
mouseID <- c("Sex","Genotype","MouseID")
order_vector <- c("Sex","Site","Genotype")
run_repeated_PERMANOVA(lt_ileum.dist,
                       nohet_lt_ileum_meta,
                       site,
                       mouseID,
                       order_vector)

# Mucosal Colon
site <- c("Site")
mouseID <- c("Sex","Genotype","MouseID")
order_vector <- c("Sex","Site","Genotype")
run_repeated_PERMANOVA(muccol.dist,
                       muccol_meta,
                       site,
                       mouseID,
                       order_vector)
site <- c("Site")
mouseID <- c("Sex","Genotype","MouseID")
order_vector <- c("Sex","Site","Genotype")
run_repeated_PERMANOVA(lt_muccol.dist,
                       nohet_lt_muccol_meta,
                       site,
                       mouseID,
                       order_vector)