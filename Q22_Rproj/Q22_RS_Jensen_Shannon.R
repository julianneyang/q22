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
muccol_meta <- metadata %>% filter(Subset=="Mucosal_Colon", SampleID %in% names(counts))
row.names(muccol_meta) <- muccol_meta$SampleID
muccol <- muccol_meta$SampleID
muccol_counts <- counts %>% select(all_of(muccol))

# Luminal Colon no HET 
nohet_lt_ileum_meta <- metadata %>%
  filter(Genotype!="HET")%>%
  filter(Subset=="Luminal_Colon", SampleID %in% names(counts))
row.names(nohet_lt_ileum_meta) <- nohet_lt_ileum_meta$SampleID
ileum <- nohet_lt_ileum_meta$SampleID
nohet_ileum_counts <- counts %>% select(all_of(ileum))

# Mucosal Colon no HET
nohet_lt_muccol_meta <- metadata %>% 
  filter(Genotype!="HET")%>%
  filter(Subset=="Mucosal_Colon", SampleID %in% names(counts))
row.names(nohet_lt_muccol_meta) <- nohet_lt_muccol_meta$SampleID
muccol <- nohet_lt_muccol_meta$SampleID
nohet_muccol_counts <- counts %>% select(all_of(muccol))

## Prevalence filter datasets -- 
# Ileum
0.15*28 #4 samples
ileum_counts <- prevalence_filter(ileum_counts,4)

# Mucosal Colon 
0.15*89
muccol_counts <- prevalence_filter(muccol_counts,13)

# Luminal Colon - no HET 
0.15*42 #6 samples
nohet_lt_lumcol_counts <- prevalence_filter(nohet_lumcol_counts,6)

# Mucosal Colon - no HET 
0.15*41 #6 samples
nohet_lt_muccol_counts <- prevalence_filter(nohet_muccol_counts,6)

## Calculate RS Jensen Shannon distance matrix -- 


muccol.dist <- calculate_rsjensen(muccol_counts)
ileum.dist <- calculate_rsjensen(ileum_counts)
lt_ileum.dist <- calculate_rsjensen(nohet_lt_ileum_counts)
lt_muccol.dist <- calculate_rsjensen(nohet_lt_muccol_counts)

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

slt_mc_pcoa <- generate_pcoA_plots(distance_matrix=muccol.dist,
                                     counts = muccol_counts,
                                     metadata = muccol_meta,
                                     title="Long Term - Mucosal Colon RS Jensen",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Long_Term/MucCol_Top_Taxa_PcoA.csv")
nohet_slt_lc_pcoa <- generate_pcoA_plots(distance_matrix=lt_ileum.dist,
                                   counts = nohet_ileum_counts,
                                   metadata = nohet_lt_ileum_meta,
                                   title="Long Term - Luminal Colon RS Jensen",
                                   colorvariable = Genotype,
                                   colorvector = cols,
                                   wa_scores_filepath = "Long_Term/nohet_ileum_Top_Taxa_PcoA.csv")

nohet_slt_mc_pcoa <- generate_pcoA_plots(distance_matrix=lt_muccol.dist,
                                   counts = nohet_muccol_counts,
                                   metadata = nohet_lt_muccol_meta,
                                   title="Long Term - Mucosal Colon RS Jensen",
                                   colorvariable = Genotype,
                                   colorvector = cols,
                                   wa_scores_filepath = "Long_Term/nohet_MucCol_Top_Taxa_PcoA.csv")

plot_grid(slt_lc_pcoa,slt_mc_pcoa, labels=c("C","D"),label_size = 20)
plot_grid(nohet_slt_lc_pcoa,nohet_slt_mc_pcoa, labels=c("C","D"),label_size = 20)


## PERMANOVA

# Luminal Colon 
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
data.dist<-muccol.dist
metadata <- muccol_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sex + Site + Genotype, data=metadata, permutations=10000)
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