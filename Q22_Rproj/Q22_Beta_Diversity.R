library(vegan) #v 2.6.2
library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(plyr)

setwd("C://Users/Jacobs Laboratory/Documents/JCYang/q22") # CHANGE to the directory containing the fastq files
here::i_am("Q22_Rproj/Q22_RS_Jensen_Shannon.R")
metadata <- read.delim("Q22_Microbiome/starting_files/Q22_Metadata.tsv",header=TRUE)
metadata$SampleID <- gsub("-",".",metadata$SampleID)
metadata$SampleID <- paste0("X","", metadata$SampleID)
metadata$Genotype <- revalue(metadata$Q22, replace = c("KO" = "Q22","WT"="WT"))
counts <- read.delim("Q22_Microbiome/starting_files/Q22_ASV.tsv",header=TRUE, row.names=1)

## Apply minimum sequencing depth threshold --
annotation$feature <- row.names(counts)
annotation$taxonomy <- counts$taxonomy
counts <- counts %>% select(-c(taxonomy))
summary(colSums(counts))
counts <- counts[colSums(counts) >= 3000]

## Split counts into colon and ileum -- 

# Ileum
ileum_meta <- metadata %>% filter(Site=="ILE", SampleID %in% names(counts))
row.names(ileum_meta) <- ileum_meta$SampleID
ileum <- ileum_meta$SampleID
ileum_counts <- counts %>% select(all_of(ileum))

# Colon
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

## Calculate Bray Curtis distance matrix
ileum_dis <- vegdist(t(ileum_counts), method = "bray")
meta <- metadata %>% filter(SampleID %in% names(ileum_counts))
meta <- meta[match(names(ileum_counts), meta$SampleID), ]
adonis2(ileum_dis ~ Litter+Sex+Q22, data = meta)

colon_dis <- vegdist(t(colon_counts), method = "bray")
meta <- metadata %>% filter(SampleID %in% names(colon_counts))
meta <- meta[match(names(colon_counts), meta$SampleID), ]
adonis2(colon_dis ~ Litter+Sex+Q22, data = meta)

## Calculate Robust Aitchison distance matrix
ileum_dis <- vegdist(t(ileum_counts), method = "robust.aitchison")
meta <- metadata %>% filter(SampleID %in% names(ileum_counts))
meta <- meta[match(names(ileum_counts), meta$SampleID), ]
set.seed(11)
adonis2(ileum_dis ~ Litter+Sex+Q22, data = meta)

colon_dis <- vegdist(t(colon_counts), method = "robust.aitchison")
meta <- metadata %>% filter(SampleID %in% names(colon_counts))
meta <- meta[match(names(colon_counts), meta$SampleID), ]
set.seed(11)
adonis2(colon_dis ~ Litter + Sex + Q22, data = meta)

## Calculate Jaccard distance matrix
ileum_dis <- vegdist(t(ileum_counts), method = "jaccard")
meta <- metadata %>% filter(SampleID %in% names(ileum_counts))
meta <- meta[match(names(ileum_counts), meta$SampleID), ]
adonis2(ileum_dis ~ Litter+Sex+Q22, data = meta)

colon_dis <- vegdist(t(colon_counts), method = "jaccard")
meta <- metadata %>% filter(SampleID %in% names(colon_counts))
meta <- meta[match(names(colon_counts), meta$SampleID), ]
adonis2(colon_dis ~ Litter+Sex+Q22, data = meta)



## Principal Coordinates Analysis -- 

cols <- c("WT"="black", "Q22"="red")


ileum_pcoa <- generate_pcoA_plots(distance_matrix=ileum.dist,
                                     counts = ileum_counts,
                                     metadata = ileum_meta,
                                     title="Ileum",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Q22_Microbiome/beta_diversity/ileum_Top_Taxa_RSJ_PcoA.csv")
ileum_pcoa



ileum_pcoa_roba <- generate_pcoA_plots(distance_matrix=ileum_dis,
                                  counts = ileum_counts,
                                  metadata = ileum_meta,
                                  title="Ileum",
                                  colorvariable = Genotype,
                                  colorvector = cols,
                                  wa_scores_filepath = "Q22_Microbiome/beta_diversity/ileum_Top_Taxa_Robust_Aitchison_PcoA.csv")
ileum_pcoa_roba

colon_pcoa <- generate_pcoA_plots(distance_matrix=colon.dist,
                                     counts = colon_counts,
                                     metadata = colon_meta,
                                     title="Colon",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Q22_Microbiome/beta_diversity/colon_Top_Taxa_PcoA.csv")
colon_pcoa 

colon_pcoa_roba <- generate_pcoA_plots(distance_matrix=colon_dis,
                                  counts = colon_counts,
                                  metadata = colon_meta,
                                  title="Colon",
                                  colorvariable = Genotype,
                                  colorvector = cols,
                                  wa_scores_filepath = "Q22_Microbiome/beta_diversity/colon_Top_Taxa_Robust_AitchisonPcoA.csv")
colon_pcoa_roba

## View taxa with greatest influence 
top_taxa <- read.csv("Q22_Microbiome/beta_diversity/ileum_Top_Taxa_Robust_Aitchison_PcoA.csv",header=TRUE)
top_taxa$feature <- top_taxa$sample
top_taxa <- merge(top_taxa, annotation, by="feature")
print(top_taxa$taxonomy)