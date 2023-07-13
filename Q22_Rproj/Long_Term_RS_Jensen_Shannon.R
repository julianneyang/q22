library(vegan) #v 2.6.2

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Long_Term_RS_Jensen_Shannon.R")

metadata <- read.table("Long_Term/starting_files/SLC_LT_metadata.tsv", header=TRUE)
counts <- read.table("Long_Term/starting_files/SLT_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Split counts into colon subsets -- 

# Luminal Colon 
lumcol_meta <- metadata %>% filter(Subset=="Luminal_Colon", SampleID %in% names(counts))
row.names(lumcol_meta) <- lumcol_meta$SampleID
lumcol <- lumcol_meta$SampleID
lumcol_counts <- counts %>% select(all_of(lumcol))

# Mucosal Colon
muccol_meta <- metadata %>% filter(Subset=="Mucosal_Colon", SampleID %in% names(counts))
row.names(muccol_meta) <- muccol_meta$SampleID
muccol <- muccol_meta$SampleID
muccol_counts <- counts %>% select(all_of(muccol))

# Luminal Colon no HET 
nohet_lt_lumcol_meta <- metadata %>%
  filter(Genotype!="HET")%>%
  filter(Subset=="Luminal_Colon", SampleID %in% names(counts))
row.names(nohet_lt_lumcol_meta) <- nohet_lt_lumcol_meta$SampleID
lumcol <- nohet_lt_lumcol_meta$SampleID
nohet_lumcol_counts <- counts %>% select(all_of(lumcol))

# Mucosal Colon no HET
nohet_lt_muccol_meta <- metadata %>% 
  filter(Genotype!="HET")%>%
  filter(Subset=="Mucosal_Colon", SampleID %in% names(counts))
row.names(nohet_lt_muccol_meta) <- nohet_lt_muccol_meta$SampleID
muccol <- nohet_lt_muccol_meta$SampleID
nohet_muccol_counts <- counts %>% select(all_of(muccol))

## Prevalence filter datasets -- 
# Luminal Colon
0.15*90 #13 samples
lumcol_counts <- prevalence_filter(lumcol_counts,13)

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
lumcol.dist <- calculate_rsjensen(lumcol_counts)
lt_lumcol.dist <- calculate_rsjensen(nohet_lt_lumcol_counts)
lt_muccol.dist <- calculate_rsjensen(nohet_lt_muccol_counts)

## Principal Coordinates Analysis -- 

cols <- c("WT"="black", "HET"= "blue", "MUT"="red")

slt_lc_pcoa <- generate_pcoA_plots(distance_matrix=lumcol.dist,
                                     counts = lumcol_counts,
                                     metadata = lumcol_meta,
                                     title="Long Term - Luminal Colon RS Jensen",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Long_Term/LumCol_Top_Taxa_PcoA.csv")

slt_mc_pcoa <- generate_pcoA_plots(distance_matrix=muccol.dist,
                                     counts = muccol_counts,
                                     metadata = muccol_meta,
                                     title="Long Term - Mucosal Colon RS Jensen",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Long_Term/MucCol_Top_Taxa_PcoA.csv")
nohet_slt_lc_pcoa <- generate_pcoA_plots(distance_matrix=lt_lumcol.dist,
                                   counts = nohet_lumcol_counts,
                                   metadata = nohet_lt_lumcol_meta,
                                   title="Long Term - Luminal Colon RS Jensen",
                                   colorvariable = Genotype,
                                   colorvector = cols,
                                   wa_scores_filepath = "Long_Term/nohet_LumCol_Top_Taxa_PcoA.csv")

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
data.dist<-lumcol.dist
metadata <- lumcol_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sex + Site + Genotype, data=metadata, permutations=10000)
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
data.dist<-lt_lumcol.dist
metadata <- nohet_lt_lumcol_meta

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
run_repeated_PERMANOVA(lumcol.dist,
                       lumcol_meta,
                       site,
                       mouseID,
                       order_vector)
site <- c("Site")
mouseID <- c("Sex","Genotype","MouseID")
order_vector <- c("Sex","Site","Genotype")
run_repeated_PERMANOVA(lt_lumcol.dist,
                       nohet_lt_lumcol_meta,
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