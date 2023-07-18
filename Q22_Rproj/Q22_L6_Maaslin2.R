### Trios - Dataset wrangling ---
### Date : 4/4/23
library(here) #v 1.0.1
library(dplyr) #v 1.0.7
library(Maaslin2) #v 1.2.0
library(funrar)

setwd("/home/julianne/Documents/q22/") # CHANGE to the directory containing the fastq files
here::i_am("Q22_Rproj/Q22_L6_Maaslin2.R")

metadata <- read.delim("Q22_Microbiome/starting_files/Q22_Metadata.tsv", header=TRUE)
metadata$SampleID <- gsub("-",".",metadata$SampleID)
metadata$SampleID <- paste0("X","", metadata$SampleID)
counts <- read.delim("Q22_Microbiome/collapsed_ASV/export_L6_Q22_ASV/feature-table.tsv", header = TRUE, row.names=1)

## Apply minimum sequencing depth threshold --
counts <- counts %>% select(-c("taxonomy"))
counts <- counts[colSums(counts) >= 3000]

## Split counts into colon subsets -- 

# Luminal Colon 
lumcol_meta <- metadata %>% filter(Site=="COL", SampleID %in% names(counts))
row.names(lumcol_meta) <- lumcol_meta$SampleID
lumcol <- lumcol_meta$SampleID
lumcol_counts <- counts %>% select(all_of(lumcol))

# Ileum
ileum_meta <- metadata %>% filter(Site=="ILE", SampleID %in% names(counts))
row.names(ileum_meta) <- ileum_meta$SampleID
ileum <- ileum_meta$SampleID
ileum_counts <- counts %>% select(all_of(ileum))


## Read in L6 aggregated counts -- 
lc_genus <- lumcol_counts

ileum_genus <- ileum_counts

run_Maaslin2 <- function(counts_filepath, metadata_filepath, subset_string, fixed_effects_vector, results_string, random_effects_vector) {
  #input_data <- read.csv(counts_filepath, header=TRUE, row.names=1) # choose filtered non rarefied csv file
  
  df_input_data <- as.data.frame(counts_filepath)
  #df_input_data <- select(df_input_data, -c("taxonomy"))
  #df_input_data <- lumcol_counts
  
  transposed_input_data <- t(df_input_data)
  transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
  df_relative_ASV <- make_relative(transposed_input_data)
  df_relative_ASV <- as.data.frame(df_relative_ASV)
  Relative_Abundance <- summarize_all(df_relative_ASV, mean)
  Relative_Abundance <- as.data.frame(t(Relative_Abundance))
  
  readr::write_rds(Relative_Abundance,paste0("Q22_Microbiome/differential_taxa/Relative_Abundance-",subset_string,"-L6.RDS"))
  
  #input_metadata <-read.delim(metadata_filepath,sep="\t",header=TRUE, row.names=1)
  input_metadata <- as.data.frame(metadata_filepath)
  #input_metadata <- lumcol_meta
  
  target <- colnames(df_input_data)
  input_metadata = input_metadata[match(target, row.names(input_metadata)),]
  target == row.names(input_metadata)
  
  df_input_metadata<-input_metadata
  df_input_metadata$Sex <- factor(df_input_metadata$Sex, levels = c("Male", "Female"))
  df_input_metadata$Litter <- factor(df_input_metadata$Litter)
  df_input_metadata$Q22 <- factor(df_input_metadata$Q22, levels=c("WT","KO"))
  sapply(df_input_metadata,levels)
  
  ?Maaslin2
  fit_data = Maaslin2(input_data=df_input_data, 
                      input_metadata=df_input_metadata, 
                      output = paste0("Q22_Microbiome/differential_taxa/",subset_string, {{results_string}}), 
                      fixed_effects = {{fixed_effects_vector}},normalization="TSS",
                      random_effects= {{random_effects_vector}},
                      min_prevalence = 0.15,
                      transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)
}
# Luminal Colon
fixed_effects <- c("Sex", "Q22")
random_effects <- c("Litter")
run_Maaslin2(lc_genus,lumcol_meta,"Colon", fixed_effects, "_L6_Maaslin2_Sex_Q22_1-Litter", random_effects)

# Ileum
fixed_effects <- c("Sex", "Q22")
random_effects <- c("Litter")
run_Maaslin2(ileum_genus,ileum_meta,"Ileum", fixed_effects, "_L6_Maaslin2_Sex_Q22_1-Litter", random_effects)

# Luminal Colon: 0 significant features by Genotype 

data<-read.table("Q22_Microbiome/differential_taxa/Colon_L6_Maaslin2_Sex_Q22_1-Litter/all_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25) %>% filter(metadata=="Q22")

# Ileum

data<-read.table("Q22_Microbiome/differential_taxa/Ileum_L6_Maaslin2_Sex_Q22_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.05) %>% filter(metadata=="Q22")

data<- (merge(data, annotation, by = 'feature'))
print(data$taxonomy)
write.csv(data, "Q22_Microbiome/differential_taxa/Ileum_ASV_Maaslin2_Sex_Q22_1-Litter/annotated_significant_results.csv")
