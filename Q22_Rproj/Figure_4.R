library(ggplot2)
library(rlang)
library(cowplot)
library(viridis)
library(tidyr)
library(dplyr)

# Read in startle data 
data<-read.csv("Q22_Behavior/Startle_PPI_Analysis - Analysis.csv", header=TRUE)
data$SLC_Genotype<-data$SLC
data$SLC_Genotype <- factor(data$SLC_Genotype, levels=c("WT", "HET", "MUT"))
data$Q22 <- factor(data$Q22, levels=c("WT","KO"))
data_slcwt <- data %>% filter(SLC_Genotype=="WT")

# Grab significantly diff abundant taxa 
significant_feature <- read.csv("Q22_Microbiome/differential_taxa/Ileum_ASV_Maaslin2_Sex_Q22_1-Litter/annotated_significant_results.csv")
streptococcus <- significant_feature$feature

#Grab Ileum counts
metadata <- read.delim("Q22_Microbiome/starting_files/Q22_Metadata.tsv", header=TRUE)
metadata$SampleID <- gsub("-",".",metadata$SampleID)
metadata$SampleID <- paste0("X","", metadata$SampleID)
counts <- read.delim("Q22_Microbiome/starting_files/Q22_ASV.tsv", header = TRUE,row.names=1)

ileum_meta <- metadata %>% filter(Site=="ILE", SampleID %in% names(counts))
row.names(ileum_meta) <- ileum_meta$SampleID
ileum <- ileum_meta$SampleID
ileum_counts <- counts %>% select(all_of(ileum))

# For Relative Abundances 
transposed_input_data <- t(ileum_counts)
transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
df_relative_ASV <- make_relative(transposed_input_data)
df_relative_ASV <- as.data.frame(df_relative_ASV)
rowSums(df_relative_ASV) #Every sample sums to 1
ileum_counts <- df_relative_ASV %>% select(streptococcus)

#For Absolute Abundances
ileum_counts <- ileum_counts[streptococcus,]
ileum_counts <- as.data.frame(t(ileum_counts))

# Continue with merging ileum counts with metadata
ileum_counts$SampleID <- row.names(ileum_counts)
ileum_counts_meta <- merge(ileum_counts, metadata, by="SampleID")


# Merge startle data and streptococcus counts 
correlate_microbiome <- merge(data_slcwt, ileum_counts_meta, by="MouseID")
colnames(correlate_microbiome)[colnames(correlate_microbiome) == "TACGTAGGTCCCGAGCGTTATCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTGGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTATGCTTTGGAAACTGTTCAACTTGAGTGCAGAAGGGGAGAGTGGAATTCCATGTGTAGCGGTGGAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGTAGCGAACAGG"] <- "Streptococcus"

# Pearson
correlation <- cor(log(correlate_microbiome$Streptococcus+0.5), 
                   correlate_microbiome$First_average_VMax, method = "pearson")
p_value <- cor.test(log(correlate_microbiome$Streptococcus+0.5), 
                    correlate_microbiome$First_average_VMax, method = "pearson")$p.value
correlation <- cor(log(correlate_microbiome$Streptococcus+0.5), 
                   correlate_microbiome$middle_average_VMax, method = "pearson")
p_value <- cor.test(log(correlate_microbiome$Streptococcus+0.5), 
                    correlate_microbiome$middle_average_VMax, method = "pearson")$p.value
correlation <- cor(log(correlate_microbiome$Streptococcus+0.5), 
                   correlate_microbiome$last_average_VMax, method = "pearson")
p_value <- cor.test(log(correlate_microbiome$Streptococcus+0.5), 
                    correlate_microbiome$last_average_VMax, method = "pearson")$p.value

#Spearman - with no log transformation
correlation <- cor(correlate_microbiome$Streptococcus, 
                   correlate_microbiome$First_average_VMax, method = "spearman")
p_value <- cor.test(correlate_microbiome$Streptococcus, 
                    correlate_microbiome$First_average_VMax, method = "spearman")$p.value
truncated_correlation <- round(correlation, 4)
truncated_p_value <- round(p_value, 4)

spearman_vs_first_startle <- ggplot(correlate_microbiome, aes(Streptococcus, First_average_VMax)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = paste("First Startle"),
       subtitle = paste(paste("correlation=", truncated_correlation), paste("p=", truncated_p_value))) + 
  ylab("Average VMax") +
  xlab("Relative Abundance") +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle =  element_text(hjust = 0.5)) 

correlation <- cor(correlate_microbiome$Streptococcus, 
                   correlate_microbiome$middle_average_VMax, method = "spearman")
p_value <- cor.test(correlate_microbiome$Streptococcus, 
                    correlate_microbiome$middle_average_VMax, method = "spearman")$p.value
truncated_correlation <- round(correlation, 4)
truncated_p_value <- round(p_value, 4)

spearman_vs_middle_startle <- ggplot(correlate_microbiome, aes(Streptococcus, middle_average_VMax)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = paste("Middle Startle"),
       subtitle = paste(paste("correlation=", truncated_correlation), paste("p=", truncated_p_value))) + 
  ylab("Average VMax") +
  xlab("Relative Abundance") +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle =  element_text(hjust = 0.5)) 

correlation <- cor(correlate_microbiome$Streptococcus, 
                   correlate_microbiome$last_average_VMax, method = "spearman")
p_value <- cor.test(correlate_microbiome$Streptococcus, 
                    correlate_microbiome$last_average_VMax, method = "spearman")$p.value
truncated_correlation <- round(correlation, 4)
truncated_p_value <- round(p_value, 4)

spearman_vs_last_startle <- ggplot(correlate_microbiome, aes(Streptococcus,last_average_VMax)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = paste("Last Startle"),
       subtitle = paste(paste("correlation=", truncated_correlation), paste("p=", truncated_p_value))) + 
  ylab("Average VMax") +
  xlab("Relative Abundance") +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle =  element_text(hjust = 0.5)) 


dat <- ggplot(correlate_microbiome, aes(Q22.x,Streptococcus, fill=Q22.x))+ 
  geom_boxplot(alpha=0.25)+
  scale_fill_viridis_d()+
  geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
  theme_cowplot(12) +
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title="Streptococcus ASV")+
  ylab("Relative Abundance")+
  xlab("")

plot_grid(dat, spearman_vs_first_startle,
          spearman_vs_middle_startle, spearman_vs_last_startle,
          labels=c("A","B","C","D"))

#Spearman - with log transformation
correlation <- cor(correlate_microbiome$Streptococcus, 
                   correlate_microbiome$First_average_VMax, method = "spearman")
p_value <- cor.test(correlate_microbiome$Streptococcus, 
                    correlate_microbiome$First_average_VMax, method = "spearman")$p.value
truncated_correlation <- round(correlation, 4)
truncated_p_value <- round(p_value, 4)

spearman_vs_first_startle <- ggplot(correlate_microbiome, aes(log(Streptococcus+pseudocount), First_average_VMax)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = paste("Spearman Correlation:", truncated_correlation),
       subtitle = paste("p-value:", truncated_p_value)) + 
  theme_cowplot(12)

correlation <- cor(correlate_microbiome$Streptococcus, 
                   correlate_microbiome$middle_average_VMax, method = "spearman")
p_value <- cor.test(correlate_microbiome$Streptococcus, 
                    correlate_microbiome$middle_average_VMax, method = "spearman")$p.value
truncated_correlation <- round(correlation, 4)
truncated_p_value <- round(p_value, 4)

spearman_vs_middle_startle <- ggplot(correlate_microbiome, aes(log(Streptococcus+pseudocount), middle_average_VMax)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = paste("Spearman Correlation:", truncated_correlation),
       subtitle = paste("p-value:", truncated_p_value)) + 
  theme_cowplot(12)

correlation <- cor(correlate_microbiome$Streptococcus, 
                   correlate_microbiome$last_average_VMax, method = "spearman")
p_value <- cor.test(correlate_microbiome$Streptococcus, 
                    correlate_microbiome$last_average_VMax, method = "spearman")$p.value
truncated_correlation <- round(correlation, 4)
truncated_p_value <- round(p_value, 4)

spearman_vs_last_startle <- ggplot(correlate_microbiome, aes(log(Streptococcus+pseudocount),last_average_VMax)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = paste("Spearman Correlation:", truncated_correlation),
       subtitle = paste("p-value:", truncated_p_value)) + 
  theme_cowplot(12)

pseudocount <- 0.5
dat <- ggplot(correlate_microbiome, aes(Q22.x,log(Streptococcus+pseudocount), fill=Q22.x))+ 
  geom_boxplot(alpha=0.25)+
  scale_fill_viridis_d()+
  geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
  theme_cowplot(12) +
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title="Streptococcus ASV")+
  ylab("log(counts+0.5)")+
  xlab("")

plot_grid(dat, spearman_vs_first_startle,
          spearman_vs_middle_startle, spearman_vs_last_startle,
          labels=c("A","B","C","D"))


