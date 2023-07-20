
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(nlme)
library(ggpubr)
setwd("/home/julianne/Documents/q22/") # CHANGE to the directory containing the fastq files
here::i_am("Q22_Rproj/Q22_Alpha_Diversity.R")

#Colon 
otus<- readr::read_delim("Q22_Microbiome/alpha_diversity/alpha_COL_Q22_ASV/otus_dir/alpha-diversity.tsv")
row.names(otus) <- otus$...1
shannon<-readr::read_delim("Q22_Microbiome/alpha_diversity/alpha_COL_Q22_ASV/shannon_dir/alpha-diversity.tsv")
row.names(shannon) <- shannon$...1
chao1<-readr::read_delim("Q22_Microbiome/alpha_diversity/alpha_COL_Q22_ASV/chao1_dir/alpha-diversity.tsv")
row.names(chao1) <- chao1$...1


data<- merge(otus,shannon, by="...1")
data<- merge(data,chao1, by="...1")
data$SampleID <- data$...1

metadata<- read.delim("Q22_Microbiome/starting_files/Q22_Metadata.tsv")
metadata$Genotype <- revalue(metadata$Q22, replace = c("KO" = "Q22","WT"="WT"))
colon_data_meta <- merge(data,metadata, by="SampleID")


#Ileum 
otus<- readr::read_delim("Q22_Microbiome/alpha_diversity/alpha_ILE_Q22_ASV/otus_dir/alpha-diversity.tsv")
row.names(otus) <- otus$...1
shannon<-readr::read_delim("Q22_Microbiome/alpha_diversity/alpha_ILE_Q22_ASV/shannon_dir/alpha-diversity.tsv")
row.names(shannon) <- shannon$...1
chao1<-readr::read_delim("Q22_Microbiome/alpha_diversity/alpha_ILE_Q22_ASV/chao1_dir/alpha-diversity.tsv")
row.names(chao1) <- chao1$...1


data<- merge(otus,shannon, by="...1")
data<- merge(data,chao1, by="...1")
data$SampleID <- data$...1

metadata<- read.delim("Q22_Microbiome/starting_files/Q22_Metadata.tsv")
metadata$Genotype <- revalue(metadata$Q22, replace = c("KO" = "Q22","WT"="WT"))
ileum_data_meta <- merge(data,metadata, by="SampleID")

### Function for plotting alpha diversity ---
generate_adiv_plots <- function(input_data, X, Y, min, max){
  #read in files
  data <- as.data.frame(input_data)

  #declare order of variables
  data$Genotype <- factor(data$Genotype, levels=c("WT", "Q22"))
  #graph plot
  ggplot(data=data,aes(x={{X}},y={{Y}}, fill=Q22)) + 
    geom_boxplot(alpha=0.25)+
    scale_fill_viridis_d()+
    #geom_line(aes(group = MouseID,color=Genotype),size=1)+
    geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    #ylim(min,max) +
    theme(legend.position = "none")
  
}

### Make and store plots ---
compare <-c("WT","Q22")

adiv_colon_shannon<- generate_adiv_plots(colon_data_meta, Genotype, shannon_entropy, 2, 7) +
  #stat_compare_means(comparisons = compare,method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)+
  ggtitle("Colon")+
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5))
adiv_colon_shannon

adiv_colon_otus<- generate_adiv_plots(colon_data_meta, Genotype, observed_features, 0, 300) +
  #stat_compare_means(comparisons = compare,method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)+
  ggtitle("Colon")+
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5))
adiv_colon_otus

adiv_colon_chao1<- generate_adiv_plots(colon_data_meta, Genotype, chao1, 0, 300) +
  stat_compare_means(comparisons = compare,method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)+
  ggtitle("chao1")+
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5))
adiv_colon_chao1


adiv_ileum_shannon<- generate_adiv_plots(ileum_data_meta, Genotype, shannon_entropy, 2, 7) +
  stat_compare_means(comparisons = compare,method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)+
  ggtitle("Ileum")+
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5))
adiv_ileum_shannon

adiv_ileum_otus<- generate_adiv_plots(ileum_data_meta, Genotype, observed_features, 0, 300) +
  stat_compare_means(comparisons = compare,method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)+
  ggtitle("Ileum")+
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5))
adiv_ileum_otus

adiv_ileum_chao1<- generate_adiv_plots(ileum_data_meta, Genotype, chao1, 0, 300) +
  stat_compare_means(comparisons = compare,method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)+
  ggtitle("chao1")+
  theme(plot.title = element_text(hjust = 0.5))
adiv_ileum_chao1

plot_grid(adiv_ileum_shannon, adiv_ileum_otus, adiv_ileum_chao1, labels=c("A","B","C"), nrow=1)

### Alpha Diversity Stats ---
data_meta$Genotype <-factor(data_meta$Genotype, levels=c("WT", "HET","MUT"))
output <- lme(fixed= shannon_entropy ~ Sequencing_Run + Sex + Site+ Genotype, random = ~1|MouseID, data=data_meta)
summary(output)
output <- lme(fixed= observed_features ~ Sequencing_Run +Sex+ Site+ Genotype, random = ~1|MouseID, data=data_meta)
summary(output)
output <- lme(fixed= chao1 ~ Sequencing_Run +Sex+ Site+ Genotype, random = ~1|MouseID, data=data_meta)
summary(output)

output <- lme(fixed= shannon_entropy ~ Sequencing_Run + Site+ Sex*Genotype, random = ~1|MouseID, data=data_meta)
summary(output)
output <- lme(fixed= observed_features ~ Sequencing_Run +Sex+ Site+ Sex*Genotype, random = ~1|MouseID, data=data_meta)
summary(output)
output <- lme(fixed= chao1 ~ Sequencing_Run +Sex+ Site+ Sex*Genotype, random = ~1|MouseID, data=data_meta)
summary(output)
