library(Maaslin2)
library(funrar)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)
#library(circlize)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Q22_Rproj/Q22_PWY_Maaslin2.R")

metadata <- read.delim("Q22_Microbiome/starting_files/Q22_Metadata.tsv", header=TRUE)
metadata$SampleID <- gsub("-",".",metadata$SampleID)
metadata$SampleID <- paste0("X","", metadata$SampleID)
metadata$Genotype <- revalue(metadata$Q22, replace = c("KO" = "Q22","WT"="WT"))
counts <- read.delim("Q22_Microbiome/starting_files/Q22_ASV.tsv", header = TRUE,row.names=1)
pathway <- read.delim("Q22_Microbiome/picrust2/export_pathway_abundance/annotated_pathway.tsv", header = TRUE,row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(pathway, "feature") %>% select(c("feature", "description"))
counts <- counts %>% select(-c("taxonomy"))
pathway <- pathway %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 3000]

## Grab samples in the pathway table that are present in counts --
names(pathway)
pathway <- pathway %>% select(c(names(counts)))

### Ileum ---

input_metadata <- metadata
df_input_data <- pathway
samples <- input_metadata %>% filter(Site =="ILE", SampleID %in% names(df_input_data)) %>% pull(SampleID)

df_input_data <- df_input_data[, samples]

row.names(input_metadata) <- input_metadata$SampleID
target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)


df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Genotype <- factor(df_input_metadata$Genotype, levels=c("WT","Q22"))
df_input_metadata$Litter <- factor(df_input_metadata$Litter)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
sapply(df_input_metadata,levels)

fit_data = Maaslin2(input_data=df_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0("Q22_Microbiome/picrust2/Ile_PWY_Maaslin2_Sex_Genotype_1-Litter"), 
                    fixed_effects = c("Sex", "Genotype"), normalization = "TSS", 
                    random_effects = c("Litter"),
                    #reference = c("Genotype,WT", "Site,Distal_Colon"),
                    min_prevalence = 0.15,
                    transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)

fit_data = Maaslin2(input_data=df_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0("Q22_Microbiome/picrust2/Ile_PWY_Maaslin2_Sex_Genotype"), 
                    fixed_effects = c("Sex", "Genotype"), normalization = "TSS", 
                    #reference = c("Genotype,WT", "Site,Distal_Colon"),
                    min_prevalence = 0.15,
                    transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)

## Add annotations to Maaslin2 results --
# Ileum

data<-read.delim("Q22_Microbiome/picrust2/Ile_PWY_Maaslin2_Sex_Genotype_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25) %>% filter(metadata=="Genotype")
data$feature <- gsub("\\.","-",data$feature)
data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Q22_Microbiome/picrust2/Ile_PWY_Maaslin2_Sex_Genotype_1-Litter/significant_results.csv")

data<-read.delim("Q22_Microbiome/picrust2/Ile_PWY_Maaslin2_Sex_Q22/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25) %>% filter(metadata=="Q22")
data$feature <- gsub("\\.","-",data$feature)
data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Q22_Microbiome/picrust2/Ile_PWY_Maaslin2_Sex_Q22/significant_results.csv")

### Colon ---

input_metadata <- metadata
df_input_data <- pathway
samples <- input_metadata %>% filter(Site =="COL", SampleID %in% names(df_input_data)) %>% pull(SampleID)

df_input_data <- df_input_data[, samples]

row.names(input_metadata) <- input_metadata$SampleID
target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)


df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Q22 <- factor(df_input_metadata$Q22, levels=c("WT","Q22"))
df_input_metadata$Litter <- factor(df_input_metadata$Litter)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
sapply(df_input_metadata,levels)

fit_data = Maaslin2(input_data=df_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0("Q22_Microbiome/picrust2/COL_PWY_Maaslin2_Sex_Q22_1-Litter"), 
                    fixed_effects = c("Sex", "Q22"), normalization = "TSS", 
                    random_effects = c("Litter"),
                    #reference = c("Genotype,WT", "Site,Distal_Colon"),
                    min_prevalence = 0.15,
                    transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)

fit_data = Maaslin2(input_data=df_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0("Q22_Microbiome/picrust2/COL_PWY_Maaslin2_Sex_Q22"), 
                    fixed_effects = c("Sex", "Q22"), normalization = "TSS", 
                    #reference = c("Genotype,WT", "Site,Distal_Colon"),
                    min_prevalence = 0.15,
                    transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)

## Add annotations to Maaslin2 results --
# Colon

data<-read.delim("Q22_Microbiome/picrust2/COL_PWY_Maaslin2_Sex_Q22_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25) %>% filter(metadata=="Q22")
data$feature <- gsub("\\.","-",data$feature)
data<- (merge(data, annotation, by = 'feature'))
write.table(data, "Q22_Microbiome/picrust2/COL_PWY_Maaslin2_Sex_Q22_1-Litter/significant_results.tsv")

data<-read.delim("Q22_Microbiome/picrust2/COL_PWY_Maaslin2_Sex_Q22/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25) %>% filter(metadata=="Q22")
data$feature <- gsub("\\.","-",data$feature)
data<- (merge(data, annotation, by = 'feature'))
write.table(data, "Q22_Microbiome/picrust2/COL_PWY_Maaslin2_Sex_Q22/significant_results.tsv")


### Visualization of the Results ---

create_bar_plot <- function(res_plot, title) {
  y <- tapply(res_plot$coef, res_plot$description, function(y) mean(y))
  y <- sort(y, FALSE)
  res_plot$description <- factor(as.character(res_plot$description), levels = names(y))
  cols <- c(WT = "black", Q22 = "firebrick")
  res_plot %>%
    arrange(coef) %>%
    filter(qval < 0.25, abs(coef) > 0) %>%
    ggplot(aes(coef, description, fill = site)) +
    geom_bar(stat = "identity") +
    theme_cowplot(12) +
    theme(axis.text.y = element_text(face = "bold")) +
    scale_fill_manual(values = cols) +
    labs(x = "Effect size", y = "", fill = "") +
    theme(legend.position = "none") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
}

## Model 1: PWY ~ Site + Sex + Genotype + (1|mouseID) --
# Ileum -
data<-read.csv("Q22_Microbiome/picrust2/Ile_PWY_Maaslin2_Sex_Genotype_1-Litter/significant_results.csv", 
               header=TRUE, row.names=1)

res_plot <- data %>%
    filter(value == "Q22") %>%
    unique() %>%
    mutate(site = ifelse(coef < 0, "WT", "Q22"))
pwy_mut <- create_bar_plot(res_plot, "WT vs Q22 q<0.25")
res_plot <- data %>%
    filter(value == "HET") %>%
    unique() %>%
    mutate(site = ifelse(coef < 0, "WT", "HET"))
pwy_het <- create_bar_plot(res_plot, "Mucosal colon: WT vs HET q<0.25")
plot_grid(pwy_mut, pwy_het, labels = c("A", "B"))

# Luminal Colon -
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Luminal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
annotation <- read.delim("Trios/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

res_plot <- data %>%
  filter(value == "MUT") %>%
  unique() %>%
  mutate(site = ifelse(coef < 0, "WT", "MUT"))
pwy_mut <- create_bar_plot(res_plot, "Mucosal Colon: WT vs MUT q<0.25")
res_plot <- data %>%
  filter(value == "HET") %>%
  unique() %>%
  mutate(site = ifelse(coef < 0, "WT", "HET"))
pwy_het <- create_bar_plot(res_plot, "Mucosal colon: WT vs HET q<0.25")
plot_grid(pwy_mut, pwy_het, labels = c("A", "B"))

## Model 2: PWY ~ Sequencing_Run + Site + Sex + Genotype + (1|mouseID) --
# Mucosal Colon -
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Mucosal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
annotation <- read.delim("Trios/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")
higher_classification <- read.csv("MetaCyc_pathwaynames_Key.csv", row.names=1, header=TRUE)
higher_classification$feature <- row.names(higher_classification)
data <- merge(data,higher_classification, by="feature")

res_plot <- data %>%
  filter(value == "MUT") %>%
  unique() %>%
  mutate(site = ifelse(coef < 0, "WT", "MUT"))
pwy_mut <- create_bar_plot(res_plot, "Mucosal Colon: WT vs MUT q<0.25")
res_plot <- data %>%
  filter(value == "HET") %>%
  unique() %>%
  mutate(site = ifelse(coef < 0, "WT", "HET"))
pwy_het <- create_bar_plot(res_plot, "Mucosal colon: WT vs HET q<0.25")
plot_grid(pwy_mut, pwy_het, labels = c("A", "B"))

# Luminal Colon -
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Luminal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
annotation <- read.delim("Trios/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

res_plot <- data %>%
  filter(value == "MUT") %>%
  unique() %>%
  mutate(site = ifelse(coef < 0, "WT", "MUT"))
pwy_mut <- create_bar_plot(res_plot, "Luminal Colon: WT vs MUT q<0.25")
res_plot <- data %>%
  filter(value == "HET") %>%
  unique() %>%
  mutate(site = ifelse(coef < 0, "WT", "HET"))
pwy_het <- create_bar_plot(res_plot, "Luminal Colon: WT vs HET q<0.25")
plot_grid(pwy_mut, pwy_het, labels = c("A", "B"))

### Visualizing results as circle plots ---

# Mucosal Colon MUT enriched 
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Mucosal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0)
annotation <- read.delim("Trios/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

#Classification strategy no 1 
higher_classification <- read.csv("MetaCyc_pathwaynames_Key.csv", row.names=1, header=TRUE)
higher_classification$feature <- row.names(higher_classification)
data <- merge(data,higher_classification, by="feature")
data <- rename(data, replace = c("L4A" = "classification"))
data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0) %>%
  select(c("feature", "classification","coef"))

#Classification strategy no 2 
# split the paths column by "|"
higher_classification <- read.delim("Huttenhower_MetaCyc_Hierarchy.txt",header=TRUE)
df <- higher_classification
df_split <- strsplit(df$annotation, "\\|")
df_new <- data.frame(do.call(rbind, df_split))
df_new$feature <- higher_classification$feature
df_new$feature <- gsub("-",".",df_new$feature)
data <- merge(data,df_new, by="feature")
data <- data %>% 
  select(c("feature", "X1","coef"))

mat <- data 
?circos.text()
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.25,
              cex=1)
}, bg.border = NA) 
circos.clear()

# Mucosal Colon MUT depleted 
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Mucosal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef<0)
annotation <- read.delim("Trios/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

# Classification strategy no 1
higher_classification <- read.csv("MetaCyc_pathwaynames_Key.csv", row.names=1, header=TRUE)
higher_classification$feature <- row.names(higher_classification)
data <- merge(data,higher_classification, by="feature")
data <- rename(data, replace = c("L4A" = "classification"))
data <- data %>% 
  select(c("feature", "classification","coef"))
data<- remove_missing(data)

#Classification strategy no 2 
# split the paths column by "|"
higher_classification <- read.delim("Huttenhower_MetaCyc_Hierarchy.txt",header=TRUE)
df <- higher_classification
df_split <- strsplit(df$annotation, "\\|")
df_new <- data.frame(do.call(rbind, df_split))
df_new$feature <- higher_classification$feature
df_new$feature <- gsub("-",".",df_new$feature)
data <- merge(data,df_new, by="feature")
data <- data %>% 
  select(c("feature", "X1","coef"))
data <- unique(data)

mat <- data 
?circos.text()
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()

# Luminal Colon MUT enriched 
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Luminal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef>0)
annotation <- read.delim("Trios/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

# split the paths column by "|"
higher_classification <- read.delim("Huttenhower_MetaCyc_Hierarchy.txt",header=TRUE)
df <- higher_classification
df_split <- strsplit(df$annotation, "\\|")
df_new <- data.frame(do.call(rbind, df_split))
df_new$feature <- higher_classification$feature
df_new$feature <- gsub("-",".",df_new$feature)

# Merge higher_classification column
data <- merge(data,df_new, by="feature")
data <- data %>% 
  select(c("description", "X1","coef"))

mat <- data 
?circos.text()
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.25,
              cex=1)
}, bg.border = NA) 
circos.clear()

# Luminal Colon MUT depleted 
data<-read.table("Trios/differential_Pathway/PICRUST2_PWY_Luminal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% 
  filter(value=="MUT") %>% 
  filter(coef<0)
annotation <- read.delim("Trios/differential_Pathway/annotated_pwy.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub("-", ".", annotation$feature)
data <- merge(data,annotation, by="feature")

# Classification strategy no 1
higher_classification <- read.csv("MetaCyc_pathwaynames_Key.csv", row.names=1, header=TRUE)
higher_classification$feature <- row.names(higher_classification)
data <- merge(data,higher_classification, by="feature")
data <- rename(data, replace = c("L4A" = "classification"))
data <- data %>% 
  select(c("feature", "classification","coef"))
data<- remove_missing(data)

#Classification strategy no 2 
# split the paths column by "|"
higher_classification <- read.delim("Huttenhower_MetaCyc_Hierarchy.txt",header=TRUE)
df <- higher_classification
df_split <- strsplit(df$annotation, "\\|")
df_new <- data.frame(do.call(rbind, df_split))
df_new$feature <- higher_classification$feature
df_new$feature <- gsub("-",".",df_new$feature)
data <- merge(data,df_new, by="feature")
plot <- data %>% 
  select(c("feature", "X1","coef"))
plot <- data %>% 
  select(c("feature", "X2","coef"))

mat <- plot
?circos.text()
circos.clear()
dev.new(width=10,height=10)
chordDiagram(mat,annotationTrack = "grid",preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
obj <- circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = 0.4,
              cex=1)
}, bg.border = NA) 
circos.clear()