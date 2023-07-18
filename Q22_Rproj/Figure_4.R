library(ggplot2)
library(rlang)
library(cowplot)
library(viridis)
library(tidyr)
library(dplyr)

#Read in startle data 
data<-read.csv("Q22_Behavior/Startle_PPI_Analysis - Analysis.csv", header=TRUE)
data$SLC_Genotype<-data$SLC
data$SLC_Genotype <- factor(data$SLC_Genotype, levels=c("WT", "HET", "MUT"))
data$Q22 <- factor(data$Q22, levels=c("WT","KO"))
data_slcwt <- data %>% filter(SLC_Genotype=="WT")

#Grab Ileum counts
metadata <- read.delim("Q22_Microbiome/starting_files/Q22_Metadata.tsv", header=TRUE)
metadata$SampleID <- gsub("-",".",metadata$SampleID)
metadata$SampleID <- paste0("X","", metadata$SampleID)
counts <- read.delim("Q22_Microbiome/starting_files/Q22_ASV.tsv", header = TRUE,row.names=1)

ileum_meta <- metadata %>% filter(Site=="ILE", SampleID %in% names(counts))
row.names(ileum_meta) <- ileum_meta$SampleID
ileum <- ileum_meta$SampleID
ileum_counts <- counts %>% select(all_of(ileum))

# Grab significantly diff abundant taxa 
significant_feature <- read.csv("Q22_Microbiome/differential_taxa/Ileum_ASV_Maaslin2_Sex_Q22_1-Litter/annotated_significant_results.csv")
streptococcus <- significant_feature$feature




generate_boxplots <- function(input_data, X, Y, min,max){
  data<-as.data.frame(input_data)
  
  ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{X}})) + 
    #geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    geom_boxplot(alpha=0.25)+
    scale_fill_viridis_d()+
    geom_point(size=1,position=position_jitter(width=0.25),alpha=0.5)+
    theme_cowplot(12) +
    ylim(min,max)+ 
    #facet_grid(~Sex)+    
    theme(legend.position = "none")
  
}
## First middle and last startle over time by Q22 
pos <- data_long %>% filter(Q22=="KO")
neg <- data_long %>% filter(Q22=="WT")

data_long_slcwt <- data_long %>% filter(SLC_Genotype=="WT")

posvneg <- ggplot(data_long, aes(x = Stage, y = Startle_Response, group = Q22, color = Q22)) +
  geom_point(aes(group = MouseID)) + 
  geom_line(aes(group = MouseID,color=Q22),size=1)+
  labs(x = "Time (minutes)", y = "FP_output") +
  scale_color_viridis_d()  +
  theme_cowplot(16) + 
  ggtitle("ASO FP output over time") + 
  theme(legend.position = "top",legend.justification = "center",legend.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))
posvneg 

posvneg_slcwt <- ggplot(data_long_slcwt, aes(x = Stage, y = Startle_Response, group = Q22, color = Q22)) +
  geom_point(aes(group = MouseID)) + 
  geom_line(aes(group = MouseID,color=Q22),size=1)+
  labs(x = "Time (minutes)", y = "FP_output") +
  scale_color_viridis_d()  +
  theme_cowplot(16) + 
  ggtitle("ASO FP output over time") + 
  theme(legend.position = "top",legend.justification = "center",legend.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))
posvneg_slcwt

## Startle PPI on full dataset --
summary(data$X70db_Percent_Inhibition)
db70 <-generate_boxplots(data, Q22, X70db_Percent_Inhibition,-30,85) +
  ggtitle("70 db Prepulse Inhibition")+
  theme(plot.title = element_text(hjust = 0.5))
db70

db75 <-generate_boxplots(data, Q22, X75db_Percent_Inhibition,-30,85) +
  ggtitle("75 db Prepulse Inhibition")+
  theme(plot.title = element_text(hjust = 0.5))
db75

db80 <-generate_boxplots(data, Q22, X75db_Percent_Inhibition,-30,85) +
  ggtitle("80 db Prepulse Inhibition")+
  theme(plot.title = element_text(hjust = 0.5))
db80

summary(data$First_average_VMax)
first_startle <-generate_boxplots(data, Q22, First_average_VMax, 0,1200) +
  ggtitle("First Startle")+
  theme(plot.title = element_text(hjust = 0.5))
first_startle

summary(data$middle_average_VMax)
middle_startle <-generate_boxplots(data, Q22, middle_average_VMax, 0,1000) +
  ggtitle("Middle Startle")+
  theme(plot.title = element_text(hjust = 0.5))
middle_startle

summary(data$last_average_VMax)
last_startle <-generate_boxplots(data, Q22, last_average_VMax, 0,1000) +
  ggtitle("last Startle")+
  theme(plot.title = element_text(hjust = 0.5))
last_startle

## Startle PPI on full dataset by SLC Genotype --
summary(data$X70db_Percent_Inhibition)
db70 <-generate_boxplots(data, SLC_Genotype, X70db_Percent_Inhibition,-30,85) +
  ggtitle("70 db Prepulse Inhibition")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~Q22)
db70

db75 <-generate_boxplots(data, SLC_Genotype, X75db_Percent_Inhibition,-30,85) +
  ggtitle("75 db Prepulse Inhibition")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~Q22)
db75

db80 <-generate_boxplots(data, SLC_Genotype, X75db_Percent_Inhibition,-30,85) +
  ggtitle("80 db Prepulse Inhibition")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~Q22)
db80

summary(data$First_average_VMax)
first_startle <-generate_boxplots(data, SLC_Genotype, First_average_VMax, 0,1200) +
  ggtitle("First Startle")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~Q22)
first_startle

summary(data$middle_average_VMax)
middle_startle <-generate_boxplots(data, SLC_Genotype, middle_average_VMax, 0,1000) +
  ggtitle("Middle Startle")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~Q22)
middle_startle

summary(data$last_average_VMax)
last_startle <-generate_boxplots(data, SLC_Genotype, last_average_VMax, 0,1000) +
  ggtitle("last Startle")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~Q22)
last_startle

## Startle PPI on only SLC WT mice
data_slcwt <- data %>% filter(SLC_Genotype=="WT")
summary(data_slcwt$X70db_Percent_Inhibition)
db70_slcwt <-generate_boxplots(data_slcwt, Q22, X70db_Percent_Inhibition,-36,85) +
  ggtitle("70 db")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="", y="% Prepulse Inhibition")
db70_slcwt

db75_slcwt <-generate_boxplots(data_slcwt, Q22, X75db_Percent_Inhibition,-36,85) +
  ggtitle("75 db")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="", y="% Prepulse Inhibition")
db75_slcwt 

db80_slcwt <-generate_boxplots(data_slcwt, Q22, X75db_Percent_Inhibition,-36,85) +
  ggtitle("80 db")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="", y="% Prepulse Inhibition")
db80_slcwt

summary(data$First_average_VMax)
startle_slcwt <-generate_boxplots(data_slcwt, Q22, First_average_VMax, 0,800) +
  ggtitle("First Startle")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="", y="Average VMax")
startle_slcwt

summary(data$middle_average_VMax)
middle_startle_slcwt <-generate_boxplots(data_slcwt, Q22, middle_average_VMax, 0,800) +
  ggtitle("Middle Startle")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="", y="Average VMax")
middle_startle_slcwt

last_startle_slcwt <-generate_boxplots(data_slcwt, Q22, last_average_VMax, 0,800) +
  ggtitle("Last Startle")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="", y="Average VMax")
last_startle_slcwt


plot_grid(startle_slcwt, middle_startle_slcwt, last_startle_slcwt,
          db70_slcwt, db75_slcwt,db80_slcwt, labels=c("A","B","C","D","E", "F"))


## Parametric Stats for SLC WT only --
lm1 <- lm(First_average_VMax ~ Sex+ Q22, data = data_slcwt)
summary(lm1)

wilcox.test(First_average_VMax~Q22,data_slcwt )

lm2 <- lm(X70db_Percent_Inhibition ~ Litter+ Sex + Q22, data = data_slcwt)
summary(lm2)

wilcox.test(X70db_Percent_Inhibition~Q22,data_slcwt )


lm3 <- lm(X75db_Percent_Inhibition ~ Litter+ Sex + Q22, data = data_slcwt)
summary(lm3)

wilcox.test(X75db_Percent_Inhibition~Q22,data_slcwt )

lm4 <- lm(X80db_Percent_Inhibition ~ Litter+ Sex + Q22, data = data_slcwt)
summary(lm4)

wilcox.test(X80db_Percent_Inhibition~Q22,data_slcwt )

lm5 <- lm(middle_average_VMax ~Sex +Q22, data = data_slcwt)
summary(lm5)
wilcox.test(middle_average_VMax~Q22,data_slcwt )


lm6 <- lm(last_average_VMax ~Sex + Q22, data = data_slcwt)
summary(lm6)
wilcox.test(last_average_VMax~Q22,data_slcwt )



# Save outputs -
sink("Buried_Food_Pellet_Stats.md")
cat("\n\nSummary for all data:\n")
print(summary(lm1))
cat("\n\nSummary for interaction:\n")
print(summary(interact))
cat("\n\nSummary for Tg Neg:\n")
print(summary(lm3))
cat("\n\nSummary for Tg Pos:\n")
print(summary(lm2))
cat("\n\nSummary for Females Tg Neg:\n")
print(summary(lm4))
cat("\n\nSummary for Females Tg Pos:\n")
print(summary(lm5))
cat("\n\nSummary for Males Tg Neg:\n")
print(summary(lm6))
cat("\n\nSummary for Males Tg Pos:\n")
print(summary(lm7))
sink()