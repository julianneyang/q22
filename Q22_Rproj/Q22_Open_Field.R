library(ggplot2)
library(rlang)
library(cowplot)
library(viridis)
library(tidyr)
library(dplyr)

data<-read.csv("../Q22 Open Field Analysis - Sheet1.csv", header=TRUE)
data$SLC_Genotype<-data$SLC
data$SLC_Genotype <- factor(data$SLC_Genotype, levels=c("WT", "HET", "MUT"))
data$Q22 <- factor(data$Q22, levels=c("WT","KO"))
generate_boxplots <- function(input_data, X, Y, min,max){
  data<-as.data.frame(input_data)
  
  ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{X}})) + 
    #geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    geom_boxplot(alpha=0.25)+
    scale_fill_viridis_d()+
    geom_point(size=1,position=position_jitter(width=0.25),alpha=0.5)+
    theme_cowplot(12) +
    #ylim(min,max)+ 
    #facet_grid(~Sex)+    
    theme(legend.position = "none")
  
}

## Open Field on full dataset --
distance <-generate_boxplots(data, Q22, Distance,-30,85) +
  ggtitle("Distance in Open Field")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="", y="Distance (cm)")
distance

centertime <-generate_boxplots(data, Q22, Center_Time,-30,85) +
  ggtitle("Center Time")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="", y="Time (s)")
centertime

plot_grid(distance, centertime, labels=c("A","B"))

## Open Field on full dataset by SLC Genotype--
distance <-generate_boxplots(data, SLC_Genotype, Distance,-30,85) +
  ggtitle("Distance in Open Field")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="", y="Distance (cm)")+
  facet_wrap(~Q22)
distance

centertime <-generate_boxplots(data, SLC_Genotype, Center_Time,-30,85) +
  ggtitle("Center Time")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="", y="Time (s)")+
  facet_wrap(~Q22)
centertime

plot_grid(distance, centertime, labels=c("A","B"))

## Startle PPI on only SLC WT mice
data_slcwt <- data %>% filter(SLC_Genotype=="WT")
distance_slcwt <-generate_boxplots(data_slcwt, Q22, Distance,-30,85) +
  ggtitle("Distance in Open Field")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="", y="Distance (cm)")
distance_slcwt

centertime_slcwt <-generate_boxplots(data_slcwt, Q22, Center_Time,-30,85) +
  ggtitle("Center Time")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="", y="Time (s)")
centertime_slcwt

plot_grid(distance_slcwt, centertime_slcwt, 
          olm_slcwt_training,olm_slcwt,
          labels=c("A","B","C","D"),
          nrow=2)

## Stratified by SLC Genotype --
b<-generate_boxplots(data, SLC_Genotype, Total_Time,0,51)+
  facet_wrap(~Q22)
b

# Stratified by SLC Genotype and Sex 
e<-generate_boxplots(data, SLC_Genotype, Total_Time,0,51)+
  facet_wrap(Sex~ASO_Tg)
e

## Parametric Stats for SLC WT only --
lm1 <- lm(Distance ~ Sex+ Q22, data = data_slcwt)
summary(lm1)

wilcox.test(Distance~Q22,data_slcwt )

lm1 <- lm(Center_Time~ Sex+ Q22, data = data_slcwt)
summary(lm1)

wilcox.test(Center_Time~Q22,data_slcwt )



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
