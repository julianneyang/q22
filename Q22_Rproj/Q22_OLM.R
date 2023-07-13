library(ggplot2)
library(rlang)
library(cowplot)
library(viridis)
library(tidyr)
library(dplyr)

data<-read.csv("../OLM_Analysis - OLM_analysis.csv", header=TRUE)
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
data_testing <- data %>% filter(Day=="Testing")
olm <-generate_boxplots(data_testing, Q22, Discrimination_Ratio,-30,85) +
  ggtitle("Testing")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="", y="Discrimination Ratio")
olm

data_training <- data %>% filter(Day=="Training")
olm_training <-generate_boxplots(data_training, Q22, Discrimination_Ratio,-30,85) +
  ggtitle("Training")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="", y="Discrimination Ratio")
olm_training

olm_long <- ggplot(data=data,aes(x=Day,y=Discrimination_Ratio, fill=Q22))+ 
  scale_fill_viridis_d() +
  geom_boxplot(alpha=0.25)+
  geom_point(aes(group = Animal)) + 
  geom_line(aes(group = Animal,color=Q22),size=1)+
  theme_cowplot(12)+
  facet_wrap(~Q22)

plot_grid(olm_training,olm,
          olm_long, labels=c("A","B","C"),
          ncol=2, nrow = 2,
          rel_widths = as.numeric(0.5,0.5,1)) 

## Open Field on full dataset by SLC Genotype--
data_testing <- data %>% filter(Day=="Testing")
olm <-generate_boxplots(data_testing, SLC_Genotype, Discrimination_Ratio,-30,85) +
  ggtitle("Testing")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="", y="Discrimination Ratio")+
  facet_wrap(~Q22)
olm

data_training <- data %>% filter(Day=="Training")
olm_training <-generate_boxplots(data_training, SLC_Genotype, Discrimination_Ratio,-30,85) +
  ggtitle("Training")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="", y="Discrimination Ratio")+
  facet_wrap(~Q22)
olm_training

olm_long <- ggplot(data=data,aes(x=Day,y=Discrimination_Ratio, fill=SLC_Genotype))+ 
  scale_fill_viridis_d() +
  geom_boxplot(alpha=0.25)+
  geom_point(aes(group = Animal)) + 
  geom_line(aes(group = Animal,color=SLC_Genotype),size=1)+
  theme_cowplot(12)+
  facet_wrap(~Q22)
olm_long

plot_grid(olm_training,olm,
          olm_long, labels=c("A","B","C"),
          ncol=2, nrow = 2,
          rel_widths = as.numeric(0.5,0.5,1)) 


## Startle PPI on only SLC WT mice
data_slcwt <- data %>% filter(SLC_Genotype=="WT")
data_testing_slcwt <- data_slcwt %>% filter(Day=="Testing")
olm_slcwt <-generate_boxplots(data_testing_slcwt, Q22, Discrimination_Ratio,-30,85) +
  ggtitle("Testing")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="", y="Discrimination Ratio")
olm_slcwt

data_training_slcwt <- data_slcwt %>% filter(Day=="Training")
olm_slcwt_training <-generate_boxplots(data_training_slcwt, Q22, Discrimination_Ratio,-30,85) +
  ggtitle("Training")+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="", y="Discrimination Ratio")
olm_slcwt_training

olm_long_slcwt <- ggplot(data=data_slcwt,aes(x=Day,y=Discrimination_Ratio, fill=Q22))+ 
  scale_fill_viridis_d() +
  geom_boxplot(alpha=0.25)+
  geom_point(aes(group = Animal)) + 
  geom_line(aes(group = Animal,color=Q22),size=1)+
  theme_cowplot(12)+
  facet_wrap(~Q22)

plot_grid(olm_slcwt_training,olm_slcwt,
          olm_long_slcwt, labels=c("A","B","C"),
          ncol=2, nrow = 2,
          rel_widths = as.numeric(0.5,0.5,1)) 

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
