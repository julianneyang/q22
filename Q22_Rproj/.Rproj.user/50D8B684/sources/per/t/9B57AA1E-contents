library(ggplot2)
library(dplyr)
library(cowplot)
library(nlme)
library(tidyr)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/pdbehavior/")
data <- read.csv("Analysis_Files/SMT/Fecal Pellet Output (SMT) - Total_FP_output.csv",header=TRUE)

data_long <- pivot_longer(data, 
                          cols = starts_with("X"), 
                          names_to = "timepoint", 
                          values_to = "FP_output")

data_long$timepoint <- as.integer(stringr::str_extract(data_long$timepoint, "\\d+"))
data_long$SLC_Genotype <- factor(data_long$SLC_Genotype, levels=c("WT", "HET", "MUT"))
data_long$timepoint <- factor(data_long$timepoint)

## Plots --
one_hour <- data_long %>% filter(timepoint==60)
generate_boxplots(one_hour,SLC_Genotype,FP_output,0, 20) + 
  facet_wrap(~Sex)+
  ggtitle("Tg Positive")

generate_boxplots(pos,SLC_Genotype,FP_output,0, 30) + 
  facet_wrap(Sex~timepoint,nrow=2)+
  ggtitle("Tg Positive by Sex")

generate_boxplots(neg,SLC_Genotype,FP_output,0, 30) + 
  facet_wrap(~timepoint)+
  ggtitle("Tg Negative")

generate_boxplots(neg,SLC_Genotype,FP_output,0, 30) + 
  facet_wrap(Sex~timepoint)+
  ggtitle("Tg Negative by Sex")


# Calculate the mean and standard error for each group

df_summary <- data_long %>%
  group_by(SLC_Genotype, timepoint) %>%
  summarise(mean = mean(FP_output), 
            sd = sd(FP_output),
            se = sd / sqrt(n()))

# Plot the graph with error bars
slc_genotype <- ggplot(df_summary, aes(x = timepoint, y = mean, group = SLC_Genotype, color = SLC_Genotype)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(x = "Time (minutes)", y = "FP_output") +
  scale_color_viridis_d()  +
  theme_cowplot(16) + 
  ggtitle("SMT: FP output over time") + 
  theme(legend.position = "top", legend.justification="center",legend.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5)) 


# Calculate the mean and standard error for each group -- Females
female <- data_long %>% filter(Sex=="Female")
female_df_summary <- female %>%
  group_by(SLC_Genotype, timepoint) %>%
  summarise(mean = mean(FP_output), 
            sd = sd(FP_output),
            se = sd / sqrt(n()))

# Plot the graph with error bars
femaleplot <- ggplot(female_df_summary, aes(x = timepoint, y = mean, group = SLC_Genotype, color = SLC_Genotype)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(x = "Time (minutes)", y = "FP_output") +
  scale_color_viridis_d()  +
  theme_cowplot(16) + 
  ggtitle("SMT Females")+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))


# Calculate the mean and standard error for each group -- Males
male <- data_long %>% filter(Sex=="Male")
male_df_summary <- male %>%
  group_by(SLC_Genotype, timepoint) %>%
  summarise(mean = mean(FP_output), 
            sd = sd(FP_output),
            se = sd / sqrt(n()))

# Plot the graph with error bars
male_plot <- ggplot(male_df_summary, aes(x = timepoint, y = mean, group = SLC_Genotype, color = SLC_Genotype)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(x = "Time (minutes)", y = "FP_output") +
  scale_color_viridis_d()  +
  theme_cowplot(16) + 
  ggtitle("SMT Males") +
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

# Calculate the mean and standard error for each group -- Females
female <- data_long %>% filter(Sex=="Female")
female_df_summary <- female %>%
  group_by(SLC_Genotype, timepoint) %>%
  summarise(mean = mean(FP_output), 
            sd = sd(FP_output),
            se = sd / sqrt(n()))

# Plot the graph with error bars
femaleplot <- ggplot(female_df_summary, aes(x = timepoint, y = mean, group = SLC_Genotype, color = SLC_Genotype)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(x = "Time (minutes)", y = "FP_output") +
  scale_color_viridis_d()  +
  theme_cowplot(16) + 
  ggtitle("SMT Females")+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))


# Calculate the mean and standard error for each group -- Males
male <- data_long %>% filter(Sex=="Male")
male_df_summary <- male %>%
  group_by(SLC_Genotype, timepoint) %>%
  summarise(mean = mean(FP_output), 
            sd = sd(FP_output),
            se = sd / sqrt(n()))

# Plot the graph with error bars
male_plot <- ggplot(male_df_summary, aes(x = timepoint, y = mean, group = SLC_Genotype, color = SLC_Genotype)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(x = "Time (minutes)", y = "FP_output") +
  scale_color_viridis_d()  +
  theme_cowplot(16) + 
  ggtitle("SMT Males") +
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))


# Calculate the mean and standard error for each group -- Females
female <- data_long %>% filter(Sex=="Female")
female_df_summary <- female %>%
  group_by(Cohort,SLC_Genotype, timepoint) %>%
  summarise(mean = mean(FP_output), 
            sd = sd(FP_output),
            se = sd / sqrt(n()))

# Plot the graph with error bars
femaleplot <- ggplot(female_df_summary, aes(x = timepoint, y = mean, group = SLC_Genotype, color = SLC_Genotype)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(x = "Time (minutes)", y = "FP_output") +
  scale_color_viridis_d()  +
  theme_cowplot(16) + 
  ggtitle("SMT Females")+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~Cohort)


# Calculate the mean and standard error for each group -- Males
male <- data_long %>% filter(Sex=="Male")
male_df_summary <- male %>%
  group_by(Cohort,SLC_Genotype, timepoint) %>%
  summarise(mean = mean(FP_output), 
            sd = sd(FP_output),
            se = sd / sqrt(n()))

# Plot the graph with error bars
male_plot <- ggplot(male_df_summary, aes(x = timepoint, y = mean, group = SLC_Genotype, color = SLC_Genotype)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(x = "Time (minutes)", y = "FP_output") +
  scale_color_viridis_d()  +
  theme_cowplot(16) + 
  ggtitle("SMT Males") +
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~Cohort)

## Final Figure --
plot_grid(slc_genotype,femaleplot,male_plot,nrow=1)

### Longitudinal Stats ---
nonpara_output <- kruskal.test(FP_output~SLC_Genotype, data=data_long)
print(nonpara_output)

lm1<- lme(fixed= FP_output ~ timepoint +Sex + SLC_Genotype, random = ~1|Mouse, data=data_long)
summary(lm1)

lm1_time_ASO<- lme(fixed= FP_output ~ timepoint*SLC_Genotype + Sex, random = ~1|Mouse, data=data_long)
summary(lm1_time_ASO)

males <- data_long %>% filter(Sex=="Male")
females <- data_long %>% filter(Sex=="Female")

lm4 <- lme(fixed= FP_output ~ timepoint*SLC_Genotype, random = ~1|Mouse, data=males)
summary(lm4)
lm5 <- lme(fixed= FP_output ~ timepoint*SLC_Genotype, random = ~1|Mouse, data=females)
summary(lm5)



# Save outputs -
sink("FP_output_Stats.md")
cat("\n\nSummary for all data:\n")
print(summary(lm1))
cat("\n\nSummary for all data, non parametric:\n")
print((nonpara_output))
cat("\n\nSummary for all data, time*ASO:\n")
print(summary(lm1_time_ASO))
cat("\n\nSummary for all data, SLC*ASO:\n")
print(summary(lm1_SLC_ASO))
cat("\n\nSummary for Tg Pos:\n")
print(summary(lm2))
cat("\n\nSummary for Tg Neg:\n")
print(summary(lm3))
cat("\n\nSummary for Males Tg Pos:\n")
print(summary(lm4))
cat("\n\nSummary for Males Tg Neg:\n")
print(summary(lm5))
cat("\n\nSummary for Males Tg Pos:\n")
print(summary(lm6))
cat("\n\nSummary for Males Tg Neg:\n")
print(summary(lm7))
sink()

### Cross-Sectional Stats ---
lm1 <- lm(FP_output ~ ASO_Tg,  data=data_long)
summary(lm1)
kruskal.test(FP_output ~ ASO_Tg,  data=data_long)

lm1 <- lm(FP_output ~ Sex + ASO_Tg +SLC_Genotype,  subset(data_long, timepoint == 60))
summary(lm1)

lm1 <- lm(FP_output ~ Sex +SLC_Genotype,  subset(data_long, timepoint == 60 & ASO_Tg=="Positive"))
summary(lm1)

lm1 <- lm(FP_output ~ Sex +SLC_Genotype,  subset(data_long, timepoint == 60 & ASO_Tg=="Negative"))
summary(lm1)

lm1 <- lm(FP_output ~SLC_Genotype,  subset(data_long, timepoint == 60 & ASO_Tg=="Positive" & Sex =="Female"))
summary(lm1)

lm1 <- lm(FP_output ~ Sex +SLC_Genotype,  subset(data_long, timepoint == 60 & ASO_Tg=="Negative" ))
summary(lm1)