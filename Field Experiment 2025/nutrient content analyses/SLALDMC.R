#---- Set-Up ----

rm(list = ls())

library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(readxl)
library(car)
library(performance)

sla <- read.csv("~/MastersThesis2025_git/Field Experiment 2025/nutrient content analyses/SLA_LDMC_summ25.csv")

#---- Data Cleaning ----

str(sla)
head(sla)
summary(sla)

sla <- sla %>% 
  select(-c(SLA, LDMC, notes)) %>%
  drop_na(fresh_weight) %>% 
  filter(round == 2) %>% 
  mutate(sla = leaf_area/fresh_weight, ldmc = dry_weight/fresh_weight, 
         trt = case_when(strip %in% c(1,3,5) ~ "b", strip %in% c(2,4,6) ~ "u"))

unique(sla$plant)

#---- Visualize ----

# regression of lead x fresh weight relationship not sure why. 

ggplot(sla, mapping = aes(x = leaf_area, y = fresh_weight)) + 
  geom_point()

# sla overall 

ggplot(sla, mapping = aes(x = as.factor(round), y = sla, color = trt)) + 
  geom_boxplot() + 
  theme_classic(base_size = 14) +
  labs(x = "Plant species", y = "Surface Leaf Area (SLA)", title = "SLA Across all Plants") +
  scale_color_manual(values = c("#E91E63","#4CAF50")) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#E8F5E9", color = NA),
    strip.text = element_text(face = "bold", color = "#1B5E20", size = 12),
    axis.text = element_text(color = "#424242"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# sla of all plant species 

ggplot(sla, mapping = aes(x = as.factor(round), y = sla, color = trt)) + 
  geom_boxplot() + 
  facet_wrap(~ plant, scales = "free_y") + 
  theme_classic(base_size = 14) +
  labs(x = "Plant species", y = "Surface Leaf Area (SLA)", title = "SLA Across Treatments") +
  scale_color_manual(values = c("#E91E63","#4CAF50")) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#E8F5E9", color = NA),
    strip.text = element_text(face = "bold", color = "#1B5E20", size = 12),
    axis.text = element_text(color = "#424242"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# ldmc overall  

ggplot(sla, mapping = aes(x = as.factor(round), y = ldmc, color = trt)) + 
  geom_boxplot() + 
  theme_classic(base_size = 14) +
  labs(x = "Plant species", y = "Leaf Dry Matter Content (LDMC)", title = "LDMC across all Plants") +
  scale_color_manual(values = c("#E91E63","#4CAF50")) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#E8F5E9", color = NA),
    strip.text = element_text(face = "bold", color = "#1B5E20", size = 12),
    axis.text = element_text(color = "#424242"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# ldmc of all plant species 

ggplot(sla, mapping = aes(x = as.factor(round), y = ldmc, fill = trt)) + 
  geom_boxplot() + 
  facet_wrap(~ plant, scales = "free_y") + 
  theme_classic(base_size = 14) +
  labs(x = "Plant species", y = "Leaf Dry Matter Content (LDMC)", title = "LDMC Across Treatments") +
  scale_color_manual(values = c("#E91E63","#4CAF50")) +
  scale_fill_manual(values = c("#E75480", "#00945C")) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#E8F5E9", color = NA),
    strip.text = element_text(face = "bold", color = "#1B5E20", size = 12),
    axis.text = element_text(color = "#424242"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14)
    ) +
  theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

#ggsave(
 # filename = "SLABvsU.jpeg",
#  plot = SLABvsU,
#  device = "jpeg",
#  dpi = 600,
#  width = 9,
#  height = 5,
#  units = "in"
# )

# ---- Models ----

# sla model 

ggplot(sla, mapping = aes(x = as.factor(round), y = sla, color = trt)) + 
  geom_boxplot() + 
  facet_wrap(~ plant) 

slam <- glmmTMB(sla ~ trt * plant + (1|strip), data = sla %>% 
                 filter(round == 2, plant != "dican"))

simulateResiduals(slam, plot = T)

summary(slam)
Anova(slam)
emmeans(slam, pairwise ~ trt|plant)

# ldmc model

ggplot(sla, mapping = aes(x = as.factor(round), y = ldmc, color = trt)) + 
  geom_boxplot() + 
  facet_wrap(~ plant) 

ldmcm <- glmmTMB(ldmc ~ trt * plant + (1|strip), data = sla %>% 
                 filter(round == 2, plant != "dican"))

plot(simulateResiduals(ldmcm))

summary(ldmcm)
Anova(ldmcm)
emmeans(ldmcm, pairwise ~ trt|plant)
emmeans(ldmcm, ~ trt|plant)



