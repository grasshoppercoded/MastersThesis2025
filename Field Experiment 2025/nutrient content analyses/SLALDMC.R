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

sla <- read.csv("SLA_LDMC_summ25.csv")

#---- Data Cleaning ----

str(sla)
head(sla)
summary(sla)

sla <- sla %>% 
  select(-c(SLA, LDMC, notes)) %>%
  drop_na(fresh_weight) %>% 
  filter(plant == "wide") %>% 
  filter(round == 2) %>% 
  mutate(sla = leaf_area/fresh_weight, ldmc = dry_weight/fresh_weight, 
         trt = case_when(strip %in% c(1,3,5) ~ "b", strip %in% c(2,4,6) ~ "u"))

unique(sla$plant)





#---- Visualilze ----
ggplot(sla, mapping = aes(x = leaf_area, y = fresh_weight)) + 
  geom_point()




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

SLABvsU <- ggplot(sla, aes(x = trt, y = sla, fill = trt)) + 
  geom_boxplot(alpha = 0.7, color = "black") + 
  geom_point(alpha = 0.7, color = "gray") + # <-- alpha gives a soft, see-through fill
  theme_classic(base_size = 18) +
  labs(
    x = "Burned vs. Unburned",
    y = "Surface Leaf Area (SLA)",
    title = "SLA in burned vs. unburned "
  ) +
  scale_color_manual(values = c("#E91E63", "#4CAF50")) +
  scale_fill_manual(values = c("#E91E63", "#4CAF50")) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#E8F5E9", color = NA),
    strip.text = element_text(face = "bold", color = "#1B5E20", size = 15),
    axis.text.y = element_text(color = "#424242", size = 16),  # smaller font for y-axis
    axis.text.x = element_text(color = "#424242", size = 24),
    axis.title = element_text(face = "bold", size = 20),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 28, face = "bold")
  )
  

ggsave(
  filename = "SLABvsU.jpeg",
  plot = SLABvsU,
  device = "jpeg",
  dpi = 600,
  width = 9,
  height = 5,
  units = "in"
)

ggplot(sla, mapping = aes(x = as.factor(round), y = ldmc, color = trt)) + 
  geom_boxplot() + 
  facet_wrap(~ plant) 

lm1 <- glmmTMB(sla ~ trt * plant + (1|strip), data = sla %>% 
                 filter(round == 2, plant != "dican"))

simulateResiduals(lm1, plot = T)

summary(lm1)
Anova(lm1)

emmeans(lm1, pairwise ~ trt|plant)



