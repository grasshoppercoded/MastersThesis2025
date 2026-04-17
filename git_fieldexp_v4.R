#---- SETUP ----

# download packages  

library(tidyverse)
library(emmeans)
library(car)
library(DHARMa)
library(glmmTMB)
library(performance)
library(easystats)
library(betareg)
library(tibble)
library(svglite)
library(patchwork)

rm(list = ls())

# download data 

slaldmc <- read.csv("MastersThesis2025/SLA_LDMC_summ25.csv")
head(slaldmc)
summary(slaldmc)
str(slaldmc)

pref_trial <- read.csv("feeding_trial_2025_final.csv")
head(pref_trial)
str(pref_trial)
summary(pref_trial)

cage_exp <- read.csv("cagestock_summ25_final_spider.csv") ### grasshopper survival data
head(cage_exp)
print(cage_exp)
summary(cage_exp)

plant_abundance <- read.csv("plants_summ25_relativeabundance_fixed.csv") ### plant abundance data 
head(plant_abundance)
summary(plant_abundance)

#---- CLEANING DATA ----

### SLA-LDMC #### 

slaldmc <- slaldmc %>% 
  select(-c(SLA, LDMC, notes)) %>%
  drop_na(fresh_weight) %>% # removes "dican" species, since we did not collect data for this 
  filter(round == 2) %>% # using this round BECAUSE
  mutate(sla = leaf_area/fresh_weight, ldmc = dry_weight/fresh_weight, 
         trt = case_when(strip %in% c(1,3,5) ~ "b", strip %in% c(2,4,6) ~ "u"))

unique(slaldmc$plant)
str(slaldmc)
summary(slaldmc)
hist(slaldmc$leaf_area)

### CHOICE-ASSAY ####

pref_trial_r4 <- pref_trial %>% 
  filter(round == 4) %>% 
  mutate(frass = poo_weight - vial_weight) %>% 
  group_by(mesocosm, trt) %>% 
  summarize(tot_cons = sum(cons_area, na.rm = TRUE),
            sex = first(sex),
            body_lg = first(body_lg),
            frass = sum(frass, na.rm = TRUE),
            .groups = "drop") %>% 
  group_by(mesocosm) %>% 
  filter(!all(tot_cons == 0)) %>% 
  ungroup()

### CAGE-EXPERIMENT SURVIVAL DATA ####

##################################.
#### First part of experiment ####
##################################.

# create days since first stocking 

cage_exp$date <- as.Date(cage_exp$date, format = "%m/%d/%Y")
stock_date <- as.Date("6/20/25", format = "%m/%d/%Y")
cage_exp$days <- as.numeric(cage_exp$date - stock_date)
head(cage_exp)

# adding burned vs. unburned, mixtures and monocultures,
# proportion survival, and simplifying
# contains all rounds
# general dataset 

cage_exp_general <- cage_exp %>% 
  select(c(strip, block, cage, round, type, alive, trt, sp, ind, days)) %>% 
  mutate(burn = case_when(strip %in% c(1,3,5) ~ "b", 
                          strip %in% c(2,4,6) ~ "u")) %>%
  mutate(dep = case_when(trt %in% c("ach_high","apt_high","apt_low", "ach_low") ~ "monoculture", 
                         trt %in% c("ach_66","ach_33") ~ "mixture")) %>% 
  mutate(high_low = case_when(trt %in% c("ach_high", "apt_high") ~ "high", 
                              trt %in% c("ach_low", "apt_low") ~ "low")) %>% 
  filter(alive %in% c(0, 1) | is.na(alive)) %>%
  mutate(alive = as.numeric(alive)) %>% 
  drop_na(strip) 

cage_exp_surv_1 <- cage_exp_general%>% 
  group_by(strip, block, cage, round, type, trt, sp, burn, dep, high_low) %>% 
  summarise(perc_survival = mean(alive),
            days = mean(days),
            dens = n_distinct(ind),
            .groups = "drop")

##################################.
#### Second part of experiment ####
##################################.

# creating starting densities for each species within trt

cage_exp_2 <- cage_exp_general %>%
  filter(round == 3, type == "surv", alive == 1) %>% 
  group_by(strip, block, cage, round, type, trt, sp, burn, dep, high_low) %>%
  summarize(density = n_distinct(ind), .groups = "drop") %>%
  mutate(density = case_when(
    cage == 23 & sp == "ach" ~ 3, # manually correct starting densities for cages where round 3 counts were underestimated,
    cage == 35 & sp == "apt" ~ 4, # which otherwise caused survival proportions to exceed 1 in later rounds
    cage == 44 & sp == "ach" ~ 5,
    TRUE ~ density))

# creating total density for each cage  

part_2_wide <- cage_exp_2 %>%
  pivot_wider(names_from = sp,
              values_from = density,
              values_fill = 0) %>%
  mutate(total_start = ach + apt,
         freq_ach = (ach / total_start),
         freq_apt = (apt / total_start))

# separating frequency and density dependence & adding survival across later rounds

cage_exp_surv_2 <- cage_exp_general %>%
  filter(round %in% c(4, 5, 6), alive == 1) %>%
  group_by(strip, block, cage, round, type, trt, sp, burn, dep, high_low) %>%
  summarize(alive_n = n_distinct(ind), .groups = "drop") %>%
  complete(strip, block, cage, sp, round = c(4, 5, 6), fill = list(alive_n = 0)) %>%
  left_join(cage_exp_2 %>%
              select(strip, block, cage, sp, density),
            by = c("strip", "block", "cage", "sp")) %>%
  left_join(part_2_wide %>%
              select(strip, block, cage, total_start, freq_ach, freq_apt),
            by = c("strip", "block", "cage")) %>%
  mutate(perc_survival = alive_n / density) %>%
  filter(density > 0)

###################################.
##### Grass Composition Data ######.
###################################.

# adding plant relative abundance data 

plant_summary <- plant_abundance %>% 
  mutate(plant = str_trim(plant),
         burn = b_u,
         veg = case_when(plant %in% c("wide", "cent", "dican", "bb") ~ "grass",
                         plant %in% c("silky", "dog", "other") ~ "forb")) %>% 
  select(-c(notes, bare, b_u)) %>% 
  filter(plant != "", round == 1) %>% 
  group_by(cage, burn, round, veg) %>% 
  summarise(total = sum(perc), .groups = "drop") %>% 
  pivot_wider(names_from = veg, values_from = total, values_fill = 0) %>% 
  mutate(grass_forb_ratio = grass / forb,
         grass_perc = (grass / (grass + forb)) * 100)

# final dataset for first part of exp 

cages_1 <- cage_exp_surv_1 %>% 
  left_join(plant_summary %>% select(cage, grass_perc), by = "cage")

# final dataset for second part of exp 

cages_2 <- cage_exp_surv_2 %>%
  left_join(plant_summary %>% select(cage, grass_perc), by = "cage") %>%
  mutate(perc_survival = alive_n / density) %>%
  filter(density > 0)

# ---- DATA VIZUALIZATION AND MODELS ----- 


## ---- PHASE 1 OF EXPERIMENT ---- 

#### H1 - Plants were nutritionally higher in burned vs. unburned ####

### Round 2 of plant collection for SLA & LDMC ###

## Visuals ##

# SLA 

ggplot(slaldmc, aes(x = trt, y = sla, fill = trt)) + 
  geom_boxplot() + 
  facet_wrap(~ plant, scales = "free_y") + 
  scale_fill_viridis_d(option = "magma", begin = 0.5, end = 0.75) +
  theme_classic(base_size = 22) +
  labs(x = "Burn Treatment",
       y = "Surface Leaf Area (SLA)",
       title = "SLA Across Treatments") 

# LDMC 

ggplot(slaldmc, aes(x = trt, y = ldmc, fill = trt)) + 
  geom_boxplot() + 
  facet_wrap(~ plant, scales = "free_y") + 
  scale_fill_viridis_d(option = "magma", begin = 0.5 , end = 0.75) +
  theme_classic(base_size = 22) +
  labs(x = "Burn Treatment",
       y = "Leaf-Dry Matter-Content (LDMC)",
       title = "LDMC Across Treatments") 

# CHOICE-ASSAY 

# consumed leaf area 
ggplot(pref_trial_r4, aes(x = trt, y = tot_cons, fill = trt)) + 
  geom_boxplot() + 
  scale_fill_viridis_d(option = "magma", begin = 0.5 , end = 0.75) +
  theme_classic(base_size = 22) +
  labs(y = "Total Herbivory", x = "Burned vs. Unburned", 
       title = "Herbivory in Burned vs. Unburned")

## Models ##  

# SLA 

sla <- glmmTMB(sla ~ trt * plant + (1|strip), data = slaldmc)

simulateResiduals(sla, plot = T)
summary(sla)
Anova(sla)
emmeans(sla, pairwise ~ trt|plant)

# LDMC 

ldmc <- glmmTMB(ldmc ~ trt * plant + (1|strip), data = slaldmc)

simulateResiduals(ldmc, plot = T)
summary(ldmc)
Anova(ldmc)
emmeans(ldmc, pairwise ~ trt|plant)

# CHOICE-ASSAY 

preference <- glmmTMB(tot_cons ~ trt + body_lg + (1|mesocosm), data = pref_trial_r4)

simulateResiduals(preference, plot=T)

summary(preference)
Anova(preference)
emmeans(preference, pairwise ~ trt)

########## No effects of burn treatment on SLA, LDMC, or feeding preferences ###

#### H2 - Density dependence in both grasshopper species would be weaker in burned than unburned plots ####

### Visuals ##

# Graph 

ggplot(cages_1 %>% 
         filter(dep == "monoculture", round == 2),
       aes(x = factor(high_low, levels = c("low", "high")), y = perc_survival)) +
  geom_jitter(aes(color = factor(high_low, levels = c("low", "high"))),
              width = 0.12, height = 0,
              size = 2.5, alpha = 0.5) +
  stat_summary(aes(group = 1),
               fun = mean,
               geom = "line",
               linewidth = 1.2,
               color = "black") +
  stat_summary(aes(color = factor(high_low, levels = c("low", "high"))),
               fun = mean,
               geom = "point",
               size = 4) +
  facet_grid(sp ~ burn) +
  scale_color_viridis_d(option = "magma", end = 0.5) +
  theme_bw(base_size = 20) +
  labs(x = "Density",
       y = "Survival Proportion",
       title = "Grasshopper monoculture survival in burned vs. unburned",
       color = "Density") +
  theme(plot.title = element_text(hjust = 0.4, face = "bold", size = 18))

### Models ##

# ACH DD
DD_ach <- glmmTMB(perc_survival ~ grass_perc * high_low * burn,
                              data = cages_1 %>% 
                                filter(round == 2, sp == "ach", dep == "monoculture"),
                  family = "ordbeta")

plot(simulateResiduals(DD_ach))
Anova(DD_ach)
summary(DD_ach)
emmeans(DD_ach,pairwise ~ high_low|burn, type = "response")
emtrends(DD_ach, ~ 1, var = "grass_perc", infer = TRUE)

# APT DD
DD_apt <- glmmTMB(perc_survival ~ grass_perc * high_low * burn,
                  data = cages_1 %>% 
                    filter(sp == "apt", dep == "monoculture", round == 2,),
                  family = "ordbeta")

plot(simulateResiduals(DD_apt))
summary(DD_apt)
Anova(DD_apt)
emmeans(DD_apt, pairwise ~ high_low|burn, type = "response", infer = T)
emtrends(DD_apt, ~ 1, var = "grass_perc", infer = TRUE)

# final graph 

prop_surv_grass <- ggplot(cages_1 %>% 
         filter(dep == "monoculture", round == 2),
       aes(x = grass_perc, y = perc_survival,
           color = sp)) +
  geom_jitter(width = 0.02, height = 0,
              size = 2.5, alpha = 0.5) +
  geom_smooth(method = "lm", se = T, alpha = 0.25) +
  facet_grid(~ sp,
             labeller = as_labeller(c(
               ach = "italic('A. carinatum')",
               apt = "italic('A. sphenarioides')"
             ), label_parsed)) +
  scale_color_viridis_d(option = "magma", end = 0.5, labels = c("A. carinatum", "A. sphenarioides")) +
  theme_bw(base_size = 23) +
  labs(x = "Percent grass",
       y = "Survival proportion of grasshoppers",
       color = "Species") +
  theme(plot.title = element_text(hjust = 0.4, face = "bold", size = 18))

ggsave("prop_surv_grass.svg", plot = prop_surv_grass, width = 10, height = 6)

########## No effects of burn treatment on ACH or APT density dependence ###

#### H3 - Frequency dependence would  be weaker in burned vs. unburned plots ####

## Visuals ##

# FD Achurum only 

ggplot(cages_1 %>% 
         filter(sp == "ach", trt != "ach_low", round == 2,),
       aes(x = trt, y = perc_survival, color = trt)) +
  geom_jitter(width = 0.12,
              size = 2.5,
              alpha = 0.5) +
  stat_summary(fun = mean,
               geom = "point",
               size = 4) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.2,
               linewidth = 1) +
  facet_grid(~ burn) + 
  scale_color_viridis_d(option = "magma", end = 0.5) +
  scale_x_discrete(labels = c("ach_33" = "33%",
                              "ach_66" = "66%",
                              "ach_high" = "100%")) +
  theme_bw(base_size = 20) + 
  labs(x = "Frequency Treatment", 
       y = "Survival Proportion", 
       title = "Achurum mixture survival in burned vs. unburned",
       color = "Treatment") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22), 
        legend.position = "none")

#model
FD_ach <- glmmTMB(perc_survival ~ grass_perc * trt * burn,
                                  data = cages_1 %>% 
                                    filter(round == 2, sp == "ach", trt != "ach_low"), family = "ordbeta")

plot(simulateResiduals(FD_ach))
summary(FD_ach)
Anova(FD_ach)
emmeans(FD_ach, pairwise ~ trt|burn, type = "response", at = list(grass_perc = 30))
emmeans(FD_ach, pairwise ~ trt|burn, type = "response", at = list(grass_perc = 60))
emmeans(FD_ach, pairwise ~ trt|burn, type = "response", at = list(grass_perc = 90))

grass_perc_emmeans_ach_fd <- as_tibble(emmeans(FD_ach, ~ trt:burn|grass_perc, type = "response", at = list(grass_perc = c(30, 60, 90))))

grass_perc_ach_fd <- ggplot(grass_perc_emmeans_ach_fd,
                            aes(x = burn, y = response, color = trt)) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE),
                position = position_dodge(width = 0.45),
                width = 0.4,
                linewidth = 0.8) +
  geom_point(position = position_dodge(width = 0.45), size = 4) +
  facet_grid(~ grass_perc) + 
  scale_color_viridis_d(option = "magma", end = 0, begin = 0.5) +
  scale_x_discrete(labels = c("33" = "33%","66" = "66%", "100" = "100%")) +
  theme_bw(base_size = 20) + 
  labs(x = "Burn Treatment",
       y = "Proportion survival",
       color = "Treatment") + 
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"))

ggsave("grass_perc_ach_fd.svg", plot = grass_perc_ach_fd, width = 10, height = 6)

### three-way effect ###

# FD Apt only 

#graph 

ggplot(cages_1 %>% 
         filter(sp == "apt", trt != "apt_low", trt != "control", round == 2),
       aes(x = trt, y = perc_survival, color = trt)) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  facet_grid(~ burn) + 
  scale_color_viridis_d(option = "magma", end = 0.5) +
  scale_x_discrete(labels = c("ach_33" = "33%", "ach_66" = "66%","ach_high" = "100%")) +
  theme_bw(base_size = 20) + 
  labs(x = "Frequency Treatment", 
       y = "Survival Proportion", 
       title = "Aptenopedes mixture survival in burned vs. unburned",
       color = "Treatment") + 
  theme(plot.title = element_text(hjust = 0.4, face = "bold",size = 22), 
        legend.position = "none")

# model 
FD_apt <- glmmTMB(perc_survival ~ grass_perc * trt * burn,
                              data = cages_1 %>% 
                                filter(round == 2, sp == "apt", trt != "apt_low"), family = "ordbeta")

plot(simulateResiduals(FD_apt))
summary(FD_apt)
Anova(FD_apt)
emmeans(FD_apt,pairwise ~ trt|burn, type = "response")
emtrends(FD_apt, ~ 1, var = "grass_perc", infer = TRUE)


########## No effects of burn treatment on apt frequency dependence ###

##--- PHASE 2 OF EXPERIMENT ###########################################


#### H2 - Density dependence in both grasshopper species would be weaker in burned than unburned plots ####

# DD in both species 

library(ggplot2)
library(dplyr)
library(patchwork)

p_ach <- ggplot(cages_2 %>% 
                  filter(dep == "monoculture", sp == "ach"),
                aes(x = total_start, y = perc_survival, color = burn)) +
  geom_point(size = 2.8,
             alpha = 0.7) +
  geom_smooth(method = "lm") +
  facet_wrap(~ round, ncol = 1) +
  scale_color_viridis_d(option = "magma", begin = 0.5, end = 0.75,
                        labels = c("u" = "Unburned", "b" = "Burned")) +
  theme_bw(base_size = 20) +
  labs(x = "Total density",
       y = "Survival Proportion",
       color = "Burn treatment",
       title = "A. carinatum") +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 22),
        strip.background = element_rect(fill = "grey95",
                                        color = "black"),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

p_apt <- ggplot(cages_2 %>% 
                  filter(dep == "monoculture", sp == "apt"),
                aes(x = total_start, y = perc_survival)) +
  geom_point(size = 2.8,
             alpha = 0.7) +
  geom_smooth(method = "lm", color = "black") +
  theme_bw(base_size = 20) +
  scale_color_viridis_d(option = "magma", begin = 0.5, end = 0.75) +
  labs(x = "Total density",
       y = "Survival Proportion",
       title = "A. sphenarioides") +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 22),
        legend.position = "none")

grass_perc_ach_apt_dd <- p_ach + p_apt + plot_layout(widths = c(1, 1))

ggsave("grass_perc_ach_apt_dd.svg", plot = grass_perc_ach_apt_dd, width = 7, height = 8)


# ACH DD model 

ach_2_DD <- glmmTMB(perc_survival ~ density * burn * as.factor(round) + (1 | block/cage),
                    family = "ordbeta",
                    data = cages_2 %>%
                      filter(dep == "monoculture", sp == "ach", trt != "control", 
                             round > 3))

simulateResiduals(ach_2_DD, plot = T)
summary(ach_2_DD)
Anova(ach_2_DD)
emtrends(ach_2_DD, pairwise ~ round:burn, var = "density", infer = T)
emmeans(ach_2_DD, pairwise ~ burn:round, type = "response")

# APT DD model 

apt_2_DD <- glmmTMB(perc_survival ~ density * burn * as.factor(round) + (1 | block/cage),
                    family = "ordbeta",
                    data = cages_2 %>%
                      filter(dep == "monoculture", sp == "apt", trt != "control", 
                             round > 3))

simulateResiduals(apt_2_DD, plot = T)
summary(apt_2_DD)
Anova(apt_2_DD)

emtrends(apt_2_DD, ~ 1, var = "density", infer = T)
densities_emmeans_apt_dd <- as_tibble(emmeans(apt_2_DD, ~ round:burn|density,
                                               type = "response",
                                               at = list(density = c(4, 6, 8)), infer = T))

emmeans(apt_2_DD, pairwise ~ density, type = "response", at = list(density = c(4, 6, 8)), infer = T)

#### H3 - Frequency dependence would  be weaker in burned vs. unburned plots ####

## Visuals ##

cage_exp_2_fd_plot <- cages_2 %>% 
  filter(trt != "ach_low", trt != "apt_low", trt != "control" ) %>% 
  mutate(freq_label = case_when(
    sp == "ach" & trt == "ach_33" ~ "33%",
    sp == "ach" & trt == "ach_66" ~ "66%",
    sp == "ach" & trt == "ach_high" ~ "100%",
    sp == "ach" & trt == "apt_high" ~ "0%",
    sp == "apt" & trt == "ach_33" ~ "66%",
    sp == "apt" & trt == "ach_66" ~ "33%",
    sp == "apt" & trt == "ach_high" ~ "0%",
    sp == "apt" & trt == "apt_high" ~ "100%" ))

# only a. carinatum with burn, round, and treatment 

ach_2_fd_graph <- ggplot(cage_exp_2_fd_plot %>% 
         filter(sp == "ach"), aes(x = factor(freq_label,
                                          levels = c("0%", "33%", "66%", "100%")),
                               y = perc_survival,
                               color = factor(freq_label, levels = c("0%", "33%", "66%", "100%")))) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 1) +
  facet_grid(round ~ burn) +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.5) +
  theme_bw(base_size = 20) +
  labs(x = "Focal species frequency", y = "Survival Proportion", color = "Frequency") +
  theme(plot.title = element_text(hjust = 0.4, face = "bold", size = 22), legend.position = "none", strip.background = element_rect(fill = "grey95", color = "black"), strip.text = element_text(face = "bold"))

ggsave("ach_2_fd_graph.svg", plot = ach_2_fd_graph, width = 8, height = 7)

# only a. sphenarioides with burn and treatment 

apt_2_fd_graph <- ggplot(cage_exp_2_fd_plot %>% 
                           filter(sp == "apt"), aes(x = factor(freq_label,
                                                               levels = c("0%", "33%", "66%", "100%")),
                                                    y = perc_survival,
                                                    color = factor(freq_label, levels = c("0%", "33%", "66%", "100%")))) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 1) +
  facet_grid(~ burn) +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.5) +
  theme_bw(base_size = 20) +
  labs(x = "Focal species frequency", y = "Survival Proportion", color = "Frequency") +
  theme(plot.title = element_text(hjust = 0.4, face = "bold", size = 22), legend.position = "none", strip.background = element_rect(fill = "grey95", color = "black"), strip.text = element_text(face = "bold"))

ggsave("apt_2_fd_graph.svg", plot = apt_2_fd_graph, width = 6, height = 8)


ggplot(cage_exp_2_fd_plot, aes(x = factor(freq_label,
                                          levels = c("0%", "33%", "66%", "100%")),
                               y = perc_survival,
                               color = factor(freq_label, levels = c("0%", "33%", "66%", "100%")))) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 1) +
  facet_grid(burn ~ sp) +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.5) +
  theme_bw(base_size = 20) +
  labs(x = "Focal species frequency", y = "Survival Proportion", color = "Frequency") +
  theme(plot.title = element_text(hjust = 0.4, face = "bold", size = 22), legend.position = "none", strip.background = element_rect(fill = "grey95", color = "black"), strip.text = element_text(face = "bold"))

# ACH FD 

ach_2_fd <- glmmTMB(perc_survival ~ trt * burn * as.factor(round) + (1|block/cage),
                    data = cages_2 %>% 
                      filter(sp == "ach", trt != "ach_low"), family = "ordbeta")

plot(simulateResiduals(ach_2_fd))
summary(ach_2_fd)
Anova(ach_2_fd)
emmeans(ach_2_fd,pairwise ~ trt|burn:round, type = "response")

# APT FD

apt_2_fd <- glmmTMB(perc_survival ~ trt * burn * as.factor(round) + (1|block/cage),
                    data = cages_2 %>% 
                      filter(sp == "apt", trt != "apt_low", trt != "control") %>% 
                      filter(round > 3), family = "ordbeta")

plot(simulateResiduals(apt_2_fd))
summary(apt_2_fd)
Anova(apt_2_fd)
emmeans(apt_2_fd, pairwise ~ trt|burn, type = "response")
emmeans(apt_2_fd, pairwise ~ burn:trt, type = "response")


########## Effects of burn treatment on frequency dependence and species ###
