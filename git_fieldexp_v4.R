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

################################.
### first part of experiment ###.
################################.

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

#################################.
### second part of experiment ###.
#################################.

# creating starting densities for each species within trt

cage_exp_2 <- cage_exp_general %>%
  filter(round == 3, type == "surv", alive == 1) %>% 
  group_by(strip, block, cage, burn, trt, dep, high_low, sp) %>%
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
  group_by(strip, block, cage, sp, round) %>%
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

############################# PART 1 OF EXPERIMENT ###########################################

# ---- DATA VIZUALIZATION AND MODELS ----- 

#### H1 - Plants were nutritionally higher in burned vs. unburned ####

### Round 2 of plant collection for SLA & LDMC ###

## Visuals ##

# SLA 

ggplot(slaldmc, aes(x = trt, y = sla, fill = trt)) + 
  geom_boxplot() + 
  facet_wrap(~ plant, scales = "free_y") + 
  scale_fill_viridis_d(option = "magma", begin = 0.5 , end = 0.9) +
  theme_classic(base_size = 22) +
  labs(x = "Burn Treatment",
       y = "Surface Leaf Area (SLA)",
       title = "SLA Across Treatments") 

# LDMC 

ggplot(slaldmc, aes(x = trt, y = ldmc, fill = trt)) + 
  geom_boxplot() + 
  facet_wrap(~ plant, scales = "free_y") + 
  scale_fill_viridis_d(option = "magma", begin = 0.5 , end = 0.9) +
  theme_classic(base_size = 22) +
  labs(x = "Burn Treatment",
       y = "Leaf-Dry Matter-Content (LDMC)",
       title = "LDMC Across Treatments") 

# CHOICE-ASSAY 

# consumed leaf area 
ggplot(pref_trial_r4, aes(x = trt, y = tot_cons, fill = trt)) + 
  geom_boxplot() + 
  scale_fill_viridis_d(option = "magma", begin = 0.5 , end = 0.9) +
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

graphics.off() 
### Models ##

# ACH DD
DD_ach <- glmmTMB(perc_survival ~ grass_perc * high_low * burn,
                              data = cages_1 %>% 
                                filter(round == 2, sp == "ach", dep == "monoculture"),
                  family = "ordbeta")

plot(simulateResiduals(DD_ach))
Anova(DD_ach)
emmeans(DD_ach,pairwise ~ high_low|burn, type = "response")
summary(emtrends(DD_ach, ~ high_low * burn, var = "grass_perc"), infer = TRUE)

# APT DD
DD_apt <- glmmTMB(perc_survival ~ grass_perc * high_low * burn,
                  data = cages_1 %>% 
                    filter(sp == "apt", dep == "monoculture", round == 2,),
                  family = "ordbeta")

plot(simulateResiduals(DD_apt))
summary(DD_apt)
Anova(DD_apt)
emmeans(DD_apt, pairwise ~ high_low|burn, type = "response", infer = T)
summary(emtrends(DD_apt, ~ high_low * burn, var = "grass_perc"), infer = TRUE)

graphics.off()
ggplot(cages_1 %>% 
         filter(dep == "monoculture", round == 2),
       aes(x = grass_perc, y = perc_survival,
           color = factor(high_low, levels = c("low", "high")))) +
  geom_jitter(width = 0.02, height = 0,
              size = 2.5, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(~ sp,
             labeller = as_labeller(c(
               ach = "italic('A. carinatum')",
               apt = "italic('A. sphenarioides')"
             ), label_parsed)) +
  scale_color_viridis_d(option = "magma", end = 0.5) +
  theme_bw(base_size = 20) +
  labs(x = "Percent grass",
       y = "Survival proportion of grasshoppers",
       title = "Grasshopper monoculture survival across grass %",
       color = "Density") +
  theme(plot.title = element_text(hjust = 0.4, face = "bold", size = 18))

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



# FD Apt only 

#graph 

ggplot(cages_1 %>% 
         filter(sp == "apt", trt != "apt_low", trt != "control", round == 2),
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
       title = "Aptenopedes mixture survival in burned vs. unburned",
       color = "Treatment") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22), 
        legend.position = "none")

# model 
FD_apt <- glmmTMB(perc_survival ~ grass_perc * trt * burn,
                              data = cages_1 %>% 
                                filter(round == 2, sp == "apt", trt != "apt_low"), family = "ordbeta")

plot(simulateResiduals(FD_apt))
summary(FD_apt)
Anova(FD_apt)
emmeans(FD_apt,pairwise ~ trt|burn, type = "response")

########## No effects of burn treatment on frequency dependence ###

############################# PART 2 OF EXPERIMENT ###########################################

#### H2 - Density dependence in both grasshopper species would be weaker in burned than unburned plots ####

ggplot(cage_exp_allrounds %>% 
         filter(dep == "monoculture"),
       aes(x = total_start, y = alive_n, color = burn)) +
  geom_point(size = 2.8,
             alpha = 0.7) +
  geom_smooth(method = "lm") +
  geom_line(data = bh_preds,
            aes(x = total_start, y = alive_pred),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 1.2) +
  facet_grid(round ~ sp) +
  scale_color_viridis_d(option = "magma", begin = 0.5, end = 0.75,
                        labels = c("u" = "Unburned",
                                   "b" = "Burned")) +
  theme_bw(base_size = 20) +
  labs(x = "Total density",
       y = "Number alive",
       title = "Survival over time in the second part of the experiment",
       color = "Burn treatment") +
  theme(plot.title = element_text(hjust = 0.4,
                                  face = "bold",
                                  size = 22),
        strip.background = element_rect(fill = "grey95",
                                        color = "black"),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

ggplot(cage_exp_allrounds %>% 
         filter(dep == "monoculture", sp == "ach"),
       aes(x = burn, y = prop_survival, color = high_low)) +
  geom_boxplot() + 
  geom_point(size = 2.8,
             alpha = 0.7) +
  facet_grid(round ~ sp) +
  scale_color_viridis_d(option = "magma", begin = 0.5, end = 0.75,
                        labels = c("u" = "Unburned",
                                   "b" = "Burned")) +
  theme_bw(base_size = 20) +
  labs(x = "Total density",
       y = "% survival",
       title = "Survival over time in the second part of the experiment",
       color = "Burn treatment") +
  theme(plot.title = element_text(hjust = 0.4,
                                  face = "bold",
                                  size = 22),
        strip.background = element_rect(fill = "grey95",
                                        color = "black"),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

ggplot(cage_exp_allrounds %>% 
         filter(dep == "monoculture", sp == "apt"),
       aes(x = burn, y = prop_survival, color = high_low)) +
  geom_boxplot() + 
  geom_point(size = 2.8,
             alpha = 0.7) +
  facet_grid(round ~ sp) +
  scale_color_viridis_d(option = "magma", begin = 0.5, end = 0.75,
                        labels = c("u" = "Unburned",
                                   "b" = "Burned")) +
  theme_bw(base_size = 20) +
  labs(x = "Total density",
       y = "% survival",
       title = "Survival over time in the second part of the experiment",
       color = "Burn treatment") +
  theme(plot.title = element_text(hjust = 0.4,
                                  face = "bold",
                                  size = 22),
        strip.background = element_rect(fill = "grey95",
                                        color = "black"),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

summary(bh_ach_r4)



# ACH GML model 

# using starting density 

ach_2_DD <- glmmTMB(cbind(alive_n, density - alive_n) ~ density * burn * as.factor(round) + (1 | block/cage),
                    family = binomial,
                    data = cage_exp_allrounds %>%
                      filter(dep == "monoculture", sp == "ach"))

simulateResiduals(ach_2_DD, plot = T)
summary(ach_2_DD)
Anova(ach_2_DD)
emtrends(ach_2_DD,pairwise ~ burn|round, var = "density", infer = T)

############ ATKE CARE OF THIS 

unique(cage_exp_allrounds$prop_survival)

apt_2_DD <- glmmTMB(perc_survival ~ high_low * burn * as.factor(round) + (1 | block/cage),
                    family = "ordbeta",
                    data = surv_plant %>%
                      filter(dep == "monoculture", sp == "apt", trt != "control", 
                             round > 3))

simulateResiduals(apt_2_DD, plot = T)
summary(apt_2_DD)
Anova(apt_2_DD)
emmeans(apt_2_DD, pairwise ~ high_low|round:burn, type = "response")
emmeans(apt_2_DD, pairwise ~ burn|round, type = "response")


###########################

# using high_low categories

ach_2_DD_hl <- glmmTMB(cbind(alive_n, density - alive_n) ~ high_low * burn * as.factor(round) + (1 | block/cage),
                    family = binomial,
                    data = cage_exp_allrounds %>%
                      filter(dep == "monoculture", sp == "ach"))

simulateResiduals(ach_2_DD_hl, plot = T)
summary(ach_2_DD_hl)
Anova(ach_2_DD_hl)
emmeans(ach_2_DD_hl,pairwise ~ high_low|round:burn, type = "response") 

# ACH DD GLM model on just round 5 

ach_2_DD_5 <- glmmTMB(cbind(alive_n, density - alive_n) ~ density * burn + (1 | block/cage),
                      family = binomial,
                      data = cage_exp_allrounds %>%
                        filter(dep == "monoculture", sp == "ach", round == 5))

simulateResiduals(ach_2_DD_5, plot = T)
summary(ach_2_DD_5)
Anova(ach_2_DD_5)
emtrends(ach_2_DD_5,pairwise ~ burn, var = "density", infer = T)

# APT DD GLM model on just round 5 

apt_2_DD_5 <- glmmTMB(cbind(alive_n, density - alive_n) ~ density * burn + (1 | block/cage),
                      family = binomial,
                      data = cage_exp_allrounds %>%
                        filter(dep == "monoculture", sp == "apt", round == 5))

simulateResiduals(apt_2_DD_5, plot = T)
summary(apt_2_DD_5)
Anova(apt_2_DD_5)
emtrends(apt_2_DD_5,pairwise ~ burn, var = "density", infer = T)

#### H3 - Frequency dependence would  be weaker in burned vs. unburned plots ####

## Visuals ##

cage_exp_2_fd_plot <- surv_plant %>% 
  filter(round == 5, trt != "ach_low", trt != "apt_low", trt != "control" ) %>% 
  mutate(freq_label = case_when(
    sp == "ach" & trt == "ach_33" ~ "33%",
    sp == "ach" & trt == "ach_66" ~ "66%",
    sp == "ach" & trt == "ach_high" ~ "100%",
    sp == "ach" & trt == "apt_high" ~ "0%",
    sp == "apt" & trt == "ach_33" ~ "66%",
    sp == "apt" & trt == "ach_66" ~ "33%",
    sp == "apt" & trt == "ach_high" ~ "0%",
    sp == "apt" & trt == "apt_high" ~ "100%" ))

ggplot(cage_exp_2_fd_plot, aes(x = factor(freq_label,
                                          levels = c("0%", "33%", "66%", "100%")),
                               y = perc_survival,
                               color = factor(freq_label, levels = c("0%", "33%", "66%", "100%")))) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  facet_grid(burn ~ sp) +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.5) +
  theme_bw(base_size = 20) +
  labs(x = "Focal species frequency", y = "Survival Proportion", title = "Grasshopper mixture survival in burned vs. unburned", color = "Frequency") +
  theme(plot.title = element_text(hjust = 0.4, face = "bold", size = 22), legend.position = "none", strip.background = element_rect(fill = "grey95", color = "black"), strip.text = element_text(face = "bold"))

# ACH FD all rounds , 4/7/26

ach_2_fd_bin <- glmmTMB(prop_survival ~ trt * burn * as.factor(round) + (1|block/cage),
                    data = cage_exp_allrounds %>% 
                      filter(sp == "ach", trt != "ach_low", prop_survival < 1.01, round != 6), family = "binomial")

plot(simulateResiduals(ach_2_fd_bin))
summary(ach_2_fd_bin)
Anova(ach_2_fd_bin)
emmeans(ach_2_fd_bin,pairwise ~ trt|burn:round, type = "response")

# ACH FD all rounds 4/7

ach_2_fd_allrounds <- glmmTMB(perc_survival ~ grass_perc + trt * burn * as.factor(round) + (1|block/cage),
                    data = surv_plant %>% 
                      filter(sp == "ach", trt != "ach_low", round > 3) %>% 
                      filter(round != 6), family = "ordbeta")

plot(simulateResiduals(ach_2_fd_allrounds))
summary(ach_2_fd_allrounds)
Anova(ach_2_fd_allrounds)
emmeans(ach_2_fd_allrounds,pairwise ~ trt|burn:round, type = "response")

# ACH FD round 5 4/7

ach_2_fd <- glmmTMB(perc_survival ~ trt * burn + (1|block),
                    data = surv_plant %>% 
                      filter(sp == "ach", trt != "ach_low", round == 5), family = "ordbeta")

plot(simulateResiduals(ach_2_fd))
summary(ach_2_fd)
Anova(ach_2_fd)
emmeans(ach_2_fd,pairwise ~ trt|burn, type = "response")

# ACH FD round 5 

ach_2_fd <- glmmTMB(perc_survival ~ grass_perc * trt * burn,
                                  data = surv_plant %>% 
                                    filter(sp == "ach", trt != "ach_low", round == 5), family = "ordbeta")

plot(simulateResiduals(ach_2_fd))
summary(ach_2_fd)
Anova(ach_2_fd)
emmeans(ach_2_fd,pairwise ~ trt|burn, type = "response")

# APT OBSERVATIONS

surv_plant_apt <- cage_exp_allrounds %>% 
  filter(sp == "apt") %>% 
  group_by(round, trt, burn) %>% 
  summarise(total = n(), .groups = "drop")

surv_plant_apt <- cage_exp_allrounds %>% 
  filter(sp == "apt") %>% 
  group_by(round, trt, burn, alive_n) %>% 
  summarise(total = n(), .groups = "drop")

surv_plant_apt <- cage_exp_allrounds %>% 
  filter(sp == "apt") %>% 
  group_by(round, trt, burn) %>% 
  summarise(avg = mean(alive_n), .groups = "drop")

#perc_survival = mean(alive),
#days = mean(days),
#dens = n_distinct(ind),
#.groups = "drop"
#)  

# APT FD rounds 4 & 5

apt_2_fd_allrounds <- glmmTMB(perc_survival ~ grass_perc + trt * burn * as.factor(round) + (1|block/cage),
                    data = surv_plant %>% 
                      filter(sp == "apt", trt != "apt_low", round != 6, trt != "control") %>% 
                      filter(round > 3), family = "ordbeta")

plot(simulateResiduals(apt_2_fd_allrounds))

summary(apt_2_fd_allrounds)
Anova(apt_2_fd_allrounds)
emmeans(apt_2_fd_allrounds,pairwise ~ trt|burn:round, type = "response")

# APT FD 

apt_2_fd <- glmmTMB(perc_survival ~ grass_perc * trt * burn,
                                  data = surv_plant %>% 
                                    filter(sp == "apt", trt != "apt_low", round == 5, trt != "control"), family = "ordbeta")

plot(simulateResiduals(apt_2_fd))

summary(apt_2_fd)
Anova(apt_2_fd)
emmeans(apt_2_fd,pairwise ~ trt|burn, type = "response")

########## Effects of burn treatment on frequency dependence and species ###

############################# OVERALL RESULTS ###########################################

anova_table <- function(model, model_name) {
  out <- as.data.frame(Anova(model))
  out$effect <- rownames(out)
  rownames(out) <- NULL
  out$model <- model_name
  out <- out[, c("model", "effect", "Chisq", "Df", "Pr(>Chisq)")]
  out
}

emm_table <- function(emm_obj, model_name) {
  if("emm_list" %in% class(emm_obj)) {
    out <- as.data.frame(summary(emm_obj[[1]]))
  } else {
    out <- as.data.frame(summary(emm_obj))
  }
  out$model <- model_name
  out
}

contrast_table <- function(emm_obj, model_name) {
  if("emm_list" %in% class(emm_obj)) {
    out <- as.data.frame(summary(emm_obj[[2]]))
    out$model <- model_name
    out
  }
}

part1_models <- list(sla = sla, ldmc = ldmc, preference = preference, DD_ach = DD_ach,
                     DD_apt = DD_apt, FD_ach = FD_ach, FD_apt = FD_apt)

part2_models <- list(ach_2_DD = ach_2_DD, ach_2_DD_5 = ach_2_DD_5, apt_2_DD_5 = apt_2_DD_5,
                     ach_2_FD = ach_2_FD, apt_2_fd = apt_2_fd)

all_models <- c(part1_models, part2_models)

part1_emm <- list(sla = emmeans(sla, pairwise ~ trt|plant),
                  ldmc = emmeans(ldmc, pairwise ~ trt|plant),
                  preference = emmeans(preference, ~ trt),
                  DD_ach = emmeans(DD_ach, pairwise ~ high_low|burn, type = "response"),
                  DD_apt = emmeans(DD_apt, pairwise ~ high_low|burn, type = "response"),
                  FD_ach = emmeans(FD_ach, pairwise ~ trt|burn, type = "response"),
                  FD_apt = emmeans(FD_apt, pairwise ~ trt|burn, type = "response"))

part2_emm <- list(ach_2_DD = emtrends(ach_2_DD, pairwise ~ burn|round, var = "density", infer = TRUE),
                  ach_2_DD_5 = emtrends(ach_2_DD_5, pairwise ~ burn, var = "density", infer = TRUE),
                  apt_2_DD_5 = emtrends(apt_2_DD_5, pairwise ~ burn, var = "density", infer = TRUE),
                  ach_2_FD = emmeans(ach_2_FD, pairwise ~ trt|burn, type = "response"),
                  apt_2_fd = emmeans(apt_2_fd, pairwise ~ trt|burn, type = "response"))

all_emm <- c(part1_emm, part2_emm)

anova_results <- bind_rows(anova_table(sla, "sla"),
                           anova_table(ldmc, "ldmc"),
                           anova_table(preference, "preference"),
                           anova_table(DD_ach, "DD_ach"),
                           anova_table(DD_apt, "DD_apt"),
                           anova_table(FD_ach, "FD_ach"),
                           anova_table(FD_apt, "FD_apt"),
                           anova_table(ach_2_DD, "ach_2_DD"),
                           anova_table(ach_2_DD_5, "ach_2_DD_5"),
                           anova_table(apt_2_DD_5, "apt_2_DD_5"),
                           anova_table(ach_2_FD, "ach_2_FD"),
                           anova_table(apt_2_fd, "apt_2_fd"))

emm_results <- bind_rows(emm_table(part1_emm$sla, "sla"),
                         emm_table(part1_emm$ldmc, "ldmc"),
                         emm_table(part1_emm$preference, "preference"),
                         emm_table(part1_emm$DD_ach, "DD_ach"),
                         emm_table(part1_emm$DD_apt, "DD_apt"),
                         emm_table(part1_emm$FD_ach, "FD_ach"),
                         emm_table(part1_emm$FD_apt, "FD_apt"),
                         emm_table(part2_emm$ach_2_DD, "ach_2_DD"),
                         emm_table(part2_emm$ach_2_DD_5, "ach_2_DD_5"),
                         emm_table(part2_emm$apt_2_DD_5, "apt_2_DD_5"),
                         emm_table(part2_emm$ach_2_FD, "ach_2_FD"),
                         emm_table(part2_emm$apt_2_fd, "apt_2_fd"))

contrast_results <- bind_rows(contrast_table(part1_emm$sla, "sla"),
                              contrast_table(part1_emm$ldmc, "ldmc"),
                              contrast_table(part1_emm$preference, "preference"),
                              contrast_table(part1_emm$DD_ach, "DD_ach"),
                              contrast_table(part1_emm$DD_apt, "DD_apt"),
                              contrast_table(part1_emm$FD_ach, "FD_ach"),
                              contrast_table(part1_emm$FD_apt, "FD_apt"),
                              contrast_table(part2_emm$ach_2_DD, "ach_2_DD"),
                              contrast_table(part2_emm$ach_2_DD_5, "ach_2_DD_5"),
                              contrast_table(part2_emm$apt_2_DD_5, "apt_2_DD_5"),
                              contrast_table(part2_emm$ach_2_FD, "ach_2_FD"),
                              contrast_table(part2_emm$apt_2_fd, "apt_2_fd"))

anova_results
emm_results
contrast_results
