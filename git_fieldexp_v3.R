#---- SETUP ----

rm(list = ls())

# download packages  

library(tidyverse)
library(emmeans)
library(car)
library(DHARMa)
library(glmmTMB)
library(performance)
library(easystats)
library(betareg)

# download data 

slaldmc <- read.csv("MastersThesis2025/SLA_LDMC_summ25.csv")
head(slaldmc)
summary(slaldmc)
str(slaldmc)

pref_trial <- read.csv("feeding_trial_2025.csv")
head(pref_trial)
str(pref_trial)
summary(pref_trial)

cage_exp <- read.csv("cagestock_summ25_final.csv") ### grasshopper survival data
head(cage_exp)
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

unique(pref_trial$tot_cons)

pref_trial <- pref_trial %>% 
  group_by(mesocosm, trt) %>% 
  summarize(tot_cons = sum(cons_area, na.rm=T)) %>% 
  filter(tot_cons > 0) # removing cases with no herbivory 


### CAGE-EXPERIMENT SURVIVAL DATA ####

# create days since first stocking 

cage_exp$date <- as.Date(cage_exp$date, format = "%m/%d/%Y")
stock_date <- as.Date("6/20/2025", format = "%m/%d/%Y")
cage_exp$days <- as.numeric(cage_exp$date - stock_date)
head(cage_exp)

# adding burned vs. unburned, mixtures and monocultures, proportion survival, and simplifying

cage_exp <- cage_exp %>% 
  mutate(burn = case_when((strip %in% c(1,3,5) ~ "b"), 
                          (strip %in% c(2,4,6) ~ "u"))) %>%
  mutate(dep = case_when((trt %in% c("ach_high","apt_high","apt_low", "ach_low") ~ "monoculture"), 
                         (trt %in% c("ach_66","ach_33") ~ "mixture"))) %>% 
  mutate(high_low = case_when((trt %in% c("ach_high", "apt_high") ~ "high"), 
                              (trt %in% c("ach_low", "apt_low") ~ "low"))) %>% 
  filter(alive %in% c(0, 1) | is.na(alive)) %>%
  mutate(alive = as.numeric(alive)) %>% 
  select(c(strip, block, cage, round, dep, trt, high_low, burn, sp, ind, alive, days)) %>% 
  drop_na(strip)

### PLANT ABUNDANCE DATA ####

# adding plant relative abundance data 

plant_abundance <- plant_abundance %>% 
  mutate(plant = str_trim(plant)) %>% 
  mutate(burn = b_u) %>% 
  mutate(
    veg = case_when(
      plant %in% c("wide", "cent", "dican", "bb") ~ "grass",
      plant %in% c("silky", "dog", "other") ~ "forb")) %>% 
  select(-c(notes, bare, b_u)) %>% 
  filter(plant != "", round == 1) %>% 
  group_by(cage, burn, veg, round) %>% 
  summarise(total = sum(perc)) %>% 
  pivot_wider(names_from = veg, values_from = total) %>% 
  mutate(grass_forb_ratio = grass / forb) %>% 
  mutate(grass_perc = (grass/(grass + forb))*100)

# make new dataframe of joined survival + relative abundance data 

surv_plant <- plant_abundance %>% 
  left_join(cage_exp, by = "cage") ###### x and y in round and strip?????

head(surv_plant)
print(surv_plant)
str(surv_plant)


# ---- DATA VIZUALIZATION AND MODELS ----- 

#### H1 - Plants were nutritionally higher in burned vs. unburned ####

### Round 2 of plant collection for SLA & LDMC ###

## Visuals ##

# SLA 

ggplot(slaldmc, aes(x = leaf_area, y = fresh_weight)) + 
  geom_point() #regression relationship 

ggplot(slaldmc, aes(x = trt, y = sla, fill = trt)) + 
  geom_boxplot() + 
  facet_wrap(~ plant, scales = "free_y") + 
  theme_classic(base_size = 20) +
  labs(x = "Burn Treatment",
       y = "Surface Leaf Area (SLA)",
       title = "SLA Across Treatments") 

# LDMC 

ggplot(slaldmc, aes(x = trt, y = ldmc, fill = trt)) + 
  geom_boxplot() + 
  facet_wrap(~ plant, scales = "free_y") + 
  theme_classic(base_size = 20) +
  labs(x = "Burn Treatment",
       y = "Leaf-Dry Matter-Content (LDMC)",
       title = "LDMC Across Treatments") 

# CHOICE-ASSAY 

ggplot(pref_trial, aes(x = trt, y = tot_cons, fill = trt)) + 
  geom_boxplot() + 
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

preference <- glmmTMB(tot_cons ~ trt + (1|mesocosm), data = pref_trial)

simulateResiduals(preference, plot=T)
ggplot(pref_trial, aes(x = tot_cons)) + 
  geom_histogram()
check_overdispersion(preference)
#plot(resid(pref_trial$tot_cons), fitted(pref_trial$tot_cons))

preference_log <- glmmTMB(log(tot_cons) ~ trt + (1|mesocosm), data = pref_trial)

simulateResiduals(preference_log, plot=T)
ggplot(pref_trial, aes(x = log(tot_cons))) + 
  geom_histogram()
check_overdispersion(preference_log)

summary(preference)
Anova(preference)
emmeans(preference, ~ trt)

########## No effects of burn treatment on SLA, LDMC, or feeding preferences ###

#### H2 - Density dependence in both grasshopper species would be weaker in burned than unburned plots ####

## Visuals ##

# Adjusting data set for DD analysis 

cage_exp_dd<- cage_exp %>% 
  group_by(strip, block, cage, trt, dep, sp, round, burn, high_low) %>% 
  summarize(perc = mean(alive),
            .groups = "drop") 

# Graph 

ggplot(cage_exp_dd %>% 
         filter(dep == "monoculture", round == 2), # round 2 is the first survey
       aes(x = factor(high_low, levels = c("low", "high")), y = perc)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(sp ~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Density", 
       y = "Survival Proportion", 
       title = "Grasshoppers monoculture survival in burned vs. unburned") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

## Models ##

DD <- glmmTMB(perc ~ burn * high_low * sp + (1|block), data = cage_exp_dd %>% 
                    filter(dep == "monoculture"), family = "ordbeta")

plot(simulateResiduals(DD))

summary(DD)
Anova(DD)
emmeans(DD,pairwise ~ high_low|burn|sp, type = "response")

########## No effects of burn treatment on density dependence ###

#### H3 - Frequency dependence would  be weaker in burned vs. unburned plots ####

## Visuals ##

# Adjusting data set for FD analysis 

cage_exp_freq <- cage_exp %>% 
  filter(round == 2, dep == "mixture") %>% # round 2 is the first survey
  group_by(strip, block, cage, trt, dep, sp, round, burn, high_low) %>% 
  summarize(perc = mean(alive),
            days = mean(days), 
            dens = n_distinct(ind),
            .groups = "drop") 

ggplot(cage_exp_freq, aes(x = trt, y = perc)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(sp ~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Frequency Treatment", 
       y = "Survival Proportion", 
       title = "Grasshoppers mixture survival in burned vs. unburned") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

## Model ##

FD <- glmmTMB(perc ~ burn * sp * trt + (1|block), data = cage_exp_freq, family = "ordbeta")

plot(simulateResiduals(FD))

summary(FD)
Anova(FD)
emmeans(FD,pairwise ~ trt|burn|sp, type = "response")

########## No effects of burn treatment on frequency dependence ###

