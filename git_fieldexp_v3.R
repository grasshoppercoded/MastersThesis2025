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

cage_exp <- read.csv("cagestock_summ25_final_spider.csv") ### grasshopper survival data
head(cage_exp)
summary(cage_exp)

plant_abundance <- read.csv("plants_summ25_relativeabundance_fixed.csv") ### plant abundance data 
head(plant_abundance)
summary(plant_abundance)

#---- CLEANING DATA ----

############################# PART 1 OF EXPERIMENT ###########################################
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
  filter(round == 4) %>% 
  mutate(frass = poo_weight - vial_weight) %>% 
  group_by(mesocosm, trt, sex) %>% 
  summarize(tot_cons = sum(cons_area, na.rm=T))
  #filter(mesocosm < 18)# removing cases with no herbivory 

pref_trial <- pref_trial %>% 
  filter(round == 4) %>% 
  mutate(frass = poo_weight - vial_weight) %>% 
  group_by(mesocosm, trt, sex) %>% 
  summarize(tot_cons = sum(cons_area, na.rm = TRUE), .groups = "drop") %>% 
  group_by(mesocosm) %>% 
  filter(!all(tot_cons == 0)) %>% 
  ungroup()

### CAGE-EXPERIMENT SURVIVAL DATA ####

# create days since first stocking 

cage_exp$date <- as.Date(cage_exp$date, format = "%m/%d/%Y")
stock_date <- as.Date("6/20/2025", format = "%m/%d/%Y")
cage_exp$days <- as.numeric(cage_exp$date - stock_date)
head(cage_exp)

# adding burned vs. unburned, mixtures and monocultures, proportion survival, and simplifying

cage_exp_surv <- cage_exp %>% 
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

plant_summary <- plant_abundance %>% 
  mutate(plant = str_trim(plant),
         burn = b_u,
         veg = case_when(
           plant %in% c("wide", "cent", "dican", "bb") ~ "grass",
           plant %in% c("silky", "dog", "other") ~ "forb"
         )) %>% 
  select(-c(notes, bare, b_u)) %>% 
  filter(plant != "", round == 1) %>% 
  group_by(cage, burn, round, veg) %>% 
  summarise(total = sum(perc), .groups = "drop") %>% 
  pivot_wider(names_from = veg, values_from = total, values_fill = 0) %>% 
  mutate(
    grass_forb_ratio = grass / forb,
    grass_perc = (grass / (grass + forb)) * 100
  )

# make new dataframe of joined survival + relative abundance data 

surv_summary <- cage_exp_surv %>% 
  filter(round == 2) %>% 
  group_by(strip, block, cage, trt, dep, sp, round, burn, high_low) %>% 
  summarise(
    perc_survival = mean(alive),
    days = mean(days),
    dens = n_distinct(ind),
    .groups = "drop"
  )  

surv_plant <- surv_summary %>% 
  left_join(
    plant_summary %>% select(cage, grass_perc, grass_forb_ratio),
    by = "cage"
  )

###### x and y in round and strip?????

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

preference <- glmmTMB(tot_cons ~ trt * sex + (1|mesocosm), data = pref_trial)

simulateResiduals(preference, plot=T)

summary(preference)
Anova(preference)
emmeans(preference, ~ trt)


########## No effects of burn treatment on SLA, LDMC, or feeding preferences ###

#### H2 - Density dependence in both grasshopper species would be weaker in burned than unburned plots ####

## Visuals ##

# Adjusting data set for DD analysis 

cage_exp_dd<- cage_exp_surv %>% 
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

#### H3a - Frequency dependence would  be weaker in burned vs. unburned plots ####

## Visuals ##

# Adjusting data set for FD analysis 

cage_exp_fd <- cage_exp_surv %>% 
  filter(round == 2, dep == "mixture") %>% # round 2 is the first survey
  group_by(strip, block, cage, trt, dep, sp, round, burn, high_low) %>% 
  summarize(perc = mean(alive),
            days = mean(days), 
            dens = n_distinct(ind),
            .groups = "drop") 

ggplot(cage_exp_fd, aes(x = trt, y = perc)) +
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

FD <- glmmTMB(perc ~ burn * sp * trt + (1|block), data = cage_exp_fd, family = "ordbeta")

plot(simulateResiduals(FD))

summary(FD)
Anova(FD)
emmeans(FD,pairwise ~ trt|burn|sp, type = "response")

########## No effects of burn treatment on frequency dependence ###

#### H3b - Grass composition and abundance may favor the grass-specialist ####

## Visuals ##

ggplot(surv_plant, aes(x = grass_perc, y = perc_survival)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_grid(sp ~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Grass Percentage", 
       y = "Survival Proportion", 
       title = "Survival across grass abundance") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

ggplot(surv_plant, aes(x = grass_perc, y = perc_survival)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_grid(sp ~ trt) + 
  theme_bw(base_size = 20) + 
  labs(x = "Grass Percentage", 
       y = "Survival Proportion", 
       title = "Survival across grass abundance") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

## Model ##

grassdominance <- glmmTMB(perc_survival ~ grass_perc * sp * burn + (1|block), data = surv_plant, family = "ordbeta")

plot(simulateResiduals(grassdominance))

summary(grassdominance)
Anova(grassdominance)
emmeans(grassdominance,pairwise ~ burn|sp, type = "response")

############### HOW ABOUT ONLY ACHARUM?

ggplot(surv_plant %>% 
         filter(sp == "ach"),
       aes(x = grass_perc, y = perc_survival)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_grid(trt ~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Frequency Treatment", 
       y = "Survival Proportion", 
       title = "Grasshoppers mixture survival in burned vs. unburned") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

## Model ##

grassdominance_ach <- glmmTMB(perc_survival ~ grass_perc * trt * burn + (1|block),
                              data = surv_plant %>% 
                                filter(sp == "ach"), family = "ordbeta")

plot(simulateResiduals(grassdominance_ach))

summary(grassdominance_ach)
Anova(grassdominance_ach)
emmeans(grassdominance_ach,pairwise ~ burn|trt, type = "response")






#### ADDITIONAL: SPIDER PREDATION DATA ####

cage_exp_spider <- cage_exp %>% 
  mutate(spider_present = if_else(!is.na(spider) & spider == 1, 1, 0)) %>% 
  group_by(cage) %>% 
  summarise(
    spider_present = max(spider_present, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    spider_present = factor(spider_present, levels = c(0, 1), labels = c("no", "yes"))
  ) %>% 
  filter(cage <= 60) %>% 
  select(c(spider_present, cage))

cage_exp_spider <- cage_exp_surv %>% 
  left_join(exp_cage_spider, by = "cage") %>% 
  filter(round == 2) %>% 
  group_by(strip, block, cage, trt, dep, sp, round, burn, high_low, spider_present) %>% 
  summarize(perc = mean(alive),
            .groups = "drop") 

## Visual ## 

# overall 

ggplot(cage_exp_spider, aes(x = spider_present, y = perc)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(sp ~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Spider Presence", 
       y = "Survival Proportion", 
       title = "Survival across Spider Presence") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

ggplot(cage_exp_spider %>% 
         filter(sp == "ach"), aes(x = spider_present, y = perc)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(trt ~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Spider Presence", 
       y = "Survival Proportion", 
       title = "Survival across Spider Presence") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

ggplot(cage_exp_spider %>% 
         filter(sp == "apt"), aes(x = spider_present, y = perc)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(trt ~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Spider Presence", 
       y = "Survival Proportion", 
       title = "Survival across Spider Presence") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

cage_exp_spider %>% 
  filter(round == 2) %>% 
  count(spider_present)

## model 

spiders <- glmmTMB(perc ~ sp * spider_present * burn + (1|block),
                              data = cage_exp_spider, family = "ordbeta")

plot(simulateResiduals(spiders))
summary(spiders)
Anova(spiders)
emmeans(spiders, pairwise ~ spider_present|sp, type = "response")




spiders_ach <- glmmTMB(perc ~ high_low * spider_present * burn + (1|block),
                   data = cage_exp_spider %>% 
                     filter(sp == "ach"), family = ordbeta())

plot(simulateResiduals(spiders_ach))
summary(spiders_ach)
Anova(spiders_ach)

emmeans(spiders_ach, pairwise ~ high_low|spider_present type = "response")





############################# PART 2 OF EXPERIMENT ###########################################

#### H2 - Density dependence in both grasshopper species would be weaker in burned than unburned plots ####

## Visuals ##

# Graph 

ggplot(cage_exp_dd %>% 
         filter(dep == "monoculture", round == 4), # round 4 is the second stocking event
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
                filter(dep == "monoculture", round == 4), family = "ordbeta")

plot(simulateResiduals(DD))

summary(DD)
Anova(DD)
emmeans(DD,pairwise ~ high_low|burn|sp, type = "response")

########## No effects of burn treatment on density dependence ###

#### H3a - Frequency dependence would  be weaker in burned vs. unburned plots ####

## Visuals ##

# Adjusting data set for FD analysis 

ggplot(cage_exp_fd %>% 
         filter(round == 4), aes(x = trt, y = perc)) +
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

FD <- glmmTMB(perc ~ burn * sp * trt + (1|block), data = cage_exp_fd, family = "ordbeta")

plot(simulateResiduals(FD))

summary(FD)
Anova(FD)
emmeans(FD,pairwise ~ trt|burn|sp, type = "response")

########## No effects of burn treatment on frequency dependence ###

#### H3b - Grass composition and abundance may favor the grass-specialist ####

## Visuals ##

surv_summary_2 <- cage_exp_surv %>% 
  filter(round == 4) %>% 
  group_by(strip, block, cage, trt, dep, sp, round, burn, high_low) %>% 
  summarise(
    perc_survival = mean(alive),
    days = mean(days),
    dens = n_distinct(ind),
    .groups = "drop"
  )  

surv_plant_2 <- surv_summary %>% 
  left_join(
    plant_summary %>% select(cage, grass_perc, grass_forb_ratio),
    by = "cage"
  )

ggplot(surv_plant_2, aes(x = grass_perc, y = perc_survival)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_grid(sp ~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Grass Percentage", 
       y = "Survival Proportion", 
       title = "Survival across grass abundance") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

ggplot(surv_plant_2, aes(x = grass_perc, y = perc_survival)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_grid(sp ~ trt) + 
  theme_bw(base_size = 20) + 
  labs(x = "Grass Percentage", 
       y = "Survival Proportion", 
       title = "Survival across grass abundance") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

## Model ##

grassdominance <- glmmTMB(perc_survival ~ grass_perc * sp * burn + (1|block), data = surv_plant_2, family = "ordbeta")

plot(simulateResiduals(grassdominance))

summary(grassdominance)
Anova(grassdominance)
emmeans(grassdominance,pairwise ~ burn|sp, type = "response")

# NO effects

############### HOW ABOUT ONLY ACHARUM?

ggplot(surv_plant_2 %>% 
         filter(sp == "ach"),
       aes(x = grass_perc, y = perc_survival)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_grid(trt ~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Frequency Treatment", 
       y = "Survival Proportion", 
       title = "Grasshoppers mixture survival in burned vs. unburned") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

## Model ##

grassdominance_ach <- glmmTMB(perc_survival ~ grass_perc * trt * burn + (1|block),
                              data = surv_plant_2 %>% 
                                filter(sp == "ach"), family = "ordbeta")

plot(simulateResiduals(grassdominance_ach))

summary(grassdominance_ach)
Anova(grassdominance_ach)
emmeans(grassdominance_ach,pairwise ~ burn|trt, type = "response")



#### ADDITIONAL: SPIDER PREDATION DATA ####

cage_exp_spider <- cage_exp %>% 
  mutate(spider_present = if_else(!is.na(spider) & spider == 1, 1, 0)) %>% 
  group_by(cage) %>% 
  summarise(
    spider_present = max(spider_present, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    spider_present = factor(spider_present, levels = c(0, 1), labels = c("no", "yes"))
  ) %>% 
  filter(cage <= 60) %>% 
  select(c(spider_present, cage))

cage_exp_spider <- cage_exp_surv %>% 
  left_join(exp_cage_spider, by = "cage") %>% 
  filter(round == 4) %>% 
  group_by(strip, block, cage, trt, dep, sp, round, burn, high_low, spider_present) %>% 
  summarize(perc = mean(alive),
            .groups = "drop") 

## Visual ## 

# overall 

ggplot(cage_exp_spider, aes(x = spider_present, y = perc)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(sp ~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Spider Presence", 
       y = "Survival Proportion", 
       title = "Survival across Spider Presence") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

ggplot(cage_exp_spider %>% 
         filter(sp == "ach"), aes(x = spider_present, y = perc)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(trt ~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Spider Presence", 
       y = "Survival Proportion", 
       title = "Survival across Spider Presence") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

ggplot(cage_exp_spider %>% 
         filter(sp == "apt"), aes(x = spider_present, y = perc)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(trt ~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Spider Presence", 
       y = "Survival Proportion", 
       title = "Survival across Spider Presence") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

cage_exp_spider %>% 
  filter(round == 4) %>% 
  count(spider_present)

## model 

spiders <- glmmTMB(perc ~ sp * spider_present * burn + (1|block),
                   data = cage_exp_spider, family = "ordbeta")

plot(simulateResiduals(spiders))
summary(spiders)
Anova(spiders)
emmeans(spiders, pairwise ~ spider_present|sp, type = "response")




spiders_ach <- glmmTMB(perc ~ high_low * spider_present * burn + (1|block),
                       data = cage_exp_spider %>% 
                         filter(sp == "ach"), family = ordbeta())

plot(simulateResiduals(spiders_ach))
summary(spiders_ach)
Anova(spiders_ach)

emmeans(spiders_ach, pairwise ~ high_low|spider_present, type = "response")
