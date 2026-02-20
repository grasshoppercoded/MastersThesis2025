#---- SETUP ----

rm(list = ls())

# download packages  

library(tidyverse)
library(emmeans)
library(car)
library(agridat)
library(DHARMa)
library(glmmTMB)
library(MuMIn)
library(performance)
library(easystats)
library(lme4)
library(survival)
library(survminer)
library(ggsurvfit)
library(readxl)
library(betareg)
library(usethis)
library(gitcreds)

# download data 

data1 <- read_csv("cagestock_summ25_final.csv") %>% 
  select(-notes)

head(data1)

summary(data1)

#---- CLEANING DATA ----

# fix days since stocking 

data1$date <- as.Date(data1$date, format = "%m/%d/%Y")

stock_date <- as.Date("6/20/2025", format = "%m/%d/%Y")
data1$days <- as.numeric(data1$date - stock_date)
head(data1)

# adding burned vs. unburned, mixtures and monocultures, proportion survival, and simplifying

cages <- data1 %>% 
  mutate(burn = case_when((strip %in% c(1,3,5) ~ "b"),
                          (strip %in% c(2,4,6) ~ "u"))) %>% 
  mutate(dep = case_when((trt %in% c("ach_high","apt_high","apt_low", "ach_low") ~ "monoculture"),
                         (trt %in% c("ach_66","ach_33") ~ "mixture"))) %>% 
  filter(alive %in% c(0, 1) | is.na(alive)) %>%
  select(c(strip, block, cage, round, dep, trt, burn, sp, ind, alive, days)) %>% 
  drop_na(strip)

unique(cages$alive)
head(cages)
print(cages)

# ---- H1: Grasshoppers survived more in burned treatment. Aptenopedes died more overall ----

## ---- 1.1 survival regression  ----

## ROUND 2 ONLY ##

### 1.1 survival in burned vs. unburned over time (round 2, after 3 initial weeks)

ggplot(cages %>% 
         filter(round %in% c(1,2)), aes(x = days, y = alive, color = sp)) +
  geom_jitter(height=.1) +
  geom_smooth(method = "lm") +
  facet_grid(~ burn) 

## ---- 1.1 survival boxplot ----

### 1.1 overall survival in burned vs. unburned using % survival (round 2)

# make data set with proportion survival 

DD <- cages %>% 
  group_by(strip, block, cage, trt, dep, sp, round, burn) %>% 
  filter(round == 2) %>% 
  summarise(perc = mean(alive),
            days = mean(days), 
            dens = n_distinct(ind),
            .groups = "drop")

ggplot(DD, 
       aes(x = sp, y = perc, color = sp)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Species", 
       y = "Survival Proportion", 
       title = "Overall survival of species in burned vs. unburned") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))
ggplot(DD, 
       aes(x = sp, y = perc, color = sp)) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = 0.2) +
  geom_jitter(width = 0.1, alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(~ burn) + 
  theme_bw(base_size = 26) + 
  labs(x = "Species", 
       y = "Survival Proportion", 
       title = "Overall survival of species in all treatments") + 
  theme(plot.title = element_text(hjust = 0.5, 
                                  face = "bold", 
                                  size = 28))



# ---- 1.2 Overall survival with all treatments ----

ggplot(DD, 
       aes(x = sp, y = perc, color = sp)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(dep ~ burn) + 
  theme_bw(base_size = 26) + 
  labs(x = "Species", 
       y = "Survival Proportion", 
       title = "Overall survival of species in all treatments") + 
  theme(plot.title = element_text(hjust = 0.5, 
                                  face = "bold", 
                                  size = 28))

ggplot(DD, 
       aes(x = sp, y = perc, color = sp)) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = 0.2) +
  geom_jitter(width = 0.1, alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(dep ~ burn) + 
  theme_bw(base_size = 26) + 
  labs(x = "Species", 
       y = "Survival Proportion", 
       title = "Overall survival of species in all treatments") + 
  theme(plot.title = element_text(hjust = 0.5, 
                                  face = "bold", 
                                  size = 28))




H1 <- glmmTMB(perc ~ sp * burn + (1|block), # can remove the 3-way interactions because its non-sig
                     data = DD,
                     family = "ordbeta") # ordbeta can handle 0 and 1 and proportions between

plot(simulateResiduals(H1))

summary(H1)
Anova(H1)
emmeans(H1, pairwise ~ sp|burn, type = "response")

# ---- H2: Achurum performed best in mixtures ----

## ---- 2.1 treatment boxplot ----

## ---- 2.1 survival in burned vs. unburned across treatments (round 2) of Achurum 

ggplot(DD %>% 
         filter(sp == "ach"), 
       aes(x = sp, y = perc, color = sp)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(dep ~ burn) +
  theme_bw(base_size = 22) + 
  labs(x = "A. carinatum", 
       y = "Survival Proportion", 
       title = "Survival of A. carinatum in all treatments") + 
  theme(plot.title = element_text(hjust = 0.5, 
                                  face = "bold", 
                                  size = 28))

ggplot(DD %>% 
         filter(sp == "ach"), 
       aes(x = sp, y = perc, color = sp)) +
  stat_summary(fun = mean, geom = "point", size = 3.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = 0.2) +
  geom_jitter(width = 0.1, alpha = 0.4) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(dep ~ burn) +
  theme_bw(base_size = 22) + 
  labs(x = "A. carinatum", 
       y = "Survival Proportion", 
       title = "Survival of A. carinatum in all treatments") + 
  theme(plot.title = element_text(hjust = 0.5, 
                                  face = "bold", 
                                  size = 28))

## model 

H2 <- glmmTMB(perc ~ dep * burn + (1|block), # can remove the 3-way interactions because its non-sig
                     data = DD %>% filter(sp == "ach"),
                     family = "ordbeta") # ordbeta can handle 0 and 1 and proportions between

summary(H2)
Anova(H2)
emmeans(H2, pairwise ~ dep|burn, type = "response") # need to back-transform

plot(simulateResiduals(H2))


# ---- H3: Grasshoppers experienced more density dependence in monocultures of the burned plot ----

## ---- 3.1 overall density dependence in mixtures and monocultures ----

ggplot(DD,aes(x = factor(dens), y = perc, color = dep)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(~ burn) +
  theme_bw(base_size = 22) + 
  labs(x = "Densities and frequencies", 
       y = "Survival Proportion", 
       title = "Density and frequency dependence of both species in all treatments") + 
  theme(plot.title = element_text(hjust = 0.2, 
                                  face = "bold", 
                                  size = 23))

## ----- 3.2 achurum density dependence in mixtures and monocultures ----

ggplot(DD %>%
         filter(sp == "ach", dep == "monoculture"),
       aes(x = factor(dens), y = perc, color = dep)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(~ burn) +
  theme_bw(base_size = 22) + 
  labs(x = "Densities of Achurum", 
       y = "Survival Proportion", 
       title = "Density dependence of A. carinatum in monocultures") + 
  theme(plot.title = element_text(hjust = 0.2, 
                                  face = "bold", 
                                  size = 22))
ggplot(DD %>%
         filter(sp == "ach", dep == "monoculture"),
       aes(x = factor(dens), y = perc, color = dep)) +
  stat_summary(fun = mean, geom = "point", size = 3.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = 0.2) +
  geom_jitter(width = 0.1, alpha = 0.4) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(~ burn) +
  theme_bw(base_size = 22) + 
  labs(x = "Densities of Achurum", 
       y = "Survival Proportion", 
       title = "Density dependence of A. carinatum in monocultures") + 
  theme(plot.title = element_text(hjust = 0.2, 
                                  face = "bold", 
                                  size = 22), 
        legend.position = "none")


## ----- 3.2 achurum % frequency in mixtures and monocultures ----

ach_DD <- DD %>%
  mutate(ach_dens = case_when((trt == "ach_high"~ "100"),
                              (trt == "ach_low"~ "100"),
                              (trt == "ach_66"~ "66"),
                              (trt == "ach_33"~ "33")),
         ach_dens = factor(ach_dens, levels = c("33", "66", "100")))

ggplot(ach_DD %>%
         filter(sp == "ach"),
       aes(x = factor(ach_dens), y = perc, color = dep)) +
  stat_summary(fun = mean, geom = "point", size = 3.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = 0.2) +
  geom_jitter(width = 0.1, alpha = 0.4) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(~ burn) +
  theme_bw(base_size = 22) + 
  labs(x = "Percent frequency of Achurum", 
       y = "Survival Proportion", 
       title = "Frequency dependence of Achurum in all treatments") + 
  theme(plot.title = element_text(hjust = 0.2, 
                                  face = "bold", 
                                  size = 22), 
        legend.position = "none")

## ----- 3.3 aptenopedes density dependence in mixtures and monocultures ----

ggplot(DD %>%
         filter(sp == "apt"),
       aes(x = factor(dens), y = perc, color = dep)) +
  stat_summary(fun = mean, geom = "point", size = 3.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = 0.2) +
  geom_jitter(width = 0.1, alpha = 0.4) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(~ burn)

ggplot(DD %>%
         filter(sp == "apt"),
       aes(x = factor(dens), y = perc, color = dep)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(~ burn) + 
  labs(x = "Density of Aptenopedes", 
       y = "Survival Proportion", 
       title = "Density dependence of Aptenopedes in all treatments") + 
  theme_bw(base_size = 22) + 
  theme(plot.title = element_text(hjust = 0.2, 
                                  face = "bold", 
                                  size = 22))

H4 <- glmmTMB(perc ~ dep * burn + (1|block), # can remove the 3-way interactions because its non-sig
              data = DD %>% filter(sp == "apt"),
              family = "ordbeta") # ordbeta can handle 0 and 1 and proportions between

summary(H4)
Anova(H4)
emmeans(H4, pairwise ~ dep|burn, type = "response") # need to back-transform

plot(simulateResiduals(H2))

## ---- 3.4 comparison of both species density dependence 

ggplot(DD, aes(x = factor(dens), y = perc, color = dep)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(sp ~ burn) +
  theme_bw(base_size = 22) + 
  labs(x = "Densities", 
       y = "Survival Proportion", 
       title = "DD of both species in all treatments") + 
  theme(plot.title = element_text(hjust = 0.5, 
                                  face = "bold", 
                                  size = 28))

## ---- H3 model ----

densxsp <- glmmTMB(perc ~ sp * dep + dens + burn + (1|block), # can remove the 3-way interactions because its non-sig
                   data = DD,
                   family = "ordbeta")

plot(simulateResiduals(densxsp))

summary(densxsp)
Anova(densxsp)

emmeans(densxsp, pairwise ~ perc|dens, type = "response") # another way to look at contrasts
emmeans(densxsp, pairwise ~ sp*dens, type = "response")

# ---- H4: Achurum performed better overall in grass-dominated plots ----

# upload nutrient/plant dataset 

data2 <- read_csv("Field Experiment 2025/nutrient content analyses/plants_summ25_relativeabundance_fixed.csv")

head(data2)
summary(data2)

ra <- data2 %>% 
  mutate(plant = str_trim(plant)) %>% 
  mutate(
    veg = case_when(
      plant %in% c("wide", "cent", "dican", "bb") ~ "grass",
      plant %in% c("silky", "dog", "other") ~ "forb")) %>% 
  select(-c(notes, bare)) %>% 
  filter(plant != "")

# clean data to be joined, grass to forb ratio created 

surv_join <- ra %>%
  filter(round == 1) %>%
  group_by(cage, b_u, veg, round) %>% 
  summarise(total = sum(perc)) %>% 
  pivot_wider(names_from = veg, values_from = total) %>% 
  mutate(grass_forb_ratio = grass / forb)
 

# add grass % 

surv_join <- surv_join %>% 
  mutate(grass_perc = (grass/(grass + forb))*100)

# left join of only round 4? 

surv_plants <- surv_join %>% 
  left_join(DD, by = "cage")

## ---- 4.1 regression of overall achurum survival and grass percentage ---- 

ggplot(surv_plants %>% 
         filter(sp == "ach"),
       aes(x = grass_perc, y = perc)) +
  #geom_jitter(width = 0.05, height = 0.05) +
  geom_point() +
  geom_smooth(method = "lm") # can add family=beta

## ---- 4.2 achurum survival and grass % based on frequency dependence regression ----

ggplot(surv_plants %>% 
         filter(sp == "ach"),
       aes(x = grass_perc, y = perc, color=burn)) +
  geom_jitter(width = 0.05, height = 0.05) +
  geom_smooth(method = "lm") + 
  facet_grid(~ dep)

## ---- 4.3 achurum survival and grass % based on density dependence regression ----

ggplot(surv_plants %>% 
         filter(sp == "ach"),
       aes(x = grass_perc, y = perc, color = dep)) +
  geom_jitter(width = 0.05, height = 0.05) +
  geom_smooth(method = "lm") + 
  facet_grid(~ dens)












#------------------------------------------------------------------------ 

## 1.2 regression of monoculture only survival in burned vs unburned 

ggplot(cages %>% 
         filter(trt %in% c("ach_high", "ach_low", "apt_high", "apt_low")), 
       mapping = aes(x = days, y = alive, color = trt)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(burn~sp) 

# okay, aptenopedes clearly had more mortality 

## model 

m1.2 <- glmmTMB(alive ~ days * sp * trt * burn + (1|block) + (1|cage),
                data = cages %>% 
                  filter(trt %in% c("ach_high", "ach_low", "apt_high", "apt_low")),
                family = "binomial")

Anova(m1.2)

plot(simulateResiduals(m1.2))


# ---- H2: DENSITY DEPENDENCE USING PROPORTIONS ----

# cleaning Data

head(cages)

DD <- cages %>% 
  group_by(strip, block, cage, trt, sp, round, burn) %>% 
  summarise(perc = mean(alive),
            days = mean(days), 
            dens = n_distinct(ind),
            .groups = "drop") %>% 
  mutate(freq = case_when((trt %in% c("ach_low", "ach_high", "apt_high", "apt_low") ~ "mono"),
                          (trt %in% c("ach_66", "ach_33") ~ "mix")))

unique(cages$ind)
## ---- plot proportions. alive using days ----

# regression of all treatments, regardless of b vs. u

ggplot(DD, aes(x = days, y = perc, color = trt)) +
  geom_point() +
  geom_smooth() +
  facet_grid(trt ~ sp)

# regression of just mix vs monocultures using round 2 only 

cages %>% 
  filter(round == 2) %>% 
  summary()

# round 2 
ggplot(DD %>% 
         drop_na(perc) %>% 
         filter(round == 1 | round == 2),
       aes(x = days, y = perc, color = burn)) +
  geom_point() +
  geom_smooth(method = "lm", method.args = list(family = "binomial"), se = T) +
  facet_grid(freq ~ sp) 

# regression of mix vs monocultures over densities 

ggplot(DD %>% 
         drop_na(perc) %>% 
         filter(round == 4), aes(x = dens, y = perc)) +
  geom_point() +
  geom_smooth(method = "lm", method.args = list(family = "binomial"), se = T) +
  facet_grid(round ~ sp)

#round 5 

ggplot(DD %>% 
         drop_na(perc) %>% 
         filter(round == 5), aes(x = dens, y = perc)) +
  geom_point() +
  geom_smooth(method = "lm", method.args = list(family = "binomial"), se = T) +
  facet_grid(round ~ sp)


# boxplot style freq

ggplot(DD %>%
         filter(round == 4), aes(x = interaction(dens,freq), y = perc, color = as.factor(freq))) +
  geom_boxplot() +
  geom_point() +
  facet_grid(round ~ sp)

ggplot(DD %>%
         filter(round == 4), aes(x = interaction(dens,freq), y = perc, color = as.factor(freq))) +
  geom_boxplot() +
  geom_point() +
  facet_grid(burn ~ sp)

ggplot(DD %>%
         filter(round == 2), aes(x = interaction(dens,freq), y = perc, color = as.factor(freq))) +
  geom_boxplot() +
  geom_point() +
  facet_grid(burn ~ sp)

#####

ggplot(DD %>%
         filter(round == 2, sp == "ach", dens == 6), aes(x = burn, y = perc, color = interaction(dens,freq))) +
  geom_boxplot() +
  geom_point() 




unique(DD$perc)
unique(DD$dens)


# ---- MODELL ----

# density 

dens_burn_mod <- glmmTMB(perc ~ freq * burn + (1|block), family = "ordbeta", data = DD %>% 
                           filter(round == 4, sp == "ach", dens == 6))

summary(dens_burn_mod)

Anova(dens_burn_mod)

emmeans(dens_burn_mod, pairwise~ freq|burn, type = "response")

emdat <- emmeans(dens_burn_mod, ~ freq|burn, type = "response") %>%  
  as.data.frame()

#jitterdodge
#geom_pointrange

# density by burn only 


dens_burn_mod <- glmmTMB(perc ~ dens * burn + (1|block), family = "ordbeta", data = DD %>% 
                           filter(round == 4, sp == "ach"))
summary(dens_burn_mod)
Anova(dens_burn_mod)
emmeans(dens_burn_mod)

emtrends(dens_burn_mod, ~ burn, var = "dens", infer = T)

# freq 

freq_mod <- glmmTMB(perc ~ dens * sp * freq + (1|block), family = "ordbeta", data = DD %>% 
                      filter(round == 4))
Anova(freq_mod)

# ---- B vs. U ---- 

# boxplot 

ggplot(DD %>%
         filter(round == 4, sp == "ach", freq == "mono"), aes(x = interaction(dens,freq), y = perc, color = as.factor(freq))) +
  geom_boxplot() +
  geom_point() +
  facet_grid(~burn)

unique(DD$dens)

# boxplot of density and burn treatment only  

ggplot(DD %>%
         filter(round == 4, sp == "ach"), aes(x = as.factor(freq), y = perc)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(~ burn)

# regression 

ggplot(DD %>% 
         drop_na(perc) %>% 
         filter(round == 4, sp == "ach")) +
  geom_point(aes(x = dens, y = perc, color = as.factor(freq))) +
  geom_smooth(aes(x = dens, y = perc), 
              method = "glm", method.args = list(family = "binomial"),
              se = T) +
  facet_grid(burn ~ sp)

######## ------------- BODY SIZE/LENGTH ##################

ggplot(DD %>%
         drop_na(sex) %>% 
         filter(sex != "u") %>% 
         filter(round == 6, sp == "ach"), aes(x = interaction(dens,freq), y = length)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(sex ~ burn)

library(zip)





