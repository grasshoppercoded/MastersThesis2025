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

# downloading git 
  
use_git_config(user.name = "Lucia", user.email = "lucia.naviasalva@ufl.edu")
  
# download data 

data1 <- read_csv("field experiment/cagestock_summ25_final.csv") %>% 
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

# ---- H1: Grasshoppers survived more in burned treatment.Aptenopedes died more overall ----

## ---- 1.1 survival regression  ----

## ROUND 2 ONLY ##

### 1.1 survival in burned vs. unburned over time (round 2, after 3 initial weeks)


ggplot(cages %>% 
         filter(round %in% c(1,2)), aes(x = days, y = alive, color = sp)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  facet_grid(~ burn) 

## ---- 1.1 survival boxplot ----

### 1.1 survival in burned vs. unburned using % survival (round 2)

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
  facet_grid(~ burn)

# ---- H2: Achurum performed best in mixtures ----

## ---- 2.1 treatment boxplot ----
## ---- 2.1 survival in burned vs. unburned across treatments (round 2)

ggplot(DD, 
       aes(x = sp, y = perc, color = sp)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(dep ~ burn)

## model 

survxburn <- glmmTMB(perc ~ sp * dep * burn + (1|block),
                data = DD,
                family = "binomial")

summary(survxburn)
Anova(survxburn)
emmeans(survxburn, pairwise ~ sp)
emmeans(survxburn, pairwise ~ sp|dep)

plot(simulateResiduals(survxburn))

# ---- H3: Achurum experienced more density dependence in monocultures of the burned plot

## ---- achurum densities in mixtures and monocultures

ggplot(DD,aes(x = factor(dens), y = perc, color = dep)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(~ burn)

# make data with percent achurum 

ach_DD <- DD %>% 
  mutate(ach_dens = case_when((trt == "ach_high"~ "100"),
                         (trt == "ach_low"~ "100"),
                         (trt == "ach_66"~ "66"),
                         (trt == "ach_33"~ "33"))) 

ggplot(ach_DD %>% 
         filter(sp == "ach"), 
       aes(x = factor(ach_dens), y = perc, color = dep)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(~ burn)







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





