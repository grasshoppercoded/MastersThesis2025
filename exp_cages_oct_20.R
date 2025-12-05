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

# adding burned vs. unburned 

cages <- data1 %>% 
  mutate(burn = case_when((strip %in% c(1,3,5) ~ "b"),
                          (strip %in% c(2,4,6) ~ "u"))) %>% 
  filter(alive %in% c(0, 1)| is.na(alive)) %>% 
  drop_na(strip)
  

unique(cages$alive)

#---- H1: GRASSHOPPERS SURVIVED MORE IN BURNED TREATMENT.APTENOPEDES DIED MORE ----

##---- plot ----

## 1.1 regression of mixes and monocultures in burned vs. unburned over time

ggplot(cages, aes(x = days, y = alive, color = sp)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(~burn) 

## model 

m1.1 <- glmmTMB(alive ~ days * sp * burn + (1|block),
                data = cages,
                family = "binomial")

Anova(m1.1)

plot(simulateResiduals(m1.1))

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
  group_by(strip, block, cage, trt, sp, date, type, round, burn) %>% 
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





