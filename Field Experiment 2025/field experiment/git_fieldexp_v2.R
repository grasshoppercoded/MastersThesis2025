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

data1 <- read.csv("cagestock_summ25_final.csv") ### grasshopper survival data
head(data1)
summary(data1)

data2 <- read.csv("plants_summ25_relativeabundance_fixed.csv") ### plant abundance data 
head(data2)
summary(data2)

#---- CLEANING DATA ----

### survival data ###

# create days since first stocking 

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
  mutate(high_low = case_when((trt %in% c("ach_high", "apt_high") ~ "high"), 
                              (trt %in% c("ach_low", "apt_low") ~ "low"))) %>% 
  filter(alive %in% c(0, 1) | is.na(alive)) %>%
  mutate(alive = as.numeric(alive)) %>% 
  select(c(strip, block, cage, round, dep, trt, high_low, burn, sp, ind, alive,days)) %>% 
  drop_na(strip)

unique(cages$alive)

### plant abundance data ###

# adding plant relative abundance data 

data2 <- data2 %>% 
  mutate(plant = str_trim(plant)) %>% 
  mutate(
    veg = case_when(
      plant %in% c("wide", "cent", "dican", "bb") ~ "grass",
      plant %in% c("silky", "dog", "other") ~ "forb")) %>% 
  select(-c(notes, bare)) %>% 
  filter(plant != "", round == 1) %>% 
  group_by(cage, b_u, veg, round) %>% 
  summarise(total = sum(perc)) %>% 
  pivot_wider(names_from = veg, values_from = total) %>% 
  mutate(grass_forb_ratio = grass / forb) %>% 
  mutate(grass_perc = (grass/(grass + forb))*100)

# make new dataframe of survival + relative abundance data 

surv_plant <- data2 %>% 
  left_join(cages, by = "cage") %>% 
  select(c(strip, block, cage, trt, dep, sp, burn, high_low, grass_perc, grass_forb_ratio)) %>%  
  group_by(strip, block, cage, trt, dep, sp, burn, high_low) %>% 
  summarize(perc = mean(alive),
            days = mean(days), 
            dens = n_distinct(ind),
            .groups = "drop") 

head(cages)
print(cages)
str(cages)

# ---- DATA VISUALIZATION FOR HYPOTHESES ----

#### Hypothesis 2 - A. carinatum survived more in the burned treatment (first 3 weeks) ####

# make data frame just for the initial 3 weeks, with proportion survival

R2 <- cages %>% 
  filter(round == 2) %>% 
  group_by(strip, block, cage, trt, dep, sp, round, burn, high_low) %>% 
  summarize(perc = mean(alive),
            days = mean(days), 
            dens = n_distinct(ind),
            .groups = "drop") 

str(cages$alive)
head(R2)

ggplot(R2 %>% 
         filter(sp == "ach"), 
       aes(x = sp, y = perc)) +
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

head(R2)

R2_ach <- glmmTMB(perc ~ burn * high_low + (1|block), data = R2 %>% 
                    filter(dep == "monoculture", sp == "ach"), family = "ordbeta")

plot(simulateResiduals(R2_ach))

summary(R2_ach)
Anova(R2_ach)
emmeans(R2_ach,pairwise ~ high_low|burn)

# grass ratio 

R2_ach <- glmmTMB(perc ~ sp * burn * high_low + (1|block), data = R2 %>% 
                    filter(dep == "monoculture"), family = "ordbeta")

plot(simulateResiduals(R2_ach))

summary(R2_ach).                Z
Anova(R2_ach)
emmeans(R2_ach,pairwise ~ high_low|sp)
