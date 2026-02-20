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

data1 <- read.csv("cagestock_summ25_final.csv")

head(data1)

summary(data1)

#---- CLEANING DATA ----

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
  filter(alive %in% c(0, 1) | is.na(alive)) %>%
  select(c(strip, block, cage, round, dep, trt, burn, sp, ind, alive, days)) %>% 
  drop_na(strip)

unique(cages$alive)
head(cages)
print(cages)

# ---- DATA VISUALIZATION FOR HYPOTHESES ----

#### Hypothesis 2 - A. carinatum survived more in the burned treatment (first 3 weeks) ####

# make data frame just for the initial 3 weeks, with proportion survival

R2 <- cages %>% 
  group_by(strip, block, cage, trt, dep, sp, round, burn) %>% 
  filter(round == 2) %>% 
  summarize(perc = mean(alive),
            days = mean(days), 
            dens = n_distinct(ind),
            .groups = "drop")

head(R2)

ggplot(R2, 
       aes(x = sp, y = perc) %>% 
         filter(sp == "ach")) +
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





