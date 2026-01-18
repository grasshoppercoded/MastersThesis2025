#---- clean, download ----

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

# download data 

setwd("C:/Users/Lucia.naviasalva/OneDrive - University of Florida/R")

data1 <- read.csv("cagestock_summ25.csv")

head(data1)

#---- data processing ----

# fix days since stocking 

data1$date <- as.Date(data1$date, format = "%m/%d/%Y")

stock_date <- as.Date("6/20/2025", format = "%m/%d/%Y")
data1$days <- as.numeric(data1$date - stock_date)
head(data1)

# clean up data 

cages <- data1 %>%
  drop_na(alive) 

head(cages)
  
#---- making survival model? ----

survival <- glmmTMB(alive ~ trt * days + (1|cage)+ (1|block)+ (1|strip), 
                    family = binomial, data = cages) 
                    
# error with alive binomial 

cages$alive <- as.numeric(cages$alive)

survival <- glmmTMB(alive ~ trt * days + (1|cage)+ (1|block)+ (1|strip), 
                    family = binomial, data = cages) 
# cant run 

cages_ <- cages %>%
  filter(alive %in% c(0, 1) | is.na(alive))

survival <- glmmTMB(alive ~ trt * days + (1|cage)+ (1|block)+ (1|strip), 
                    family = binomial, data = cages_) 
summary(survival)

ggplot(cages_, mapping = aes(x = ))

