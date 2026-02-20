#---- SETUP ----

rm(list = ls())

# download packages  

library(tidyverse)
library(emmeans)
library(car)
library(DHARMa)
library(glmmTMB)
library(MuMIn)
library(performance)
library(easystats)
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
