#---- SETUP ----

rm(list = ls())

# packages 

library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(grid)
library(ca)
library(stringr)

# loading data 

ra <- read.csv("plants_summ25_relativeabundance_fixed.csv")

# explore 

str(ra)

# check for typos 

unique(ra$plant)

# fix burned vs. unburned, create forb and grass, and clean data - LONG FORMAT

ra <- ra %>% 
  mutate(plant = str_trim(plant)) %>% 
  mutate(
    veg = case_when(
      plant %in% c("wide", "cent", "dican", "bb") ~ "grass",
      plant %in% c("silky", "dog", "other") ~ "forb")) %>% 
  select(-c(notes, bare)) %>% 
  filter(plant != "")

# check for duplicate columns 

ra %>% 
  summarise(n = n(), .by = c(cage, b_u, veg, plant)) %>% 
  filter(n > 1) 

count(ra, plant, veg) #mkae sure all species are there 

# make data into wide format 

wide_ra <- ra %>% 
  pivot_wider(
    id_cols = c(cage, b_u, round),
    names_from = plant,
    values_from = perc
  )

wide_ra

#---- ATTEMPT NMDS? -----

# create dissimilarity matrix

vegra <- wide_ra %>% 
  select(-c(b_u, cage, round))

radist <- vegdist(vegra, "bray")

nmdsra <- metaMDS(radist, k = 2, trace = T)
nmdsra

#or nmdsra <- metaMDS(vegra, distance = "bray", k = 2, trymax = 200)

stressplot(nmdsra)

# grouping variables 

# NMDS with ggplot 

data.scores <- as.data.frame(nmdsra$points) # axes 
data.scores$burned <- wide_ra$b_u
data.scores$cage <- wide_ra$cage
data.scores$round <- factor(wide_ra$round)

str(data.scores) #making sure round is categorical

# plot nmds 

ggplot(data.scores, aes (x = MDS1, y = MDS2, color = burned, shape = round)) + 
  geom_point(size = 3) + 
  geom_text_repel(aes(label = cage)) + 
  stat_ellipse() + 
  labs(title = "NMDS of plants species across burned sites") + 
  theme_classic()

# chatgpt ggplot 

ggplot(data.scores, aes(x = MDS1, y = MDS2)) +
  
  # points colored by burn
  geom_point(aes(color = burned), size = 3) +
  
  # ellipses filled by round, outlined by burn
  stat_ellipse(aes(color = burned, fill = round), 
               alpha = 0.1, 
               geom = "polygon") +
  
  geom_text_repel(aes(label = cage)) +
  
  theme_classic() +
  labs(title = "NMDS with burn (color) and round (ellipse fill)")

#---- ATTEMPT PoCA?

# species-only matrix 
# vegra <- wide_ra %>% select(-c(cage, burn, round))

vegdist_ra <- vegdist(vegra, "bray")
vegdist_ra

# similarity matrix 

cmd_ra <- cmdscale(vegdist_ra, k = 5, eig = TRUE)
cmd_ra

# extract eigenvalues 

eigenvalues <- cmd_ra$eig[1:5]
propVar <- eigenvalues / sum(eigenvalues)
cumVar <- cumsum(propVar)

PCoA_Table <- tibble(
  Axis = 1:5,
  Eigenvalue = eigenvalues,
  Proportion = propVar,
  Cumulative = cumVar
)

print(PCoA_Table)

# scree plot 

ggplot(PCoA_Table, aes(x = Axis, y = Eigenvalue)) +
  geom_point(size = 3) +
  geom_line() +
  labs(title = "PCoA Scree Plot",
       x = "Axis",
       y = "Eigenvalue") +
  theme_minimal()

# create scores 

scores_df <- as.data.frame(cmd_ra$points[, 1:2])
scores_df$cage <- wide_ra$cage
scores_df$burn <- wide_ra$b_u
scores_df$round <- factor(wide_ra$round)   # round as factor

# plot! 

ggplot(scores_df, aes(x = V1, y = V2,
                      color = burn,
                      shape = round)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = cage)) +
  labs(
    title = "PCoA of Plant Communities",
    x = paste0("PCoA 1 (", round(propVar[1] * 100), "%)"),
    y = paste0("PCoA 2 (", round(propVar[2] * 100), "%)")
  ) +
  theme_classic()

# species scores 

species_ra <- wascores(cmd_ra$points[, 1:2], vegra)
species_df <- as.data.frame(species_ra)
species_df$species <- rownames(species_df)

# plot with this 

ggplot() +
  geom_point(data = scores_df, aes(x = V1, y = V2, color = burn, shape = round), size = 3) +
  geom_text_repel(data = scores_df, aes(x = V1, y = V2, label = cage)) +
  geom_point(data = species_df, aes(x = V1, y = V2), color = "red") +
  geom_text_repel(data = species_df, aes(x = V1, y = V2, label = species), color = "red") +
  labs(
    title = "PCoA with Species Scores",
    x = paste0("PCoA 1 (", round(propVar[1] * 100), "%)"),
    y = paste0("PCoA 2 (", round(propVar[2] * 100), "%)")
  ) +
  theme_minimal()

# facet-wrap 

ggplot(scores_df, aes(x = V1, y = V2, color = burn)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = cage)) +
  labs(
    title = "PCoA of Plant Communities by Round",
    x = paste0("PCoA 1 (", round(propVar[1] * 100), "%)"),
    y = paste0("PCoA 2 (", round(propVar[2] * 100), "%)")
  ) +
  facet_wrap(~ round)

## chatgpt ggplot with convex hulls OF NMDS !!!

# compute convex hulls for each burned category

hulls <- data.scores %>% 
  group_by(burned) %>% 
  slice(chull(MDS1, MDS2))

# plot with convex hulls

ggplot(data.scores, aes(x = MDS1, y = MDS2, color = burned)) +
  geom_polygon(data = hulls, aes(fill = burned, group = burned), alpha = 0.15, color = NA) +
  geom_point(size = 2) +
  geom_text(aes(label = cage), vjust = -0.5) +
  theme_bw() +
  labs(title = "NMDS of Plant Communities by Burn Treatment") +
  scale_color_manual(values = c("b" = "tomato3", "u" = "seagreen")) +
  scale_fill_manual(values = c("b" = "tomato3", "u" = "seagreen"))

#---- PERMANOVA ---- 

# bray curtis distance using your species matrix (vegra)

dis_traits <- vegdist(vegra, "bray")

# permanova testing burned vs unburned

set.seed(11)
perm_traits <- adonis2(dis_traits ~ data.scores$burned, permutations = 1000)

# view results

perm_traits

# histogram of permuted F values

F_df <- data.frame(F_perm = permustats(perm_traits)$permutations)

ggplot(F_df, aes(x = F_perm)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  geom_vline(xintercept = perm_traits$F[1], color = "red", linewidth = 1.2) +
  labs(title = "Histogram of F statistics", x = "F", y = "Frequency") +
  theme_minimal()

# ---- TWO-WAY PERMANOVA ---- 

# make sure round is a factor
data.scores$round <- factor(data.scores$round)

# bray-curtis distance
dis_traits <- vegdist(vegra, "bray")

# run permanova with burn, round, and their interaction
set.seed(11)
perm_traits <- adonis2(dis_traits ~ burned + round + burned:round, data = data.scores, permutations = 1000, by = "terms")

# view results
perm_traits

# testing homogeneity of dispersion 

# test dispersion for burned
disp_burn <- betadisper(dis_traits, data.scores$burned)
permutest(disp_burn)

# test dispersion for round
disp_round <- betadisper(dis_traits, data.scores$round)
permutest(disp_round)

#---- DOING GRASS VS. FORB DOMINANCE ----

# define grass and forb species
grass_species <- c("wide", "cent", "dican", "bb")
forb_species  <- c("silky", "dog", "other")

# calculate total grass and forb abundance per cage

wide_ra_dom <- wide_ra %>% 
  mutate(
    grass_total = rowSums(across(all_of(grass_species)), na.rm = TRUE),
    forb_total  = rowSums(across(all_of(forb_species)), na.rm = TRUE)
  )

# add the dominance 

wide_ra_dom <- wide_ra_dom %>% 
  mutate(
    dominance = ifelse(grass_total > forb_total, "grass_dom", "forb_dom")
  )

wide_ra_dom$dominance <- factor(wide_ra_dom$dominance)

# adding dominance to my scores 

data.scores$dominance <- wide_ra_dom$dominance

# compute hulls by burn × dominance

hulls_dom <- data.scores %>% 
  group_by(burned, dominance) %>% 
  slice(chull(MDS1, MDS2))

ggplot(data.scores, aes(x = MDS1, y = MDS2)) +
  geom_polygon(
    data = hulls_dom,
    aes(fill = dominance, color = burned, group = interaction(burned, dominance)),
    alpha = 0.15
  ) +
  geom_point(aes(color = burned, shape = dominance), size = 3) +
  geom_text_repel(aes(label = cage)) +
  theme_bw() +
  labs(title = "NMDS of Plant Communities: Grass vs Forb Dominance × Burn Treatment")

#---- PERMANOVA DOMINANCE -----

data.scores$burned   <- factor(data.scores$burned)
data.scores$dominance <- factor(data.scores$dominance)

dist_dom_burn <- vegdist(vegra, method = "bray")

set.seed(11)
perm_dom_burn <- adonis2(
  dist_dom_burn ~ burned + dominance,
  data = data.scores,
  permutations = 1000, 
  by = "margin"
)

# use by margin! 

perm_dom_burn

# testing homogeinity 

disp_burn_dom <- betadisper(dist_dom_burn, data.scores$burned)
permutest(disp_burn_dom)

disp_dom <- betadisper(dist_dom_burn, data.scores$dominance)
permutest(disp_dom)
