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

## ---- Data Exploration ----

str(ra)

# check for typos 

unique(ra$plant)

##---- Long Format ----

# create forb and grass, and clean data 

ra <- ra %>% 
  mutate(plant = str_trim(plant)) %>% 
  mutate(
    veg = case_when(
      plant %in% c("wide", "cent", "dican", "bb") ~ "grass",
      plant %in% c("silky", "dog", "other") ~ "forb")) %>% 
  select(-c(notes, bare)) %>% 
  filter(plant != "")

count(ra, plant, veg) #make sure all species are there 

##---- Wide Format ---- 

wide_ra <- ra %>% 
  pivot_wider(
    id_cols = c(cage, b_u, round),
    names_from = plant,
    values_from = perc
  )

wide_ra

#---- NMDS -----

##---- Create Dissimilarity Matrix ----

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

##---- Ploting NMDS ----

# plot nmds with burn vs. unburn and round 1 & round 2

ggplot(data.scores, aes (x = MDS1, y = MDS2, color = burned, shape = round)) + 
  geom_point(size = 3) + 
  stat_ellipse() + 
  labs(title = "NMDS of plants species across burned sites and rounds") + 
  theme_classic()

# simplified, no rounds

ggplot(data.scores %>% 
         filter(round == 2),
       aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = burned), size = 3) +
  stat_ellipse(aes(color = burned), alpha = 0.05, geom = "polygon") +
  theme_classic() +
  labs(title = "NMDS of Plant Community Across Burned Sites",
       color = "Burn Treatment") +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  )


#---- PCoA ----

# species-only matrix 
# vegra <- wide_ra %>% select(-c(cage, burn, round))

vegdist_ra <- vegdist(vegra, "bray")
vegdist_ra

# dissimilarity matrix 

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

##---- Plotting PCoA ----

# plot PCoA with treatment and rounds

ggplot(scores_df, aes(x = V1, y = V2,
                      color = burn,
                      shape = round)) +
  geom_point(size = 3) +
  labs(
    title = "PCoA of Plant Communities",
    x = paste0("PCoA 1 (", round(propVar[1] * 100), "%)"),
    y = paste0("PCoA 2 (", round(propVar[2] * 100), "%)")
  ) +
  theme_minimal()

# species scores 

species_ra <- wascores(cmd_ra$points[, 1:2], vegra)
species_df <- as.data.frame(species_ra)
species_df$species <- rownames(species_df)

# plot with plant species names

ggplot() +
  geom_point(data = scores_df, aes(x = V1, y = V2, color = burn, shape = round), size = 3) +
  geom_point(data = species_df, aes(x = V1, y = V2), color = "red") +
  geom_text_repel(data = species_df, aes(x = V1, y = V2, label = species), color = "red") +
  labs(
    title = "PCoA of Plant Communities with Species Scores",
    x = paste0("PCoA 1 (", round(propVar[1] * 100), "%)"),
    y = paste0("PCoA 2 (", round(propVar[2] * 100), "%)")
  ) +
  theme_minimal()

# facet-wrap of PCoA by round

ggplot(scores_df, aes(x = V1, y = V2, color = burn)) +
  geom_point(size = 3) +
  labs(
    title = "PCoA of Plant Communities by Round",
    x = paste0("PCoA 1 (", round(propVar[1] * 100), "%)"),
    y = paste0("PCoA 2 (", round(propVar[2] * 100), "%)")
  ) +
  facet_wrap(~ round) 

ggplot(scores_df, aes(x = V1, y = V2, color = burn)) +
  geom_point(size = 3) +
  facet_wrap(~ round) +
  labs(
    title = "PCoA of Plant Communities Across Burned and Unburned Plots and Survey Rounds",
    x = paste0("PCoA 1 (", round(propVar[1] * 100), "%)"),
    y = paste0("PCoA 2 (", round(propVar[2] * 100), "%)")
  ) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )


## ggplot with convex hulls OF NMDS !!!

hulls <- data.scores %>% 
  group_by(burned) %>% 
  slice(chull(MDS1, MDS2))

# plot with convex hulls

ggplot(data.scores, aes(x = MDS1, y = MDS2, color = burned)) +
  geom_polygon(data = hulls, aes(fill = burned, group = burned), alpha = 0.15, color = NA) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "NMDS of Plant Communities by Burn Treatment") +
  scale_color_manual(values = c("b" = "tomato3", "u" = "seagreen")) +
  scale_fill_manual(values = c("b" = "tomato3", "u" = "seagreen"))

#---- PERMANOVA ---- 

# bray curtis distance using your species matrix (vegra)

dis_traits <- vegdist(vegra, "bray")

# permanova testing burned vs unburned

set.seed(11)
perm_traits_burn <- adonis2(dis_traits ~ data.scores$burned,
                       permutations = 1000)

# view results

perm_traits_burn

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
perm_traits_2 <- adonis2(dis_traits ~ burned + round + burned:round,
                       data = data.scores, permutations = 1000,
                       by = "terms")

# view 2WAY-PERMANOVA results
perm_traits_2

# testing homogeneity of dispersion 

# test dispersion for burned
disp_burn <- betadisper(dis_traits, data.scores$burned)
permutest(disp_burn)

# test dispersion for round
disp_round <- betadisper(dis_traits, data.scores$round)
permutest(disp_round)

#---- GRASS VS. FORB DOMINANCE (just for visualization) ----

# define grass and forb species
grass_species <- c("wide", "cent", "dican", "bb")
forb_species  <- c("silky", "dog", "other")

# calculate total grass and forb abundance per cage
wide_ra_dom <- wide_ra %>% 
  mutate(
    grass_total  = rowSums(across(all_of(grass_species)),  na.rm = TRUE),
    forb_total   = rowSums(across(all_of(forb_species)),   na.rm = TRUE),
    total_abund  = grass_total + forb_total,
    grass_percent = grass_total / total_abund
  )

# classify dominance: >50% grass cover = grass dominant
wide_ra_dom <- wide_ra_dom %>% 
  mutate(
    dominance = ifelse(grass_percent > 0.5, "grass_dom", "forb_dom")
  )

wide_ra_dom$dominance <- factor(wide_ra_dom$dominance)

# add dominance to NMDS scores for plotting only
data.scores$dominance <- wide_ra_dom$dominance

# compute convex hulls for burn × dominance groups
hulls_dom <- data.scores %>% 
  group_by(burned, dominance) %>% 
  slice(chull(MDS1, MDS2))

# NMDS plot with dominance groups
ggplot(data.scores, aes(x = MDS1, y = MDS2)) +
  geom_polygon(
    data = hulls_dom,
    aes(fill = dominance, color = burned, 
        group = interaction(burned, dominance)),
    alpha = 0.15
  ) +
  geom_point(aes(color = burned, shape = dominance), size = 3) +
  theme_bw() +
  labs(title = "NMDS of Plant Communities: Grass vs Forb Dominance × Burn Treatment")

