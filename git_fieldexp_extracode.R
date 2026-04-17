######### PART 1 OF EXP ###############

# H2 - DD 

# overall for DD in all species, basic graph  

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


# old and original DD model for all species, insignificant 

DD <- glmmTMB(perc ~ burn * high_low + sp + (1|block), data = cage_exp_dd %>% 
                filter(dep == "monoculture"), family = "ordbeta")

plot(simulateResiduals(DD))

summary(DD)
Anova(DD)
emmeans(DD,pairwise ~ high_low|burn|sp, type = "response")

# H3 - FD 

# overall basic FD Graph 

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
# overall nice graph FD 


# FD ACH with grass as a continous variable 

ggplot(surv_plant %>% 
         filter(sp == "ach", trt != "ach_low"),
       aes(x = grass_perc, y = perc_survival, color = trt)) +
  geom_point(size = 2.5,
             alpha = 0.7) +
  geom_smooth(method = "lm",
              se = TRUE,
              linewidth = 1.2) +
  facet_grid(trt ~ burn) + 
  scale_color_viridis_d(option = "magma", end = 0.85) +
  theme_bw(base_size = 20) + 
  labs(x = "Grass Percent", 
       y = "Survival Proportion", 
       title = "Ach mixture survival in burned vs. unburned",
       color = "Treatment") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))
# OVERALL species, fd, old model 

FD <- glmmTMB(perc ~ burn * trt + sp + (1|block), data = cage_exp_fd,
              family = "ordbeta")

plot(simulateResiduals(FD))

summary(FD)
Anova(FD)
emmeans(FD,pairwise ~ trt|burn|sp, type = "response")

# grass dependence on both species 
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

# same thing, but with each treatment 
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
# grass dominance on survival of both species 

grassdominance <- glmmTMB(perc_survival ~ grass_perc * sp * burn + (1|block),
                          data = surv_plant, family = "ordbeta")

plot(simulateResiduals(grassdominance))

summary(grassdominance)
Anova(grassdominance)
emmeans(grassdominance,pairwise ~ burn|sp, type = "response")

#grass dominance on survival of achurum, in all treatments

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


######## PART TWO OF EXP ###############

################ ORIGINAL MODEL for DD in both sp, in part 2 of exp.##

DD <- glmmTMB(perc ~ burn * high_low * sp * (1|block), data = cage_exp_dd %>% 
                filter(dep == "monoculture", round == 4), family = "ordbeta")

plot(simulateResiduals(DD))

summary(DD)
Anova(DD)
emmeans(DD,pairwise ~ high_low|burn|sp, type = "response")


#DD simple boxplot with proportion

ggplot(cage_exp_dd %>% 
         filter(dep == "monoculture", round == 5), # round 4 is the second stocking event
       aes(x = factor(high_low, levels = c("low", "high")), y = perc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.3) +
  facet_grid(sp ~ burn) + 
  theme_bw(base_size = 20) + 
  labs(x = "Density", 
       y = "Survival Proportion", 
       title = "Grasshoppers monoculture survival in burned vs. unburned") + 
  theme(plot.title = element_text(hjust = 0.4, face = "bold", size = 22))
graphics.off()

 #DD original plot with number alive 

ggplot(cage_exp_allrounds %>% #find non-linear model to fit into graphs 
         filter(dep == "monoculture"),
       aes(x = density, y = alive_n, color = burn)) + ###
  geom_point() +
  geom_smooth(alpha = 0.15, formula = y ~ 0 + x) +
  facet_grid(round ~ sp) +
  theme_bw(base_size = 20) +
  labs(x = "total_start",
       y = "alive_n",
       title = "Survival over time in the second part of the experiment") +
  theme(plot.title = element_text(hjust = 0.4,
                                  face = "bold",
                                  size = 22))

###### Bevertion & Ricker model for DD Ach 

# beverton-holt model

ach_mono <- cage_exp_allrounds %>%
  filter(dep == "monoculture", sp == "ach")

bh_r4 <- nls(alive_n ~ density * s0 / (1 + a * density),
             data = ach_mono %>% filter(round == 4),
             start = list(s0 = 0.8, a = 0.05))
bh_r5 <- nls(alive_n ~ density * s0 / (1 + a * density),
             data = ach_mono %>% filter(round == 5),
             start = list(s0 = 0.8, a = 0.05))
bh_r6 <- nls(alive_n ~ density * s0 / (1 + a * density),
             data = ach_mono %>% filter(round == 6),
             start = list(s0 = 0.8, a = 0.05))

new_r4 <- data.frame(density = seq(min(filter(ach_mono, round == 4)$density),
                                   max(filter(ach_mono, round == 4)$density),
                                   length.out = 100), round = 4)
new_r5 <- data.frame(density = seq(min(filter(ach_mono, round == 5)$density),
                                   max(filter(ach_mono, round == 5)$density),
                                   length.out = 100), round = 5)
new_r6 <- data.frame(density = seq(min(filter(ach_mono, round == 6)$density),
                                   max(filter(ach_mono, round == 6)$density),
                                   length.out = 100), round = 6)

new_r4$pred <- predict(bh_r4, newdata = new_r4)
new_r5$pred <- predict(bh_r5, newdata = new_r5)
new_r6$pred <- predict(bh_r6, newdata = new_r6)

bh_preds <- bind_rows(new_r4, new_r5, new_r6)

ggplot(ach_mono, aes(x = density, y = alive_n)) +
  geom_point() +
  geom_line(data = bh_preds, aes(x = density, y = pred)) +
  facet_wrap(~round) +
  theme_bw()

# Ricker model 

ach_mono <- cage_exp_allrounds %>%
  filter(dep == "monoculture", sp == "ach")

rk_r4 <- nls(alive_n ~ density * exp(r - a * density),
             data = ach_mono %>% filter(round == 4),
             start = list(r = -0.5, a = 0.05))

rk_r5 <- nls(alive_n ~ density * exp(r - a * density),
             data = ach_mono %>% filter(round == 5),
             start = list(r = -0.5, a = 0.05))

rk_r6 <- nls(alive_n ~ density * exp(r - a * density),
             data = ach_mono %>% filter(round == 6),
             start = list(r = -0.5, a = 0.05))

summary(rk_r4)
summary(rk_r5)
summary(rk_r6)

new_r4 <- data.frame(density = seq(min(filter(ach_mono, round == 4)$density),
                                   max(filter(ach_mono, round == 4)$density),
                                   length.out = 100),
                     round = 4)

new_r5 <- data.frame(density = seq(min(filter(ach_mono, round == 5)$density),
                                   max(filter(ach_mono, round == 5)$density),
                                   length.out = 100),
                     round = 5)

new_r6 <- data.frame(density = seq(min(filter(ach_mono, round == 6)$density),
                                   max(filter(ach_mono, round == 6)$density),
                                   length.out = 100),
                     round = 6)

new_r4$pred <- predict(rk_r4, newdata = new_r4)
new_r5$pred <- predict(rk_r5, newdata = new_r5)
new_r6$pred <- predict(rk_r6, newdata = new_r6)

rk_preds <- bind_rows(new_r4, new_r5, new_r6)

ggplot(ach_mono, aes(x = density, y = alive_n)) +
  geom_point() +
  geom_line(data = rk_preds, aes(x = density, y = pred)) +
  facet_wrap(~round) +
  theme_bw()

# model discussed we'd use with Dr. Hahn, ach DD, which we were going to go with 

grassdominance2_ach_dd <- glmmTMB(perc_survival ~ grass_perc + high_low + burn * (1|block),
                                  data = surv_plant_2 %>% 
                                    filter(sp == "ach", dep == "monoculture"), family = "ordbeta")

plot(simulateResiduals(grassdominance2_ach_dd))

summary(grassdominance2_ach_dd)
Anova(grassdominance2_ach_dd)
emmeans(grassdominance2_ach_dd,pairwise ~ high_low|burn, type = "response")

#AIC 

AIC(bh_r4, rk_r4)
AIC(bh_r5, rk_r5)
AIC(bh_r6, rk_r6)

#### frequency dependence, old models 

#model containing both species 
########## IGNORE Model ##

FD <- glmmTMB(perc ~ burn * sp * trt + (1|block), data = cage_exp_fd, family = "ordbeta")

plot(simulateResiduals(FD))

summary(FD)
Anova(FD)
emmeans(FD,pairwise ~ trt|burn|sp, type = "response")

# unrelated to FD, but survival overall of both species thruoghout rounds 

ggplot(cage_exp_allrounds %>% 
filter(dep == "mixture"),
aes(x = round, y = prop_survival, color = burn)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = 0.15) +
  facet_grid(~ sp) +
  theme_bw(base_size = 20) +
  labs(x = "Round",
       y = "Proportion survival",
       title = "Survival over time in the second part of the experiment") +
  theme(plot.title = element_text(hjust = 0.4,
                                  face = "bold",
                                  size = 22))

#ach fd orginanl model 


ach_2_FD_ignore <- glmmTMB(prop_survival ~ trt * burn * as.factor(round) + (1|block/cage), 
                           family = "betareg",
                           data = cage_exp_allrounds %>% 
                             filter(sp == "ach", trt != "ach_low"))

# ach FD chatgpt model, but did not work and not sure what's going on 

# ACH FD without prop survival 

cage_exp_allrounds <- cage_exp_allrounds %>%
  select(-density) %>%
  left_join(start_density, by = c("strip", "block", "cage", "sp")) %>%
  rename(density = density_true)
cage_exp_allrounds %>%
  filter(sp == "ach", trt %in% c("ach_33", "ach_66", "ach_high"), alive_n > density)

ach_2_FD_binom <- glmmTMB(cbind(alive_n, density - alive_n) ~ trt * burn * as.factor(round) + (1 | block/cage),
                          family = binomial,
                          data = cage_exp_allrounds %>%
                            filter(sp == "ach", trt != "ach_low"))
simulateResiduals(ach_2_FD_binom, plot = T)
r2(ach_2_FD_binom)
summary(ach_2_FD_binom)
Anova(ach_2_FD_binom)
emmeans(ach_2_FD_binom, pairwise ~ trt|round)

ach_2_FD_binom <- glmmTMB(cbind(alive_n, density - alive_n) ~ trt * burn * as.factor(round) + (1 | block/cage),
                          family = binomial,
                          data = cage_exp_allrounds %>%
                            filter(sp == "ach", trt %in% c("ach_33", "ach_66", "ach_high")))

# Grass dominance graph on overall species survival, no treatmets, simple graph 


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
  #geom_smooth(method = "glm", method.args = list(family = "betareg")) +
  geom_point(aes(color = burn)) +
  facet_grid( ~ sp) + 
  theme_bw(base_size = 20) + 
  labs(x = "Grass Percentage", 
       y = "Survival Proportion", 
       title = "Survival across grass abundance") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

#basic overall fd on both species graph 

ggplot(cage_exp_fd %>% 
         filter(round == 5), aes(x = trt, y = perc)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(burn ~ sp) + 
  theme_bw(base_size = 20) + 
  labs(x = "Frequency Treatment", 
       y = "Survival Proportion", 
       title = "Grasshoppers mixture survival in burned vs. unburned") + 
  theme(plot.title = element_text(hjust = 0.4, 
                                  face = "bold", 
                                  size = 22))

# model I think we had for FD, discussed with hahn to use, but no longer runs 

grassdominance2_apt_fd <- glmmTMB(perc_survival ~ grass_perc * trt * burn,
                                  data = surv_plant %>% 
                                    filter(sp == "apt", trt != "apt_low", round == 5), family = "ordbeta")

plot(simulateResiduals(grassdominance2_apt_fd))

summary(grassdominance2_apt_fd)
Anova(grassdominance2_apt_fd)
emmeans(grassdominance2_apt_fd,pairwise ~ trt, type = "response")


####################################################
######### NEW PART 2 DD MODELS ##################
##################################################

# using starting density 

ach_2_DD_ignore <- glmmTMB(cbind(alive_n, density - alive_n) ~ density * burn * as.factor(round) + (1 | block/cage),
                    family = binomial,
                    data = cages_2 %>%
                      filter(dep == "monoculture", sp == "ach"))

simulateResiduals(ach_2_DD_ignore, plot = T)
check_autocorrelation(ach_2_DD_ignore)
summary(ach_2_DD_ignore)
Anova(ach_2_DD_ignore)
emtrends(ach_2_DD_ignore,pairwise ~ burn|round, var = "density", infer = T)

# Proportion survival of A. carinatum  

ggplot(cages_2 %>% 
         filter(dep == "monoculture", sp == "ach"),
       aes(x = burn, y = perc_survival, color = high_low)) +
  geom_boxplot() + 
  geom_point(size = 2.8,
             alpha = 0.7) +
  facet_grid(round ~ sp) +
  scale_color_viridis_d(option = "magma", begin = 0.5, end = 0.75,
                        labels = c("u" = "Unburned",
                                   "b" = "Burned")) +
  theme_bw(base_size = 20) +
  labs(x = "Burn treatment",
       y = "Survival Proportion",
       color = "Burn treatment") +
  theme(plot.title = element_text(hjust = 0.4,
                                  face = "bold",
                                  size = 22),
        strip.background = element_rect(fill = "grey95",
                                        color = "black"),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

# Proportion survival of A. sphenarioides  

ggplot(cages_2 %>% 
         filter(dep == "monoculture", sp == "apt"),
       aes(x = burn, y = perc_survival, color = high_low)) +
  geom_boxplot() + 
  geom_point(size = 2.8,
             alpha = 0.7) +
  facet_grid(round ~ sp) +
  scale_color_viridis_d(option = "magma", begin = 0.5, end = 0.75,
                        labels = c("u" = "Unburned",
                                   "b" = "Burned")) +
  theme_bw(base_size = 20) +
  labs(x = "Burn treatment",
       y = "Survival Proportion",
       color = "Burn treatment") +
  theme(plot.title = element_text(hjust = 0.4,
                                  face = "bold",
                                  size = 22),
        strip.background = element_rect(fill = "grey95",
                                        color = "black"),
        strip.text = element_text(face = "bold"),
        legend.position = "top")


###########################

# using high_low categories

ach_2_DD_hl <- glmmTMB(cbind(alive_n, density - alive_n) ~ high_low * burn * as.factor(round) + (1 | block/cage),
                       family = binomial,
                       data = cages_2 %>%
                         filter(dep == "monoculture", sp == "ach"))

simulateResiduals(ach_2_DD_hl, plot = T)
summary(ach_2_DD_hl)
Anova(ach_2_DD_hl)
emmeans(ach_2_DD_hl,pairwise ~ high_low|round:burn, type = "response") 

# ACH DD GLM model on just round 5 

ach_2_DD_5 <- glmmTMB(cbind(alive_n, density - alive_n) ~ density * burn + (1 | block/cage),
                      family = binomial,
                      data = cage_exp_allrounds %>%
                        filter(dep == "monoculture", sp == "ach", round == 5))

simulateResiduals(ach_2_DD_5, plot = T)
summary(ach_2_DD_5)
Anova(ach_2_DD_5)
emtrends(ach_2_DD_5,pairwise ~ burn, var = "density", infer = T)

# APT DD GLM model on just round 5 

apt_2_DD_5 <- glmmTMB(cbind(alive_n, density - alive_n) ~ density * burn + (1 | block/cage),
                      family = binomial,
                      data = cage_exp_allrounds %>%
                        filter(dep == "monoculture", sp == "apt", round == 5))

simulateResiduals(apt_2_DD_5, plot = T)
summary(apt_2_DD_5)
Anova(apt_2_DD_5)
emtrends(apt_2_DD_5,pairwise ~ burn, var = "density", infer = T)


# APT FD with grass percentage (crazy effects)

apt_2_fd <- glmmTMB(perc_survival ~ grass_perc * trt * burn,
                    data = cages_2 %>% 
                      filter(sp == "apt", trt != "apt_low", round == 5, trt != "control"), family = "ordbeta")

plot(simulateResiduals(apt_2_fd))

summary(apt_2_fd)
Anova(apt_2_fd)
emmeans(apt_2_fd,pairwise ~ trt|burn, type = "response")

#### ADDITIONAL: SPIDER PREDATION DATA ####

cage_exp_spider <- cage_exp %>% 
  mutate(spider_present = if_else(!is.na(spider) & spider == 1, 1, 0)) %>% 
  group_by(cage, round) %>% 
  summarise(
    spider_present = max(spider_present, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    spider_present = factor(spider_present, levels = c(0, 1), labels = c("no", "yes"))
  ) %>% 
  filter(cage <= 60) %>% 
  select(c(spider_present, cage, round))

cage_exp_spider <- cage_exp_surv %>% 
  left_join(cage_exp_spider, by = c("cage","round")) %>% 
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
  group_by(burn) %>% 
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


