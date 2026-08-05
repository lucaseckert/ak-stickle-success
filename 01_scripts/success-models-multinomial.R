#### Source Success Models ####

## Description: exploring source population success from ancestry data
## Author: Lucas Eckert (lucas.eckert@mail.mcgill.ca)

## packages 
library(tidyverse)
library(brms)

## importing data
data<-left_join(
  read.csv("02_tidy-data/success-data-multinomial.csv"),
  rename(read.csv("02_tidy-data/env-dists.csv"),
         lake=recipient_lake, source=source_lake)) %>% 
  pivot_longer(cols=starts_with("dist"),
               values_to = "env_dist",
               names_to = "env",
               names_prefix = "dist_") %>% 
  mutate(count=if_else(count==0 & introduced==1,1,count),
         n_total=1000)

data_wide<-data %>% 
  pivot_wider(names_from = source, values_from = c(count, init_prop, env_dist, introduced))




fit <- brm(
  bf(count | trials(n_total) ~ 0 + source +
       source:env_dist + source:init_prop +
       source:introduced + (1 | lake)),
  family = dirichlet(),
  data = filter(data,time=="net",env=="1"),
  chains = 4, cores = 4
)

# The model formula structure for brms
model_formula <- brms::bf(
  count | trials(n_total) ~ source + init_prop + env_dist + (1 | lake)
)

# Fit the model using the Binomial family
model_binom <- brms::brm(
  formula = model_formula,
  data = filter(data,time=="net",env=="1",introduced==1),
  family = binomial(),
  chains = 4, cores = 4
)


# One formula per source population
f_SR <- bf(count_SR | trials(n_total) ~ 1 + env_dist_SR + init_prop_SR + introduced_SR + (1 | lake))
f_WK <- bf(count_WK | trials(n_total) ~ 1 + env_dist_WK + init_prop_WK + introduced_WK + (1 | lake))
f_LG <- bf(count_LG | trials(n_total) ~ 1 + env_dist_LG + init_prop_LG + introduced_LG + (1 | lake))
f_SL <- bf(count_SL | trials(n_total) ~ 1 + env_dist_SL + init_prop_SL + introduced_SL + (1 | lake))
f_WB <- bf(count_WB | trials(n_total) ~ 1 + env_dist_WB + init_prop_WB + introduced_WB + (1 | lake))
f_FG <- bf(count_FG | trials(n_total) ~ 1 + env_dist_FG + init_prop_FG + introduced_FG + (1 | lake))
f_WT <- bf(count_WT | trials(n_total) ~ 1 + env_dist_WT + init_prop_WT + introduced_WT + (1 | lake))
f_TL <- bf(count_TL | trials(n_total) ~ 1 + env_dist_TL + init_prop_TL + introduced_TL + (1 | lake))

fit_multinom <- brm(
  mvbind(count_SR, count_WK, count_LG, count_SL,
         count_WB, count_FG, count_WT, count_TL) | trials(n_total) ~
    1 +
    env_dist_SR + env_dist_WK + env_dist_LG + env_dist_SL +
    env_dist_WB + env_dist_FG + env_dist_WT + env_dist_TL +
    init_prop_SR + init_prop_WK + init_prop_LG + init_prop_SL +
    init_prop_WB + init_prop_FG + init_prop_WT + init_prop_TL +
    introduced_SR + introduced_WK + introduced_LG + introduced_SL +
    introduced_WB + introduced_FG + introduced_WT + introduced_TL +
    (1 | lake),
  family = multinomial(),
  data   = filter(data_wide, time == "net", env == "1"),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95)
)

dat <- filter(data_wide, time == "net", env == "1")
str(dat)
summary(dat[, grep("^count_", names(dat))])


dat <- filter(data_wide, time == "net", env == "1")
fit_test <- brm(
  mvbind(count_SR, count_WK, count_LG, count_SL,
         count_WB, count_FG, count_WT, count_TL) | trials(n_total) ~ 1,
  family = multinomial(),
  data = dat,
  chains = 2, cores = 2, iter = 1000
)
