#### SOURCE SUCCESS MODELS ####

## Description: exploring source population success from ancestry data
## Author: Lucas Eckert (lucas.eckert@mail.mcgill.ca)

## packages 
source(here::here("02_tidy-data/constants.R"))
library(here)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(performance) 

## ancestry data
count_data<-read.csv(here("02_tidy-data","source-success-counts.csv"))

## combining with environmental data
data<-left_join(
  count_data,
  select(read.csv(here("02_tidy-data","env-dists-new.csv")),
         lake=recipient_lake, source=source_lake, dist_1, dist_2, dist_3=dist_4)) %>% 
  ## adjusting zero proportions
  mutate(across(starts_with("prop_"), ~ ifelse(.x == 0, .x + 0.01, .x), .names = "{.col}_adj")) %>% 
  ## adjusting counts for the F1 and F2
  mutate(count_t1_adj = count_t1*2,
         count_t2_adj = count_t2*4) %>% 
  ## computing logratios
  mutate(logratio_1 = log(prop_t1_adj/prop_t0_adj),
         logratio_2 = log(prop_t2_adj/prop_t1_adj),
         logratio_net = log(prop_t2_adj/prop_t0_adj))

## glm olre models
data$obs_id <- 1:nrow(data)
glmm_olre_models<-list(
  glmm_1_1 = glmmTMB(count_t1_adj ~ 0 + source + scale(dist_1) + offset(log(prop_t0_adj)) + (1|lake) + (1|obs_id), family=poisson, data=data),
  glmm_1_2 = glmmTMB(count_t1_adj ~ 0 + source + scale(dist_2) + offset(log(prop_t0_adj)) + (1|lake) + (1|obs_id), family=poisson, data=data),
  glmm_1_3 = glmmTMB(count_t1_adj ~ 0 + source + scale(dist_3) + offset(log(prop_t0_adj)) + (1|lake) + (1|obs_id), family=poisson, data=data),
  glmm_2_1 = glmmTMB(count_t2_adj ~ 0 + source + scale(dist_1) + offset(log(prop_t1_adj)) + (1|lake) + (1|obs_id), family=poisson, data=data),
  glmm_2_2 = glmmTMB(count_t2_adj ~ 0 + source + scale(dist_2) + offset(log(prop_t1_adj)) + (1|lake) + (1|obs_id), family=poisson, data=data),
  glmm_2_3 = glmmTMB(count_t2_adj ~ 0 + source + scale(dist_3) + offset(log(prop_t1_adj)) + (1|lake) + (1|obs_id), family=poisson, data=data),
  glmm_3_1 = glmmTMB(count_t2_adj ~ 0 + source + scale(dist_1) + offset(log(prop_t0_adj)) + (1|lake) + (1|obs_id), family=poisson, data=data), 
  glmm_3_2 = glmmTMB(count_t2_adj ~ 0 + source + scale(dist_2) + offset(log(prop_t0_adj)) + (1|lake) + (1|obs_id), family=poisson, data=data),
  glmm_3_3 = glmmTMB(count_t2_adj ~ 0 + source + scale(dist_3) + offset(log(prop_t0_adj)) + (1|lake) + (1|obs_id), family=poisson, data=data)
)

#### DIAGNOSTICS ####

## checks
sapply(glmm_olre_models, function(x) check_convergence(x)) # all converged 
sapply(glmm_olre_models, function(x) check_singularity(x)) # some singular 
lapply(glmm_olre_models, VarCorr) # singular models show collapse of lake level random effect
sapply(glmm_olre_models, function(x) check_zeroinflation(x)) # all good
sapply(glmm_olre_models, performance::check_overdispersion) # over dispersion fixed

## composition checks - even singular models show good agreement between predicted and actual counts
composition_check1 <- imap_dfr(glmm_olre_models[1:3], function(fit, model_name) {
  data %>%
    filter(!is.na(count_t1_adj)) %>%
    mutate(pred = predict(fit, type = "response")) %>%
    group_by(lake) %>%
    summarise(
      pred_total = sum(pred),
      actual_N   = first(N_t1)*2,
      diff       = pred_total - actual_N,
      pct_diff   = 100 * diff / actual_N,
      .groups    = "drop"
    ) %>%
    mutate(model = model_name)
}) %>%
  select(model, lake, pred_total, actual_N, diff, pct_diff)
view(composition_check1) # all within 5%

composition_check2 <- imap_dfr(glmm_olre_models[4:6], function(fit, model_name) {
  data %>%
    filter(!is.na(count_t1_adj)) %>%
    mutate(pred = predict(fit, type = "response")) %>%
    group_by(lake) %>%
    summarise(
      pred_total = sum(pred),
      actual_N   = first(N_t2)*4,
      diff       = pred_total - actual_N,
      pct_diff   = 100 * diff / actual_N,
      .groups    = "drop"
    ) %>%
    mutate(model = model_name)
}) %>%
  select(model, lake, pred_total, actual_N, diff, pct_diff)
view(composition_check2) # all within 2.5%

composition_check3 <- imap_dfr(glmm_olre_models[7:9], function(fit, model_name) {
  data %>%
    mutate(pred = predict(fit, type = "response")) %>%
    group_by(lake) %>%
    summarise(
      pred_total = sum(pred),
      actual_N   = first(N_t2)*4,
      diff       = pred_total - actual_N,
      pct_diff   = 100 * diff / actual_N,
      .groups    = "drop"
    ) %>%
    mutate(model = model_name)
}) %>%
  select(model, lake, pred_total, actual_N, diff, pct_diff)
view(composition_check3) # all within 1%

#### ENV EFFECTS ####

## env dist effects
dist_effects <- map_dfr(glmm_olre_models, function(mod){
  coefs <- tidy(mod, effects = "fixed")
  coefs %>%
    filter(str_detect(term, "^scale\\(dist_\\d\\)$")) %>%
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      ci = paste0("(", round(conf.low, 3), ", ", round(conf.high, 3), ")")
    ) %>%
    select(estimate, std.error, ci, statistic, p.value)
}, .id = "model_name") %>% 
  mutate(
    timespan = case_when(
      str_detect(model_name, "1_") ~ "1",
      str_detect(model_name, "2_") ~ "2",
      str_detect(model_name, "3_") ~ "net"
    ),
    env_dist = substr(model_name,8,8)
  ) %>% 
  select(model_name, timespan, env_dist, everything()) %>% 
  group_by(timespan) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>% 
  ungroup()

## export
dist_effects %>%
  mutate(p.value = if_else(p.value < 0.001, "<0.001", as.character(round(p.value,3))),
         p_adj = if_else(p_adj < 0.001, "<0.001", as.character(round(p_adj,3)))) %>%
  mutate(across(c(estimate, std.error, statistic), round, 3)) %>%
  flextable() %>% 
  save_as_docx(path="03_figures/tables/env-effects.docx")

#### DIAGNOSTICS SUMMARY ####

## model diag df
model_diagnostics <- imap_dfr(glmm_olre_models, function(mod, name) {
  
  # singularity
  singular <- performance::check_singularity(mod)
  
  # R2
  if (!singular) {
    r2 <- performance::r2(mod)
    r2_marginal <- r2$R2_marginal
    r2_conditional <- r2$R2_conditional
  } else {
    r2_marginal <- NA_real_
    r2_conditional <- NA_real_
  }
  
  # overdispersion
  od <- performance::check_overdispersion(mod)
  
  # predicted vs actual counts
  constant <- if (str_detect(name, "1_")) 2 else 4
  pred_data <- if (str_detect(name, "3_")) {
    data
  } else {
    filter(data, !is.na(count_t1_adj))
  }
  mean_pct_diff <- pred_data %>%
    mutate(pred = predict(mod, type = "response")) %>%
    group_by(lake) %>%
    summarise(
      pred_total = sum(pred),
      actual_N   = first(N_t2) * constant,
      pct_diff   = 100 * (pred_total - actual_N) / actual_N,
      .groups = "drop"
    ) %>%
    summarise(mean_pct_diff = mean(pct_diff)) %>%
    pull(mean_pct_diff)
  
  tibble(
    model_name = name,
    AIC = AIC(mod),
    #BIC = BIC(mod),
    logLik = as.numeric(logLik(mod)),
    singular = singular,
    R2_marginal = r2_marginal,
    R2_conditional = r2_conditional,
    dispersion_ratio = od$dispersion_ratio,
    mean_pct_diff = mean_pct_diff
  )
}) %>% mutate(timespan=case_when(
  str_detect(model_name, "1_") ~ "1",
  str_detect(model_name, "2_") ~ "2",
  str_detect(model_name, "3_") ~ "net"),
  env_dist = substr(model_name,8,8)) %>% 
  select(model_name, timespan, env_dist, everything())
view(model_diagnostics)

## export
model_diagnostics %>%
  mutate(across(c(AIC, logLik, R2_marginal, R2_conditional, dispersion_ratio, mean_pct_diff), round, 2)) %>%
  flextable() %>% 
  save_as_docx(path="03_figures/tables/diagnostics.docx")

#### SOURCE EFFECTS ####

## LRTs for source effect
source_LRT <- lapply(glmm_olre_models, function(mod) {
  full <- update(mod, . ~ . - 0 + 1)  # explicitly remove 0, add intercept
  null <- update(full, . ~ . - source) # remove source, keep intercept
  anova(null, full, test = "LRT")
})

lrt_df <- imap_dfr(source_LRT, function(lrt, model_name) {
  tibble(
    model   = model_name,
    logLik_null = lrt$logLik[1],
    logLik_full = lrt$logLik[2],
    chi_sq  = lrt$Chisq[2],
    df      = lrt$Df[2] - lrt$Df[1],
    p_value = lrt$`Pr(>Chisq)`[2]
  )
})

## export
lrt_df %>% 
  mutate(timespan=case_when(
    str_detect(model, "1_") ~ "1",
    str_detect(model, "2_") ~ "2",
    str_detect(model, "3_") ~ "net"),
    env_dist = substr(model,8,8)) %>% 
  group_by(timespan) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         logLik_null = round(logLik_null,1),
         logLik_full = round(logLik_full,1),
         chi_sq = round(chi_sq, 2)) %>% 
  ungroup() %>% 
  select(model, timespan, env_dist, everything()) %>% 
  mutate(p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value,3))),
         p_adj = if_else(p_adj < 0.001, "<0.001", as.character(round(p_adj,3)))) %>% 
  flextable() %>% 
  save_as_docx(path="03_figures/tables/source-lrt.docx")

#### EFFECT SIZE COMPARISON ####

## get source effects
source_effects<-emmeans::emmeans(glmm_olre_models$glmm_3_3, ~source, type="response") %>% 
  as.data.frame()
max_source_effect<-max(source_effects$rate)/min(source_effects$rate)
mid_source_effect<-max(source_effects$rate)/median(source_effects$rate)

## get dist effects
dist_range <- diff(range(data$dist_3))
dist_sd <- sd(data$dist_3)
range_in_sd <- dist_range / dist_sd
beta_dist <- fixef(glmm_olre_models$glmm_3_3)$cond["scale(dist_3)"]
max_dist_effect<-1/exp(beta_dist * range_in_sd)

#### SOURCE EFFECT EXPORT ####

all_source_effects <- map_dfr(glmm_olre_models, function(mod){
  emmeans::emmeans(mod, ~source, type = "response") %>% 
    as.data.frame()
}, .id = "model_name") %>% 
  mutate(
    timespan = case_when(
      str_detect(model_name, "1_") ~ "1",
      str_detect(model_name, "2_") ~ "2",
      str_detect(model_name, "3_") ~ "net"
    ),
    env_dist = substr(model_name,8,8)
  ) %>% 
  mutate(
    conf.low = rate - 1.96 * SE,
    conf.high = rate + 1.96 * SE,
    ci = paste0("(", round(conf.low, 3), ", ", round(conf.high, 3), ")")
  ) %>% 
  select(model_name, timespan, env_dist, everything()) %>% 
  group_by(model_name) %>% 
  mutate(
    rate_sum = sum(rate),
    rate_prop = rate / rate_sum
  ) %>% 
  ungroup()

## export
all_source_effects %>% 
  select(model_name, timespan, env_dist, source, rate, SE, ci, rate_prop) %>%
  arrange(timespan, source, env_dist) %>%
  mutate(across(c(rate, SE, rate_prop), round, 3)) %>%
  flextable() %>% 
  save_as_docx(path="03_figures/tables/source-emmeans.docx")

mean_source_effects<-all_source_effects %>% 
  filter(timespan != "2") %>% 
  group_by(timespan, source) %>% 
  summarize(mean_rate=mean(rate_prop), .groups="drop")

## save
write.csv(mean_source_effects, "02_tidy-data/mean-source-effects.csv", row.names = FALSE, quote = FALSE)


#### WITHOUT OLRE ####

## glmms without olre
glmm_no_olre<-list(
  glmm_1_1 = glmmTMB(count_t1_adj ~ 0 + source + scale(dist_1) + offset(log(prop_t0_adj)) + (1|lake), family=poisson, data=data),
  glmm_1_2 = glmmTMB(count_t1_adj ~ 0 + source + scale(dist_2) + offset(log(prop_t0_adj)) + (1|lake), family=poisson, data=data),
  glmm_1_3 = glmmTMB(count_t1_adj ~ 0 + source + scale(dist_3) + offset(log(prop_t0_adj)) + (1|lake), family=poisson, data=data),
  glmm_2_1 = glmmTMB(count_t2_adj ~ 0 + source + scale(dist_1) + offset(log(prop_t1_adj)) + (1|lake), family=poisson, data=data),
  glmm_2_2 = glmmTMB(count_t2_adj ~ 0 + source + scale(dist_2) + offset(log(prop_t1_adj)) + (1|lake), family=poisson, data=data),
  glmm_2_3 = glmmTMB(count_t2_adj ~ 0 + source + scale(dist_3) + offset(log(prop_t1_adj)) + (1|lake), family=poisson, data=data),
  glmm_3_1 = glmmTMB(count_t2_adj ~ 0 + source + scale(dist_1) + offset(log(prop_t0_adj)) + (1|lake), family=poisson, data=data), 
  glmm_3_2 = glmmTMB(count_t2_adj ~ 0 + source + scale(dist_2) + offset(log(prop_t0_adj)) + (1|lake), family=poisson, data=data),
  glmm_3_3 = glmmTMB(count_t2_adj ~ 0 + source + scale(dist_3) + offset(log(prop_t0_adj)) + (1|lake), family=poisson, data=data)
)

## disperson ratio
sapply(glmm_no_olre, function(x) performance::check_overdispersion(x))

## env effects
dist_effects_no_olre <- map_dfr(glmm_no_olre, function(mod){
  coefs <- tidy(mod, effects = "fixed")
  coefs %>%
    filter(str_detect(term, "^scale\\(dist_\\d\\)$")) %>%
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      ci = paste0("(", round(conf.low, 3), ", ", round(conf.high, 3), ")")
    ) %>%
    select(estimate, std.error, ci, statistic, p.value)
}, .id = "model_name") %>% 
  mutate(
    timespan = case_when(
      str_detect(model_name, "1_") ~ "1",
      str_detect(model_name, "2_") ~ "2",
      str_detect(model_name, "3_") ~ "net"
    ),
    env_dist = substr(model_name,8,8)
  ) %>% 
  select(model_name, timespan, env_dist, everything()) %>% 
  group_by(timespan) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>% 
  ungroup()

## export
dist_effects_no_olre %>%
  select(model_name, timespan, env_dist, estimate, std.error, p.value) %>% 
  mutate(across(c(estimate, std.error, p.value), round, 3)) %>%
  flextable() %>% 
  save_as_docx(path="03_figures/tables/env-effects-no-olre.docx")


