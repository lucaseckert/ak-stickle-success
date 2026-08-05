#### Source Success Models ####

## Description: exploring source population success from ancestry data
## Author: Lucas Eckert (lucas.eckert@mail.mcgill.ca)

## packages 
library(tidyverse)
library(lme4)
library(nlme)
library(lmerTest)
library(emmeans)
library(effectsize)
library(purrr)
library(broom.mixed)
library(car)

## importing data
data<-left_join(
  read.csv("02_tidy-data/success-data-logratio.csv"),
  rename(read.csv("02_tidy-data/env-dists-new.csv"),
         lake=recipient_lake, source=source_lake))

## scale data
dist_vars<-c(paste0("dist_",1:5))
data_scaled <- data %>%
  group_by(time) %>%                     
  mutate(across(all_of(dist_vars),     
                ~ scale(.)[,1])) %>% 
  ungroup() %>% 
  group_by(time,lake) %>% 
  mutate(log_init_prop=scale(log_init_prop)[,1]) %>% 
  ungroup()
  
## models
models_ml<-list()
models_reml<-list()
times <- c("logratio_1","logratio_2","logratio_net")
dists <- c("rlacalc","preypre","preypost","allpre","allpost")

## loop for ML and REML models
n <- 0
for (i in seq_along(times)) {
  time_label <- stringr::str_remove(times[[i]], "logratio_") 
  
  for (j in seq_along(dists)) {
    n <- n + 1
    dist_col <- paste0("dist_", j)
    model_label <- paste("model", time_label, dists[j], sep = "_")
    formula_str <- paste0("logratio ~ source + ", dist_col, " + log_init_prop")
    
    ## subset data
    current_data <- filter(data_scaled, time == time_label, !is.na(logratio))
    
    ## ML fit attempt
    ml_fit <- try(
      lme(
        fixed = as.formula(formula_str),
        random = ~1 | lake,
        correlation = corCompSymm(form = ~1 | lake),
        data = current_data,
        method = "ML", 
        control = lmeControl(maxIter = 1000, opt = "optim")
      ),
      silent = TRUE
    )
    
    ## store ML results 
    if (inherits(ml_fit, "try-error")) {
      message("ML failed for ", model_label, ". Storing NULL.")
      models_ml[[n]] <- NULL 
    } else {
      models_ml[[n]] <- ml_fit
      names(models_ml)[n] <- model_label
      message("Successfully fitted and stored ML model: ", model_label)
    }
    
    ## REML fit attempt
    reml_fit <- try(
      lme(
        fixed = as.formula(formula_str),
        random = ~1 | lake,
        correlation = corCompSymm(form = ~1 | lake),
        data = current_data,
        method = "REML",
        control = lmeControl(maxIter = 1000, opt = "optim")
      ),
      silent = TRUE
    )
    
    ## store REML results
    if (inherits(reml_fit, "try-error")) {
      message("REML failed for ", model_label, ". Storing NULL.")
      models_reml[[n]] <- NULL
    } else {
      models_reml[[n]] <- reml_fit
      names(models_reml)[n] <- model_label
      message("Successfully fitted and stored REML model: ", model_label)
    }
  }
}

## summarize effect sizes
expected_terms<-c("log_init_prop",dist_vars)

## function to get results
get_results<-function(models_list) {
  
  ## expected terms
  expected_terms<-c("log_init_prop",dist_vars)
  
  ## get results
  results <- map_dfr(names(models_list), function(mn) {
    mod <- models_list[[mn]]
    
    ## model failed
    if (is.null(mod) || inherits(mod, "try-error")) {
      return(tibble(
        model = mn,
        Parameter = expected_terms,
        Std_Estimate = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        p_value = NA_real_
      ))
    }
    
    ## model exists
    broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE) %>%
      filter(term %in% expected_terms) %>%
      mutate(model = mn) %>%
      select(
        model,
        Parameter = term,
        Std_Estimate = estimate,
        Std_Error = std.error,   # New: Standard Error
        df = df,                 # New: Degrees of Freedom
        t_value = statistic,     # New: The t-statistic
        CI_low = conf.low,
        CI_high = conf.high,
        p_value = p.value
      )
  }) %>% 
    ## filtering 
    filter(!is.na(Std_Estimate)) %>%
    mutate(
      timespan = map_chr(str_split(model, "_"), 2),
      env_dist = map_chr(str_split(model, "_"), 3)
    )
  
  ## return results
  return(results)
}

## get REML results
reml_results<-get_results(models_reml) %>% 
  rename(fixed_effect=Parameter) %>% 
  mutate(fixed_effect = str_replace(fixed_effect, "^dist_[0-9]+$", "env_dist"))

## correct for multiple testing
reml_results <- reml_results %>%
  group_by(fixed_effect) %>%
  mutate(p_value_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()
view(reml_results)

## subset
reml_subset<-reml_results %>% 
  filter(env_dist %in% c("rlacalc", "preypost", "allpost")) %>% 
  group_by(timespan, fixed_effect) %>% 
  mutate(p_value_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()
view(reml_subset)

## get ML results
ml_results<-get_results(models_ml) %>% 
  rename(fixed_effect=Parameter) %>% 
  mutate(fixed_effect = str_replace(fixed_effect, "^dist_[0-9]+$", "env_dist"))

## correct for multiple testing
ml_results <- ml_results %>%
  group_by(fixed_effect) %>%
  mutate(p_value_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()
view(ml_results)

## subset
ml_subset<-ml_results %>% 
  filter(env_dist %in% c("rlacalc", "preypre", "allpre")) %>% 
  group_by(timespan, fixed_effect) %>% 
  mutate(p_value_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()
view(ml_subset)

## copy ml_subset as table that I can paste into word
ml_subset %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  mutate(ci=paste0(CI_low, " - ", CI_high)) %>%
  select(time=timespan,
         env=env_dist,
         parameter=fixed_effect,
         estimate_std=Std_Estimate,
         std_error=Std_Error,
         ci,
         t_value,
         df,
         p_value,
         p_value_adj) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>% 
  flextable::flextable() %>% 
  flextable::save_as_docx(path="03_figures/tables/ml-subset-results.docx")

## LRT for source effect
get_source_LRT <- function(model) {
  if (is.null(model)) return(data.frame(LR = NA, df = NA, p_val = NA, delta_AIC = NA))
  
  # 1. Fit reduced model
  df_orig <- model$data
  reduced <- update(model, . ~ . - source, data = df_orig)
  
  # 2. Get ANOVA/LRT results
  lr_comp <- anova(reduced, model)
  
  # 3. Calculate Delta AIC (Reduced - Full)
  # Positive value means the Full model is better (lower AIC)
  d_aic <- AIC(reduced) - AIC(model)
  
  # Return as a row
  data.frame(
    LR = lr_comp$L.Ratio[2],
    df = lr_comp$df[2] - lr_comp$df[1], # Difference in df
    p_val = lr_comp$`p-value`[2],
    delta_AIC = d_aic
  )
}

## LRT for source effect
LRT_results<- map_dfr(models_ml, get_source_LRT, .id = "model_name") %>% 
  separate(model_name, into = c("model_prefix", "timespan", "env_dist"), sep = "_") %>% 
  select(-model_prefix) %>% 
  filter(env_dist %in% c("rlacalc", "preypre", "allpre")) %>% 
  group_by(timespan) %>% 
  mutate(p_value_adj = p.adjust(p_val, method = "BH")) %>% 
  select(time=timespan,env_dist,LR,df,delta_AIC,p_val,p_value_adj)

## export
LRT_results %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  flextable::flextable() %>% 
  flextable::save_as_docx(path="03_figures/tables/source-lrt-results.docx")

  ## get source effects
get_source_emm <- function(model) {
  if (is.null(model)) return(NULL)
  
  # Use emmeans to compute marginal means for 'source'
  emm <- emmeans(model, specs = "source")
  
  # Convert to data frame with estimates and 95% CI
  as.data.frame(emm)
}

## get means
source_emms_list <- lapply(models_ml, get_source_emm)
source_emms_list <- Map(cbind, model = names(models_ml), source_emms_list)

## clean up
source_emms_df <- bind_rows(source_emms_list) %>% 
  separate(model, into = c("model_prefix", "timespan", "env_dist"), sep = "_") %>%
  select(-model_prefix)

## subset
source_emms_subset<- source_emms_df %>% 
  filter(env_dist %in% c("rlacalc", "preypre", "allpre"))

## export
source_emms_subset %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  mutate(ci=paste0(lower.CL, " - ", upper.CL)) %>%
  select(timespan,
         env_dist,
         source,
         emmean,
         SE,
         df,
         ci) %>%
  arrange(timespan, source) %>%
  flextable::flextable() %>% 
  flextable::save_as_docx(path="03_figures/tables/source-emm-results.docx")

## means
source_means<-source_emms_subset %>% 
  group_by(source, timespan) %>% 
  summarize(mean_effect=mean(emmean)) %>% 
  filter(timespan!="2") %>% 
  pivot_wider(names_from = timespan, values_from = mean_effect,
              names_prefix = "mean_effect_")

## saving
write.csv(source_means, file="02_tidy-data/source-effect-means-new.csv",
          quote = F, row.names = F)

