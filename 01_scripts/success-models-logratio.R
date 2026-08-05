#### Source Success Models ####

## Description: exploring source population success from ancestry data
## Author: Lucas Eckert (lucas.eckert@mail.mcgill.ca)

## packages 
library(tidyverse)
library(lme4)
library(lmerTest)
library(effectsize)
library(emmeans)
library(purrr)
library(cowplot)
library(flextable)

## importing data
data<-left_join(
  read.csv("02_tidy-data/success-data-logratio.csv"),
  select(read.csv("02_tidy-data/env-dists-new.csv"),
         lake=recipient_lake, source=source_lake, dist_1, dist_2, dist_3=dist_4))

## scale data
dist_vars<-c(paste0("dist_",1:3))
data_scaled <- data %>%
  group_by(time) %>%                     
  mutate(across(all_of(dist_vars),     
                ~ scale(.)[,1])) %>% 
  ungroup() %>% 
  group_by(time,lake) %>% 
  mutate(log_init_prop=scale(log_init_prop)[,1]) %>% 
  ungroup()

## models
models<-list()
times<-c("logratio_1","logratio_2","logratio_net")
dists <-c("rlacalc","preypre","allpre")
for (i in 1:length(times)){
  time_label<-str_remove(times[[i]],"logratio_")
  for(j in 1:length(dists)){
    dist_col<-paste0("dist_",j)
    formula_str <- paste0("logratio ~ source + ", dist_col, " + log_init_prop + (1 | lake)")
    model<-lmer(as.formula(formula_str), data=filter(data_scaled, time == time_label, !is.na(logratio)))
    model_label<-paste("model",time_label,dists[j],sep="_")
    n<-((i-1)*3)+j
    models[[n]]<-model
    names(models)[n]<-model_label
  }
}

## results
model_summary_df <- map_df(
  names(models),
  function(mname) {
    model <- models[[mname]]
    
    # Standardized coefficients (without refitting)
    std_df <- tryCatch(
      as.data.frame(standardize_parameters(model, method = "basic", ci = 0.95)),
      error = function(e) {
        message("Skipping ", mname, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(std_df)) return(NULL)
    
    # Normalize column names
    std_df <- std_df %>%
      rename_with(~str_replace_all(., "\\s+", "_")) %>%
      rename_with(~str_replace_all(., "\\.", "_"))
    
    coef_col <- intersect(names(std_df), c("Std_Coefficient", "Std_Estimate", "Coefficient"))
    ci_low_col <- intersect(names(std_df), c("CI_low", "CI_Lower", "CI_low_95"))
    ci_high_col <- intersect(names(std_df), c("CI_high", "CI_Upper", "CI_high_95"))
    
    std_df <- std_df %>%
      rename(
        Std_Coefficient = all_of(coef_col),
        CI_low = all_of(ci_low_col),
        CI_high = all_of(ci_high_col)
      ) %>%
      mutate(
        Std_SE = (CI_high - CI_low) / (2 * 1.96)
      )
    
    # Extract coefficients with t, df, p-value using lmerTest
    coef_df <- as.data.frame(summary(model)$coefficients) %>%
      rownames_to_column("Parameter_name") %>%
      rename(
        Estimate = Estimate,
        Std_Error = `Std. Error`,
        t_value = `t value`,
        df = df,
        p_value = `Pr(>|t|)`
      )
    
    # Merge standardized effects with coefficient table
    merged <- std_df %>%
      select(Parameter, Std_Coefficient, CI_low, CI_high) %>%
      left_join(
        coef_df,
        by = c("Parameter" = "Parameter_name")
      ) %>%
      mutate(model_name = mname)
    
    merged
  }
)

## cleaning up df
model_summary_full<-model_summary_df %>%
  mutate(timespan = map_chr(str_split(model_name, "_"), 2),
         env_dist = map_chr(str_split(model_name, "_"), 3)) %>% 
  select(model=model_name, timespan, env_dist, parameter=Parameter, estimate_std=Std_Coefficient, se=Std_Error,
         ci_low=CI_low, ci_high=CI_high, df,t_value, p_value) %>% 
  mutate(parameter = str_replace(parameter, "^dist_[0-9]+$", "dist")) %>% 
  group_by(timespan) %>% 
  mutate(p_adj=p.adjust(p_value, method = "BH")) %>% 
  ungroup()

## essential summary
model_summary<-model_summary_full %>% 
  filter(parameter %in% c("dist","log_init_prop")) %>% 
  separate(model, into = c("model_prefix","timespan","env_dist"), sep = "_") %>%
  select(-model_prefix) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  mutate(ci=paste0(ci_low," - ",ci_high)) %>% 
  select(timespan, env_dist, parameter, estimate=estimate_std, se, ci, df, t_value, p_value, p_adj)
view(model_summary)

model_summary %>% 
  flextable() %>% 
  save_as_docx(path="03_figures/tables/full-model-summary.docx")

## source summary
emmeans_summary_df <- map_df(
  names(models),
  function(mname) {
    model <- models[[mname]]
    
    ## emmeans for source levels
    em <- emmeans(model, ~ source)
    
    ## df
    em_df <- as.data.frame(em) %>%
      rename(
        source_level = source,
        estimate = emmean,
        se = SE,
        df = df,
        ci_lower = lower.CL,
        ci_upper = upper.CL
      ) %>%
      mutate(model_name = mname)
    em_df
  })


## summary
source_summary<-emmeans_summary_df %>%
  separate(model_name, into = c("model_prefix","timespan","env_dist"), sep = "_") %>%
  select(-model_prefix) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  mutate(ci=paste0(ci_lower," - ",ci_upper)) %>%
  select(timespan, env_dist, source=source_level, estimate, se, df, ci)

##epxort
source_summary %>% 
  arrange(timespan,source) %>% 
  flextable() %>% 
  save_as_docx(path="03_figures/tables/source-emm-summary.docx")

## source means
source_means<-source_summary %>% 
  group_by(source, timespan) %>% 
  summarize(mean_effect=mean(estimate)) %>% 
  filter(timespan!="2") %>% 
  pivot_wider(names_from = timespan, values_from = mean_effect,
              names_prefix = "mean_effect_")

## saving
write.csv(source_means, file="02_tidy-data/source-effect-means-new.csv",
          quote = F, row.names = F)

## source p values

source_effects <- map_df(names(models), function(mname) {
  
  model <- models[[mname]]
  
  # Type III ANOVA using lmerTest
  a <- tryCatch(
    anova(model, type = 3),
    error = function(e) {
      message("Failed ANOVA for ", mname, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(a)) return(NULL)
  
  # Extract 'source' row
  if (!"source" %in% rownames(a)) {
    message("'source' not found in ANOVA for ", mname)
    return(NULL)
  }
  
  source_row <- as.data.frame(a["source", ]) %>%
    rownames_to_column("Effect") %>%
    mutate(model_name = mname) %>%
    select(model_name, Effect, everything())
  
  source_row
})

## clean
source_effects<-source_effects %>% 
  separate(model_name, into = c("model_prefix","timespan","env_dist"), sep = "_") %>%
  select(-model_prefix) %>%
  select(timespan, env_dist, sum_sq='Sum Sq', mean_sq='Mean Sq', df='NumDF', den_df=DenDF, f_value='F value', p_value='Pr(>F)') %>% 
  group_by(timespan) %>%
  mutate(p_adj=p.adjust(p_value, method = "BH")) %>% 
  ungroup()

source_effects %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>% 
  flextable() %>% 
  save_as_docx(path="03_figures/tables/source-effects-table.docx")
