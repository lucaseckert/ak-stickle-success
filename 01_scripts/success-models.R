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
  read.csv("02_tidy-data/success-data.csv"),
  rename(read.csv("02_tidy-data/env-dists.csv"),
         lake=recipient_lake, source=source_lake))

## models
models<-list()
times<-c("delta_1","delta_2","delta_net")
dists<-c("rla","calc","both","prey","all")
for (i in 1:length(times)){
  time_label<-str_remove(times[[i]],"delta_")
  for(j in 1:length(dists)){
    dist_col<-paste0("dist_",j)
    formula_str <- paste0("delta_norm ~ source + ", dist_col, " + initial_prop_norm + (1 | lake)")
    model<-lmer(as.formula(formula_str), data=filter(data, time == times[i]))
    model_label<-paste("model",time_label,dists[j],sep="_")
    n<-((i-1)*5)+j
    models[[n]]<-model
    names(models)[n]<-model_label
  }
}

## results
model_summary_df <- map_df(
  names(models),
  function(mname) {
    model <- models[[mname]]
    
    ## standardized coefficients without refitting
    std_df <- tryCatch(
      as.data.frame(standardize_parameters(model, method = "basic", ci = 0.95)),
      error = function(e) {
        message("Skipping ", mname, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(std_df)) return(NULL)
    
    ## normalize column names
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
      )
    
    ## extract p-values
    coef_df <- as.data.frame(summary(model)$coefficients)
    
    ## merge standardized effects and p-values
    merged <- std_df %>%
      select(Parameter, Std_Coefficient, CI_low, CI_high) %>%
      left_join(
        coef_df %>%
          rownames_to_column("Parameter") %>%
          select(Parameter, Estimate, `Pr(>|t|)`),
        by = "Parameter"
      ) %>%
      mutate(model_name = mname)
    
    merged
  }
)


## cleaning up df
model_summary_full<-model_summary_df %>%
  mutate(timespan = map_chr(str_split(model_name, "_"), 2),
         env_dist = map_chr(str_split(model_name, "_"), 3)) %>% 
  select(model=model_name, timespan, env_dist, parameter=Parameter, estimate=Estimate, estimate_std=Std_Coefficient,
         ci_low=CI_low, ci_high=CI_high, p_value='Pr(>|t|)') %>% 
  mutate(parameter = str_replace(parameter, "^dist_[0-9]+$", "dist"))

## essential summary
model_summary<-model_summary_full %>% 
  filter(parameter %in% c("dist","initial_prop_norm"))
view(model_summary)

## export
model_summary %>% mutate(across(where(is.numeric), ~ round(., 4))) %>% 
  select(-estimate) %>% 
  pivot_wider(names_from = parameter, values_from = c(estimate_std, ci_high, ci_low, p_value)) %>% 
  mutate(b_env=paste0(estimate_std_dist," (",ci_low_dist," \u2013 ",ci_high_dist,")"),
         b_freq=paste0(estimate_std_initial_prop_norm," (",ci_low_initial_prop_norm," \u2013 ",ci_high_initial_prop_norm,")")) %>% 
  select(time=timespan,
         env=env_dist,
         b_env,
         p_env=p_value_dist,
         b_freq,
         p_freq=p_value_initial_prop_norm) %>%
  flextable() %>% 
  save_as_docx(path="03_figures/tables/model-table.docx")

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

source_summary<-emmeans_summary_df %>%
  mutate(timespan = map_chr(str_split(model_name, "_"), 2),
         env_dist = map_chr(str_split(model_name, "_"), 3)) %>% 
  select(model_name, timespan, env_dist, source_level, estimate, se, df, ci_lower, ci_upper)
view(source_summary)

## net means
source_means_net<-source_summary %>% filter(timespan=="net") %>% 
  group_by(source_level) %>% 
  summarize(mean_effect_net=mean(estimate))
source_means<-source_summary %>% filter(timespan=="1") %>% 
  group_by(source_level) %>% 
  summarize(mean_effect_1=mean(estimate)) %>% 
  left_join(source_means_net)

## saving
write.csv(source_means, file="02_tidy-data/source-effect-means.csv",
          quote = F, row.names = F)

## source p values
source_effects <- map_df(names(models), function(mname) {
  model <- models[[mname]]
  a <- anova(model, type = 3)
  
  # Extract the row for 'source'
  source_row <- as.data.frame(a["source", ]) %>%
    mutate(model_name = mname) %>%
    select(model_name, everything())
  source_row
})

row.names(source_effects)<-NULL
source_effects %>% 
  mutate(across(where(is.numeric), ~ round(., 4))) %>% 
  flextable() %>% 
  save_as_docx(path="03_figures/tables/source-effects-table.docx")


source_order<-arrange(source_means, mean_effect_net)$source_level
source_summary <- source_summary %>%
  mutate(source_level = factor(source_level, levels = source_order))
source_means <- source_means %>%
  mutate(source_level = factor(source_level, levels = source_order))

source_cols_alt_named <- setNames(source_cols_alt, sources)
source_cols_alt_order <- source_cols_alt_named[as.character(source_order)]

source_effect_plot<-source_summary %>% filter(timespan=="net") %>% 
ggplot(aes(x=source_level, y=estimate, fill=source_level))+
  geom_point(shape=21, size=4, alpha=0.75)+
  geom_point(data = source_means, aes(x=source_level, fill=source_level, y=mean_effect_net),
             size=8, shape=3)+
  scale_fill_manual(values=source_cols_alt_order)+
  labs(x="Source Population", y="Effect Size")+
  annotate("text", label="A", x=1.15, y=0.22, size=12)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))

freq_dep_plot<-readRDS("03_figures/figure-parts/freq-dep-plot.RDS")

## saving
fig<- plot_grid(source_effect_plot,freq_dep_plot, ncol = 2)
ggsave(fig, filename = "03_figures/model-effects-fig.tiff",
       device = "tiff",
       width = 1440, height=600, units = "px",
       dpi=300, scale=2.5)


