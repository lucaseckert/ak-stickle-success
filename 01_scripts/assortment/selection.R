#### SELECTION ####

## packages
library(tidyverse)
library(broom)
library(ggh4x)

## constants
source(here::here("02_tidy-data/constants.R"))

## data
data<-read.csv("02_tidy-data/assortment/final-f1-data.csv") %>% 
  filter(!is.na(std_length), !is.na(FG_round)) %>% 
  select(fish_id, lake, year, std_length, ends_with("round"))

## adding mass
data<-data %>% 
  left_join(select(read.csv("02_tidy-data/other-data/lethal-data.csv"), fish_id=fishID, mass))

## hybrids
data<-data %>% 
  mutate(hybrid=if_any(ends_with("round"), ~ .x == 0.5),
         hybrid=factor(hybrid, levels=c(T,F), labels=c("Admixed","Pure")),
         lake=factor(lake, levels=c(limnetic_lakes, benthic_lakes, both_lakes)),
         year=factor(year, levels=c(1,2), labels=c("Age 1","Age 2")),
         k=100000*mass/std_length^3)

means<-data %>% group_by(lake, year, hybrid) %>% 
  summarize(mean_length=mean(std_length, na.rm = T),
            mean_k=mean(k, na.rm = T))

## plotting
source_plot<-ggplot(data, aes(x=hybrid, y=std_length, fill=hybrid))+
  geom_jitter(shape=21, size=2, alpha=0.75, width=0.15)+
  geom_point(
    data = means,
    aes(x = hybrid, y = mean_length),
    shape = 3, size = 3, stroke = 1,
    show.legend = FALSE) +
  geom_line(
    data = means,
    aes(x = hybrid, y = mean_length, group = interaction(lake, year)),
    linewidth = 0.5,
    color = "black",
    linetype = "dashed") +
  scale_fill_manual(values = c("grey55","grey90"), labels=geno_labs)+
  labs(y="Standard Length (mm)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=11),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size=10, face="bold"),
        strip.background = element_rect(fill="grey95"))+
  guides(
    fill = guide_legend(
      override.aes = list(
        size = 3.5,
        alpha = 1,
        shape = 21)))+
  facet_grid2(rows=vars(year),
             cols=vars(lake),
             strip = strip_themed(
               background_x = elem_list_rect(fill = all_lake_cols[-1])))

data %>% filter(LG_round>0, lake=="RL") %>% 
ggplot(aes(x=hybrid, y=std_length, fill=hybrid))+
  geom_jitter(shape=21, size=2, alpha=0.75, width=0.15)+
  #geom_point(
    #data = means,
    #aes(x = hybrid, y = mean_length),
    #shape = 3, size = 3, stroke = 1,
    #show.legend = FALSE) +
  #geom_line(
    #data = means,
    #aes(x = hybrid, y = mean_length, group = interaction(lake, year)),
    #linewidth = 0.5,
    #color = "black",
    #linetype = "dashed") +
  scale_fill_manual(values = c("grey55","grey90"), labels=c("Admixed","Pure"))+
  labs(y="Standard Length (mm)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=11),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size=10, face="bold"),
        strip.background = element_rect(fill="grey95"))+
  guides(
    fill = guide_legend(
      override.aes = list(
        size = 3.5,
        alpha = 1,
        shape = 21)))+
  facet_grid2(rows=vars(year),
              cols=vars(lake),
              strip = strip_themed(
                background_x = elem_list_rect(fill = all_lake_cols[-1])))

ggsave(filename = "03_figures/assortment/lg-selection.tiff",
       width = 500, height=400, units = "px",
       dpi=300, scale=3)

## eco hybrids
eco_data<-data %>% 
  filter(lake %in% c("GL","LO")) %>% 
  mutate(ecotype=case_when(
    FG_round + TL_round + WB_round + WT_round == 1 ~ "benthic",
    LG_round + SL_round + SR_round + WK_round == 1 ~ "limnetic",
    TRUE ~ "hybrid"
  ))

## eco means
eco_means<-eco_data %>% group_by(lake, year, ecotype) %>% 
  summarize(mean_length=mean(std_length, na.rm = T))

## plotting
eco_plot<-ggplot(eco_data, aes(x=ecotype, y=std_length, fill=ecotype))+
  geom_jitter(shape=21, size=2.5, alpha=0.75, width=0.15)+
  geom_point(
    data = eco_means,
    aes(x = ecotype, y = mean_length),
    shape = 3, size = 3, stroke = 1,
    show.legend = FALSE) +
  geom_line(
    data = eco_means,
    aes(x = ecotype, y = mean_length, group = interaction(lake, year)),
    linewidth = 0.5,
    color = "black",
    linetype = "dashed") +
  scale_fill_manual(values = c("#74C476","#70B9A6","#6BAED6"), labels=c("Benthic","Hybrid","Limnetic"))+
  scale_x_discrete(labels=c("benthic"="Benthic","limnetic"="Limnetic","hybrid"="Hybrid"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=11),
        axis.title = element_blank(),
        axis.text.y = element_text(size=8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size=10, face="bold"),
        strip.background = element_rect(fill="grey95"))+
  guides(
    fill = guide_legend(
      override.aes = list(
        size = 3.5,
        alpha = 1,
        shape = 21)))+
  facet_grid2(rows=vars(year),
             cols=vars(lake),
             strip = strip_themed(
               background_x = elem_list_rect(fill = "#70B9A6")))

selection_plot<-plot_grid(source_plot, eco_plot, ncol=2, rel_widths = c(3,1))

ggsave(selection_plot,
       filename = "03_figures/assortment/selection-fig.tiff",
       device = "tiff",
       width = 1200, height=400, units = "px",
       dpi=300, scale=3)

## models for each lake
results_df<-data.frame()
for (lake_name in unique(data$lake)) {
  lake_data <- filter(data, lake == lake_name)

  has_year <- n_distinct(lake_data$year) > 1
  
  model <- if (has_year) {
    lm(std_length ~ hybrid + year, data = lake_data)
  } else {
    lm(std_length ~ hybrid, data = lake_data)
  }
  
  tidy_res <- broom::tidy(model) %>%
    filter(term != "(Intercept)") %>%   # adjust if hybrid has multiple levels
    mutate(
      lake = lake_name
    )
  
  results_df <- bind_rows(results_df, tidy_res)
}

## cleaning
results_clean<-results_df %>% 
  select(lake, term, estimate, std.error, statistic, p.value) %>% 
  mutate(term = case_when(
    term=="hybridPure" ~ "genoytpe",
    term=="yearAge 2" ~ "year",)) %>% 
  group_by(term) %>% 
  mutate(p_adj=p.adjust(p.value, method="BH")) %>% 
  ungroup()

## export
results_clean %>% 
  mutate(across(c(estimate, std.error, statistic), round, 3),
         p.value=if_else(p.value<0.001, "<0.001", as.character(round(p.value, 3))),
         p_adj=if_else(p_adj<0.001, "<0.001", as.character(round(p_adj, 3)))) %>%
  mutate(lake=factor(lake, levels=c(limnetic_lakes, benthic_lakes, both_lakes))) %>% 
  arrange(lake) %>% 
  flextable::flextable() %>% 
  print(preview = "docx")

## LG RL
lg_rl_model<-lm(std_length ~ hybrid + year, data = filter(data, lake %in% "RL", LG_round>0))
summary(lg_rl_model)

## ecotype models
eco_results_df<-data.frame()
for (lake_name in unique(eco_data$lake)) {
  lake_data <- filter(eco_data, lake == lake_name) %>% 
    mutate(ecotype = factor(ecotype),
           ecotype = relevel(ecotype, ref = "hybrid"))
  has_year <- n_distinct(lake_data$year) > 1
  model <- if (has_year) {
    lm(std_length ~ ecotype + year, data = lake_data)
  } else {
    lm(std_length ~ ecotype, data = lake_data)
  }
  tidy_res <- broom::tidy(model) %>%
    filter(term != "(Intercept)") %>%   # adjust if ecotype has multiple levels
    mutate(
      lake = lake_name
    )
  eco_results_df <- bind_rows(eco_results_df, tidy_res)
}

## cleaning
eco_results_clean<-eco_results_df %>%
  #filter(term != "yearAge 2") %>%
  select(lake, term, estimate, std.error, statistic, p.value) %>%
  group_by(term) %>%
  mutate(p_adj=p.adjust(p.value, method="BH")) %>% 
  ungroup()

eco_results_clean %>% 
  mutate(across(c(estimate, std.error, statistic), round, 3),
         p.value=if_else(p.value<0.001, "<0.001", as.character(round(p.value, 3))),
         p_adj=if_else(p_adj<0.001, "<0.001", as.character(round(p_adj, 3)))) %>%
  mutate(lake=factor(lake, levels=both_lakes)) %>% 
  arrange(lake) %>% 
  flextable::flextable() %>% 
  print(preview = "docx")


eco_results_df2<-data.frame()
for (lake_name in unique(eco_data$lake)) {
  lake_data <- filter(eco_data, lake == lake_name)
  has_year <- n_distinct(lake_data$year) > 1
  model <- if (has_year) {
    lm(std_length ~ ecotype + year, data = lake_data)
  } else {
    lm(std_length ~ ecotype, data = lake_data)
  }
  tidy_res <- broom::tidy(model) %>%
    filter(term != "(Intercept)") %>%   # adjust if ecotype has multiple levels
    mutate(
      lake = lake_name
    )
  eco_results_df2 <- bind_rows(eco_results_df2, tidy_res)
}

eco_results2_clean<-eco_results_df2 %>%
  filter(term != "yearAge 2") %>%
  select(lake, term, estimate, std.error, statistic, p.value) %>%
  group_by(term) %>%
  mutate(p_adj=p.adjust(p.value, method="BH")) %>% 
  ungroup()

#### CONDITION ####

## plotting
condition_plot<-ggplot(data, aes(x=hybrid, y=k, fill=hybrid))+
  geom_jitter(shape=21, size=2, alpha=0.75, width=0.15)+
  geom_point(
    data = means,
    aes(x = hybrid, y = mean_k),
    shape = 3, size = 3, stroke = 1,
    show.legend = FALSE) +
  geom_line(
    data = means,
    aes(x = hybrid, y = mean_k, group = interaction(lake, year)),
    linewidth = 0.5,
    color = "black",
    linetype = "dashed") +
  scale_fill_manual(values = c("grey55","grey90"), labels=geno_labs)+
  labs(y="Condition (K)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=11),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size=10, face="bold"),
        strip.background = element_rect(fill="grey95"))+
  guides(
    fill = guide_legend(
      override.aes = list(
        size = 3.5,
        alpha = 1,
        shape = 21)))+
  facet_grid2(
              cols=vars(lake),
              strip = strip_themed(
                background_x = elem_list_rect(fill = all_lake_cols[-1])))

## saving
save_fig(
  plot = condition_plot,
  filename = "03_figures/assortment/condition-fig.tiff",
  width = 1, height=0.3)

## models for each lake
condition_df<-data.frame()
for (lake_name in lakes[-8]) {
  lake_data <- filter(data, lake == lake_name)
  model <- lm(scale(k) ~ hybrid, data = lake_data)
  tidy_res <- broom::tidy(model) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      lake = lake_name
    )
  
  condition_df <- bind_rows(condition_df, tidy_res)
}

## cleaning
condition_clean<-condition_df %>% 
  select(lake, term, estimate, std.error, statistic, p.value) %>% 
  mutate(term = case_when(
    term=="hybridPure" ~ "genoytpe",
    term=="yearAge 2" ~ "year",)) %>% 
  group_by(term) %>% 
  mutate(p_adj=p.adjust(p.value, method="BH")) %>% 
  ungroup()

condition_clean %>% 
  mutate(
    across(c(estimate, std.error, statistic), round, digits=3),
         p.value=if_else(p.value<0.001, "<0.001", as.character(round(p.value, 3))),
         p_adj=if_else(p_adj<0.001, "<0.001", as.character(round(p_adj, 3)))) %>%
  mutate(lake=factor(lake, levels=c(limnetic_lakes, benthic_lakes, both_lakes))) %>% 
  arrange(lake) %>% 
  flextable::flextable() %>% 
  print(preview = "docx")

## LG
lm(scale(k) ~ hybrid, data = filter(data, lake %in% "RL", LG_round>0)) %>% summary()

## eco model
eco_cond_data<-filter(eco_data, lake=="LO") %>% 
  mutate(ecotype = factor(ecotype),
         ecotype = relevel(ecotype, ref = "hybrid"))

lm(k ~ ecotype, data = eco_cond_data) %>% summary()